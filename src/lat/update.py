"""Automatic update notices for lat.

This module only *notifies*: it never installs or executes remote content.

* The cache lives under ``$XDG_CACHE_HOME/lat/update.json``.
* A cached "new version" notice is printed at most once per interval.
* The cache is refreshed in a daemon thread with a short timeout, so a slow or
  offline network never delays the command.
* Several sources are tried and the highest version wins.  A configured JSON
  endpoint, a PEP 503 simple index, ``PIP_INDEX_URL`` and the public PyPI JSON
  endpoint are all accepted.
"""

from __future__ import annotations

import json
import os
import re
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote
from urllib.request import Request, urlopen

from . import __version__

PACKAGE = "lammps-analysis-tools"
TOOL = "lat"
DEFAULT_SOURCE = f"https://pypi.org/pypi/{PACKAGE}/json"
DEFAULT_INTERVAL_HOURS = 24.0
DEFAULT_TIMEOUT_S = 1.5
MAX_BYTES = 1024 * 1024
_TRUE = {"1", "true", "yes", "on"}
_FALSE = {"0", "false", "no", "off"}
_THREAD_LOCK = threading.Lock()
_THREAD_STARTED = False


@dataclass(slots=True)
class CheckResult:
    latest: str | None
    source: str | None
    errors: list[str]


def _truthy(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in _TRUE


def _cfg_get(cfg, key, default=None):
    if cfg is None:
        return default
    try:
        value = cfg.get(key, default=default)
    except TypeError:
        value = cfg.get(key)
    except Exception:
        return default
    return default if value is None else value


def _version_key(value):
    try:
        from packaging.version import Version

        return Version(str(value))
    except Exception:
        core = str(value).split("+", 1)[0]
        match = re.match(r"^(\d+(?:\.\d+)*)", core)
        numbers = tuple(int(part) for part in match.group(1).split(".")) if match else (0,)
        pre = -1 if re.search(r"(?:a|b|rc|dev)", core) else 0
        return numbers + (pre,)


def is_newer(latest, current) -> bool:
    try:
        return _version_key(latest) > _version_key(current)
    except Exception:
        return False


def _is_prerelease(value) -> bool:
    try:
        from packaging.version import Version

        return Version(str(value)).is_prerelease
    except Exception:
        return bool(re.search(r"(?:a|b|rc|dev)", str(value)))


def cache_path() -> Path:
    base = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
    return Path(base) / TOOL / "update.json"


def enabled(cfg=None) -> bool:
    if _truthy(os.environ.get(f"{TOOL.upper()}_NO_UPDATE_CHECK")):
        return False
    if not _truthy(os.environ.get(f"{TOOL.upper()}_UPDATE_CHECK"), default=True):
        return False
    return _truthy(_cfg_get(cfg, "update.enabled", True), default=True)


def interval_seconds(cfg=None) -> float:
    raw = os.environ.get(f"{TOOL.upper()}_UPDATE_INTERVAL_HOURS")
    if raw is None:
        raw = _cfg_get(cfg, "update.interval_hours", DEFAULT_INTERVAL_HOURS)
    try:
        return max(float(raw), 0.0) * 3600.0
    except (TypeError, ValueError):
        return DEFAULT_INTERVAL_HOURS * 3600.0


def timeout_seconds(cfg=None) -> float:
    raw = os.environ.get(f"{TOOL.upper()}_UPDATE_TIMEOUT_S")
    if raw is None:
        raw = _cfg_get(cfg, "update.timeout_s", DEFAULT_TIMEOUT_S)
    try:
        return max(float(raw), 0.1)
    except (TypeError, ValueError):
        return DEFAULT_TIMEOUT_S


def _as_project_url(base: str) -> str:
    base = base.rstrip("/")
    if base.endswith("/" + PACKAGE) or base.endswith("/" + PACKAGE.replace("-", "_")):
        return base + "/"
    return base + "/" + PACKAGE + "/"


def sources(cfg=None) -> list[str]:
    out: list[str] = []
    raw_json = os.environ.get(f"{TOOL.upper()}_UPDATE_URL")
    if raw_json:
        out.append(raw_json)
    raw_index = os.environ.get(f"{TOOL.upper()}_UPDATE_INDEX")
    if raw_index:
        out.append(_as_project_url(raw_index))
    pip_index = os.environ.get("PIP_INDEX_URL")
    if pip_index:
        out.append(_as_project_url(pip_index))
    configured = _cfg_get(cfg, "update.sources", None)
    if isinstance(configured, (list, tuple)):
        out.extend(str(item) for item in configured)
    elif isinstance(configured, str) and configured.strip():
        out.append(configured.strip())
    single = _cfg_get(cfg, "update.url", None)
    if single:
        out.append(str(single))
    out.append(DEFAULT_SOURCE)
    seen: set[str] = set()
    unique: list[str] = []
    for item in out:
        item = item.strip()
        if item and item not in seen:
            seen.add(item)
            unique.append(item)
    return unique


def _fetch(url: str, timeout: float) -> str:
    request = Request(url, headers={"User-Agent": f"{TOOL}/{__version__}"})
    with urlopen(request, timeout=timeout) as response:
        payload = response.read(MAX_BYTES + 1)
    if len(payload) > MAX_BYTES:
        raise ValueError(f"response from {url} exceeds {MAX_BYTES} bytes")
    return payload.decode("utf-8", "replace")


def _version_from_filename(filename: str) -> str | None:
    name = unquote(filename.rsplit("/", 1)[-1])
    for prefix in (PACKAGE.replace("-", "_") + "-", PACKAGE + "-"):
        if not name.startswith(prefix):
            continue
        rest = name[len(prefix):]
        if rest.endswith(".tar.gz"):
            return rest[:-7]
        if rest.endswith(".zip"):
            return rest[:-4]
        return rest.split("-", 1)[0]
    return None


def _latest_from_text(text: str) -> str | None:
    stripped = text.strip()
    try:
        payload = json.loads(stripped)
    except Exception:
        payload = None
    if isinstance(payload, dict):
        info = payload.get("info")
        if isinstance(info, dict) and info.get("version"):
            return str(info["version"])
        releases = payload.get("releases")
        if isinstance(releases, dict):
            versions = [str(key) for key, files in releases.items() if files]
            if versions:
                return max(versions, key=_version_key)
    if re.fullmatch(r"v?\d+(?:\.\d+)*(?:[a-z0-9.+-]*)?", stripped):
        return stripped.lstrip("v")
    versions = []
    for match in re.finditer(r'href="([^"]+)"', text):
        version = _version_from_filename(match.group(1))
        if version:
            versions.append(version)
    for match in re.finditer(r">([^<>]+\.(?:whl|tar\.gz|zip))<", text):
        version = _version_from_filename(match.group(1).strip())
        if version:
            versions.append(version)
    if not versions:
        return None
    return max(versions, key=_version_key)


def check_sources(source_list=None, *, cfg=None, timeout=None, opener=None) -> CheckResult:
    source_list = list(source_list or sources(cfg))
    timeout = timeout if timeout is not None else timeout_seconds(cfg)
    opener = opener or _fetch
    latest = None
    latest_source = None
    errors: list[str] = []
    for url in source_list:
        try:
            version = _latest_from_text(opener(url, timeout))
        except Exception as exc:
            errors.append(f"{url}: {exc}")
            continue
        if not version:
            errors.append(f"{url}: no version found")
            continue
        if latest is None or is_newer(version, latest):
            latest, latest_source = version, url
    return CheckResult(latest=latest, source=latest_source, errors=errors)


def _read_cache(path: Path | None = None) -> dict:
    path = Path(path) if path is not None else cache_path()
    try:
        data = json.loads(path.read_text())
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _write_cache(data: dict, path: Path | None = None) -> None:
    path = Path(path) if path is not None else cache_path()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
        tmp.write_text(json.dumps(data, ensure_ascii=False, sort_keys=True))
        os.replace(tmp, path)
    except Exception:
        pass


def _refresh(cfg=None, *, now=None, cache_file=None, opener=None) -> CheckResult:
    now = time.time() if now is None else now
    result = check_sources(cfg=cfg, opener=opener)
    latest = result.latest
    if latest is None:
        latest = __version__
    if _is_prerelease(latest) and not _is_prerelease(__version__):
        latest = __version__
    data = {
        "checked_at": now,
        "latest": latest,
        "source": result.source,
        "error": result.errors[-1] if result.errors else None,
    }
    _write_cache(data, cache_file)
    return result


def _start_refresh(cfg=None, *, cache_file=None) -> None:
    global _THREAD_STARTED
    with _THREAD_LOCK:
        if _THREAD_STARTED:
            return
        _THREAD_STARTED = True

    def worker():
        try:
            _refresh(cfg, cache_file=cache_file)
        except BaseException:
            pass

    threading.Thread(target=worker, name=f"{TOOL}-update-check", daemon=True).start()


def _notice(latest: str) -> str:
    return (f"[{TOOL} NOTICE] new version {latest} available (current {__version__}); "
            f"upgrade with: pip install -U {PACKAGE}")


def maybe_notify(cfg=None, *, interactive=None, stream=None, now=None,
                 cache_file=None, start_thread=True) -> str | None:
    """Print a cached update notice and refresh the cache in the background."""
    if not enabled(cfg):
        return None
    if interactive is None:
        interactive = bool(getattr(stream or sys.stderr, "isatty", lambda: False)())
        if os.environ.get("CI"):
            interactive = False
    if not interactive:
        return None
    now = time.time() if now is None else now
    cache = _read_cache(cache_file)
    latest = cache.get("latest")
    interval = interval_seconds(cfg)
    if latest and is_newer(latest, __version__) and not (
            _is_prerelease(latest) and not _is_prerelease(__version__)):
        notified = float(cache.get("notified_at", 0) or 0)
        if now - notified >= interval:
            print(_notice(str(latest)), file=stream or sys.stderr)
            cache["notified_at"] = now
            _write_cache(cache, cache_file)
        return str(latest)
    checked = float(cache.get("checked_at", 0) or 0)
    if start_thread and now - checked >= interval:
        _start_refresh(cfg, cache_file=cache_file)
    return None


def force_check(cfg=None, *, stream=None, timeout=None, opener=None,
                cache_file=None) -> int:
    """Synchronous check used by the explicit ``--check-update`` command."""
    stream = stream or sys.stderr
    timeout = max(timeout if timeout is not None else timeout_seconds(cfg), 3.0)
    result = check_sources(cfg=cfg, timeout=timeout, opener=opener)
    if result.latest is None:
        detail = f" ({result.errors[-1]})" if result.errors else ""
        print(f"[{TOOL} NOTICE] update check failed{detail}", file=stream)
        return 1
    now = time.time()
    _write_cache({"checked_at": now, "latest": result.latest,
                  "source": result.source, "error": None}, cache_file)
    if is_newer(result.latest, __version__) and not (
            _is_prerelease(result.latest) and not _is_prerelease(__version__)):
        print(_notice(result.latest), file=stream)
    else:
        print(f"[{TOOL}] {__version__} is up to date (source: {result.source})", file=stream)
    return 0
