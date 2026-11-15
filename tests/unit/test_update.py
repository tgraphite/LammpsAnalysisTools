"""Tests for the opt-in update notice (no network is touched)."""

from __future__ import annotations

import io
import json

import pytest

from lat import __version__, update


def test_version_order():
    assert update.is_newer("1.1.0", "1.0.0")
    assert update.is_newer("1.0.1", "1.0.0")
    assert not update.is_newer("1.0.0", "1.0.0")
    assert not update.is_newer("0.9.9", "1.0.0")
    assert update.is_newer("1.1.0rc1", "1.0.0")
    assert update.is_newer("1.0.0", "1.0.0rc1")


def test_parse_pypi_json():
    payload = json.dumps({"info": {"version": "9.9.9"}, "releases": {}})
    assert update._latest_from_text(payload) == "9.9.9"


def test_parse_simple_index_wheel_and_sdist():
    html = (
        '<a href="../../packages/lammps_analysis_tools-1.0.0-py3-none-any.whl">'
        "lammps_analysis_tools-1.0.0-py3-none-any.whl</a>"
        '<a href="../../packages/lammps_analysis_tools-1.2.0.tar.gz">'
        "lammps_analysis_tools-1.2.0.tar.gz</a>"
    )
    assert update._latest_from_text(html) == "1.2.0"


def test_check_sources_takes_highest():
    def opener(url, timeout):
        if "pypi" in url:
            return json.dumps({"info": {"version": "1.1.0"}})
        return (
            '<a href="x/lammps_analysis_tools-1.2.0-py3-none-any.whl">'
            "lammps_analysis_tools-1.2.0-py3-none-any.whl</a>"
        )

    result = update.check_sources(["https://pypi.test/json", "https://mirror.test/simple/"],
                                  timeout=0.1, opener=opener)
    assert result.latest == "1.2.0"
    assert result.source == "https://mirror.test/simple/"
    assert result.errors == []


def test_check_sources_reports_errors():
    def opener(url, timeout):
        raise OSError("offline")

    result = update.check_sources(["https://example.test"], timeout=0.1, opener=opener)
    assert result.latest is None
    assert result.errors


def test_disabled_env(monkeypatch, tmp_path):
    monkeypatch.setenv("LAT_NO_UPDATE_CHECK", "1")
    stream = io.StringIO()
    assert update.maybe_notify(interactive=True, stream=stream,
                               cache_file=tmp_path / "cache.json") is None
    assert stream.getvalue() == ""


def test_cached_notice_is_printed_once(tmp_path):
    cache = tmp_path / "cache.json"
    now = 1_000_000.0
    cache.write_text(json.dumps({"checked_at": now, "latest": "9.9.9",
                                 "notified_at": 0}))
    stream = io.StringIO()
    result = update.maybe_notify(interactive=True, stream=stream, now=now,
                                 cache_file=cache, start_thread=False)
    assert result == "9.9.9"
    assert "new version 9.9.9" in stream.getvalue()
    first = cache.read_text()
    assert json.loads(first)["notified_at"] == now

    stream2 = io.StringIO()
    assert update.maybe_notify(interactive=True, stream=stream2, now=now + 60,
                               cache_file=cache, start_thread=False) == "9.9.9"
    assert stream2.getvalue() == ""


def test_stale_cache_starts_refresh(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(update, "_start_refresh", lambda cfg=None, cache_file=None: calls.append(1))
    cache = tmp_path / "cache.json"
    cache.write_text(json.dumps({"checked_at": 0, "latest": __version__}))
    update.maybe_notify(interactive=True, stream=io.StringIO(), now=10_000_000,
                        cache_file=cache, start_thread=True)
    assert calls == [1]


def test_non_interactive_does_not_check(monkeypatch, tmp_path):
    called = []
    monkeypatch.setattr(update, "_start_refresh", lambda cfg=None, cache_file=None: called.append(1))
    assert update.maybe_notify(interactive=False, stream=io.StringIO(),
                               cache_file=tmp_path / "cache.json") is None
    assert called == []


def test_force_check_new_version(tmp_path):
    def opener(url, timeout):
        return json.dumps({"info": {"version": "99.0.0"}})

    stream = io.StringIO()
    code = update.force_check(opener=opener, stream=stream, timeout=0.1,
                              cache_file=tmp_path / "cache.json")
    assert code == 0
    assert "new version 99.0.0" in stream.getvalue()


def test_force_check_up_to_date(tmp_path):
    def opener(url, timeout):
        return json.dumps({"info": {"version": __version__}})

    stream = io.StringIO()
    code = update.force_check(opener=opener, stream=stream, timeout=0.1,
                              cache_file=tmp_path / "cache.json")
    assert code == 0
    assert "up to date" in stream.getvalue()


def test_force_check_failure_is_non_fatal(tmp_path):
    def opener(url, timeout):
        raise OSError("offline")

    stream = io.StringIO()
    code = update.force_check(opener=opener, stream=stream, timeout=0.1,
                              cache_file=tmp_path / "cache.json")
    assert code == 1
    assert "update check failed" in stream.getvalue()
