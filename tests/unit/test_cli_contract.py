"""CLI contract tests surfaced by the docs/UX review."""

import subprocess
import sys


def lat(repo_root, args):
    return subprocess.run([sys.executable, str(repo_root / "lat"), *args],
                          capture_output=True, text=True)


def test_stdout_clean_stderr_has_path(repo_root, sample_dir, tmp_path):
    base = sample_dir / "sys1_reduced"
    p = lat(repo_root, ["rdf", "-s", str(base) + ".data", "-f", str(base) + ".lammpstrj",
                        "-sel1", "type 2", "-sel2", "type 5", "-cutoff", "8",
                        "-o", str(tmp_path / "r")])
    assert p.returncode == 0, p.stderr
    assert p.stdout == ""                       # stdout never polluted
    assert "wrote" in p.stderr and ".rdf.out.csv" in p.stderr


def test_warnings_shown_on_screen(repo_root, sample_dir, tmp_path):
    base = sample_dir / "sys1_reduced"
    # sample dump has no image flags -> msd must WARN on stderr, not only in header
    p = lat(repo_root, ["msd", "-s", str(base) + ".data", "-f", str(base) + ".lammpstrj",
                        "-sel", "type 2", "-o", str(tmp_path / "m")])
    assert p.returncode == 0, p.stderr
    assert "image flags" in p.stderr.lower()


def test_csv_header_has_plot_contract(repo_root, sample_dir, tmp_path):
    base = sample_dir / "sys1_reduced"
    p = lat(repo_root, ["rdf", "-s", str(base) + ".data", "-f", str(base) + ".lammpstrj",
                        "-sel1", "type 2", "-sel2", "type 5", "-cutoff", "8",
                        "-o", str(tmp_path / "r")])
    assert p.returncode == 0, p.stderr
    header = (tmp_path / "r.rdf.out.csv").read_text()
    for field in ("# program: LAT version", "# module: rdf", "# title: RDF",
                  "# xlabel: r (Angstrom)", "# ylabel: g(r)",
                  "# selections: sel1=type 2", "# definitions: cutoff=8 A",
                  "# columns: r, g_"):
        assert field in header, f"missing header field: {field}"


def test_plot1d_uses_dash_f_only(repo_root, sample_dir, tmp_path):
    base = sample_dir / "sys1_reduced"
    lat(repo_root, ["rdf", "-s", str(base) + ".data", "-f", str(base) + ".lammpstrj",
                    "-sel1", "type 2", "-sel2", "type 5", "-cutoff", "8",
                    "-o", str(tmp_path / "r")])
    csv = tmp_path / "r.rdf.out.csv"
    ok = lat(repo_root, ["plot1d", "-f", str(csv), "-o", str(tmp_path / "fig.png")])
    assert ok.returncode == 0, ok.stderr
    assert (tmp_path / "fig.png").is_file()
    # the old positional form must no longer be accepted
    bad = lat(repo_root, ["plot1d", str(csv)])
    assert bad.returncode != 0


def test_plot1d_smooth_is_a_plot_time_option(repo_root, sample_dir, tmp_path):
    base = sample_dir / "sys1_reduced"
    lat(repo_root, ["rdf", "-s", str(base) + ".data", "-f", str(base) + ".lammpstrj",
                    "-sel1", "type 2", "-sel2", "type 5", "-cutoff", "8",
                    "-o", str(tmp_path / "r")])
    csv = tmp_path / "r.rdf.out.csv"
    for method in ("none", "gaussian", "savgol"):
        p = lat(repo_root, ["plot1d", "-f", str(csv), "-smooth", method,
                            "-o", str(tmp_path / f"{method}.png")])
        assert p.returncode == 0, p.stderr
        assert (tmp_path / f"{method}.png").is_file()
    bad = lat(repo_root, ["plot1d", "-f", str(csv), "-smooth", "boxcar"])
    assert bad.returncode != 0


def _write_lat_csv(path):
    from lat import __version__

    path.write_text(
        f"# program: LAT version {__version__}\n"
        "# module: test\n"
        "# columns: frame, A, B, C\n"
        "# --- data ---\n"
        "0,1,2,3\n"
        "1,2,3,4\n"
        "2,3,4,5\n"
    )


def test_plot1d_columns_selects_y_series(repo_root, tmp_path):
    csv = tmp_path / "d.out.csv"
    _write_lat_csv(csv)
    ok = lat(repo_root, ["plot1d", "-f", str(csv), "-columns", "C,A",
                         "-o", str(tmp_path / "sub.png")])
    assert ok.returncode == 0, ok.stderr
    assert (tmp_path / "sub.png").is_file()
    # unknown name and x-only selections must fail
    assert lat(repo_root, ["plot1d", "-f", str(csv), "-columns", "Z"]).returncode != 0
    assert lat(repo_root, ["plot1d", "-f", str(csv), "-columns", "frame"]).returncode != 0


def test_plot1d_bar_and_stack(repo_root, tmp_path):
    csv = tmp_path / "d.out.csv"
    _write_lat_csv(csv)
    for flag, name in (("-bar", "bar.png"), ("-stack", "stack.png")):
        p = lat(repo_root, ["plot1d", "-f", str(csv), flag,
                            "-o", str(tmp_path / name)])
        assert p.returncode == 0, p.stderr
        assert (tmp_path / name).is_file()
    # mutually exclusive, and smoothing is meaningless for bars
    both = lat(repo_root, ["plot1d", "-f", str(csv), "-bar", "-stack"])
    assert both.returncode != 0 and "mutually exclusive" in both.stderr
    sm = lat(repo_root, ["plot1d", "-f", str(csv), "-stack", "-smooth", "gaussian"])
    assert sm.returncode != 0
    # plot2d does not accept the new options
    p2 = lat(repo_root, ["plot2d", "-f", str(csv), "-bar"])
    assert p2.returncode != 0


def test_help_is_self_documenting(repo_root):
    p = lat(repo_root, ["density", "-h"])
    assert p.returncode == 0
    assert "{x,y,z}" in p.stdout and "[default:" in p.stdout
    assert "describe -module density" in p.stdout
