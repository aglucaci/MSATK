from msatk.cli import main


def test_help_command():
    try:
        main(["--help"])
    except SystemExit as exc:
        assert exc.code == 0


def test_profile_help_command():
    try:
        main(["profile", "--help"])
    except SystemExit as exc:
        assert exc.code == 0


def test_demo_command(tmp_path):
    assert main(["demo", "--out", str(tmp_path), "--force"]) == 0
    assert (tmp_path / "report.html").exists()


def test_tables_only(data_dir, tmp_path):
    assert (
        main(
            [
                "profile",
                str(data_dir / "dna" / "tiny_dna.fasta"),
                "--out",
                str(tmp_path),
                "--tables-only",
                "--force",
            ]
        )
        == 0
    )
    assert (tmp_path / "summary.json").exists()
    assert not (tmp_path / "report.html").exists()
