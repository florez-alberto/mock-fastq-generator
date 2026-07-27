"""Unit tests for the command-line interface (cli.py)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from mock_fastq_generator.cli import _resolve_parameters, build_parser, main


class TestCLI:
    """Test CLI argument parsing and main entry point execution."""

    def test_build_parser_defaults(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["--output_file", "test.fastq"])
        assert args.output_file == "test.fastq"
        assert args.total_length == 500
        assert args.number_of_sequences == 50
        assert args.decay_model == "gaussian"

    def test_resolve_parameters_from_args(self) -> None:
        parser = build_parser()
        args = parser.parse_args(
            [
                "--output_file",
                "test.fastq",
                "--number_of_sequences",
                "10",
                "--decay_model",
                "sigmoidal",
                "--decay_rate",
                "1.0",
            ]
        )
        params = _resolve_parameters(args)
        assert params["number_of_sequences"] == 10
        assert params["decay_model"] == "sigmoidal"
        assert params["decay_rate"] == 1.0

    def test_resolve_parameters_from_json_file(self, tmp_path: Path) -> None:
        json_file = tmp_path / "config.json"
        json_file.write_text(
            json.dumps(
                {
                    "number_of_sequences": 25,
                    "decay_model": "exponential",
                    "min_val": 10,
                    "max_val": 35,
                }
            )
        )

        parser = build_parser()
        args = parser.parse_args(
            [
                "--output_file",
                "test.fastq",
                "--parameters_file",
                str(json_file),
            ]
        )
        params = _resolve_parameters(args)
        assert params["number_of_sequences"] == 25
        assert params["decay_model"] == "exponential"
        assert params["min_val"] == 10
        assert params["max_val"] == 35

    def test_main_single_end(self, tmp_path: Path) -> None:
        out_file = tmp_path / "single.fastq"
        parser = build_parser()
        args = parser.parse_args(
            [
                "--output_file",
                str(out_file),
                "--number_of_sequences",
                "5",
                "--total_length",
                "100",
            ]
        )
        main(args)

        assert out_file.exists()
        lines = out_file.read_text().strip().split("\n")
        assert len(lines) == 20  # 5 records * 4 lines

    def test_main_paired_end(self, tmp_path: Path) -> None:
        out_prefix = tmp_path / "paired"
        parser = build_parser()
        args = parser.parse_args(
            [
                "--output_file",
                str(out_prefix),
                "--paired-end",
                "--number_of_sequences",
                "5",
                "--total_length",
                "100",
            ]
        )
        main(args)

        r1_path = tmp_path / "paired_R1.fastq"
        r2_path = tmp_path / "paired_R2.fastq"
        assert r1_path.exists()
        assert r2_path.exists()

    def test_main_gzip_output(self, tmp_path: Path) -> None:
        out_file = tmp_path / "single.fastq.gz"
        parser = build_parser()
        args = parser.parse_args(
            [
                "--output_file",
                str(out_file),
                "--gzip",
                "--number_of_sequences",
                "5",
                "--total_length",
                "100",
            ]
        )
        main(args)
        assert out_file.exists()

    def test_main_stdout(self, capsys: pytest.CaptureFixture[str]) -> None:
        parser = build_parser()
        args = parser.parse_args(
            [
                "--stdout",
                "--number_of_sequences",
                "2",
                "--total_length",
                "50",
            ]
        )
        main(args)
        captured = capsys.readouterr()
        lines = captured.out.strip().split("\n")
        assert len(lines) == 8  # 2 records * 4 lines

    def test_main_paired_end_stdout(self, capsys: pytest.CaptureFixture[str]) -> None:
        parser = build_parser()
        args = parser.parse_args(
            [
                "--stdout",
                "--paired-end",
                "--number_of_sequences",
                "2",
                "--total_length",
                "50",
            ]
        )
        main(args)
        captured = capsys.readouterr()
        lines = captured.out.strip().split("\n")
        assert len(lines) == 16  # 2 pairs * 2 * 4 lines

    def test_main_default_template(self, tmp_path: Path) -> None:
        out_file = tmp_path / "default_tmpl.fastq"
        parser = build_parser()
        args = parser.parse_args(
            [
                "--output_file",
                str(out_file),
                "--number_of_sequences",
                "2",
                "--total_length",
                "50",
            ]
        )
        main(args)
        assert out_file.exists()

    def test_main_paired_end_suffix_stripping(self, tmp_path: Path) -> None:
        out_file = tmp_path / "sample.fq"
        parser = build_parser()
        args = parser.parse_args(
            [
                "--output_file",
                str(out_file),
                "--paired-end",
                "--number_of_sequences",
                "2",
                "--total_length",
                "50",
            ]
        )
        main(args)
        assert (tmp_path / "sample_R1.fastq").exists()
        assert (tmp_path / "sample_R2.fastq").exists()

    def test_main_single_gzip_append_extension(self, tmp_path: Path) -> None:
        out_file = tmp_path / "sample.fastq"
        parser = build_parser()
        args = parser.parse_args(
            [
                "--output_file",
                str(out_file),
                "--gzip",
                "--number_of_sequences",
                "2",
                "--total_length",
                "50",
            ]
        )
        main(args)
        assert (tmp_path / "sample.fastq.gz").exists()

    def test_main_sys_argv(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        out_file = tmp_path / "argv_test.fastq"
        monkeypatch.setattr(
            "sys.argv",
            [
                "mock-fastq-generator",
                "--output_file",
                str(out_file),
                "--number_of_sequences",
                "1",
            ],
        )
        main(None)
        assert out_file.exists()
