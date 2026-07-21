"""Tests for reading FASTA and writing FASTQ formats."""

from __future__ import annotations

from io import StringIO
from pathlib import Path
from unittest.mock import patch

import pytest

from mock_fastq_generator.io import read_fasta, write_fastq


class TestReadFasta:
    """Verifies robust FASTA parsing capabilities."""

    def test_single_line_fasta(self, tmp_path: Path) -> None:
        """It should parse a simple FASTA with one sequence line."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">Header\nATCGATCG\n")
        assert read_fasta(str(fasta_file)) == "ATCGATCG"

    def test_multi_line_fasta(self, tmp_path: Path) -> None:
        """It should concatenate multi-line sequences."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">Header\nATCG\nATCG\n")
        assert read_fasta(str(fasta_file)) == "ATCGATCG"

    def test_handles_blank_lines(self, tmp_path: Path) -> None:
        """It should skip empty lines."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">Header\n\nATCG\n\n\nATCG\n")
        assert read_fasta(str(fasta_file)) == "ATCGATCG"

    def test_uppercases_output(self, tmp_path: Path) -> None:
        """It should convert lower-case sequences to uppercase."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">Header\natcgATcg\n")
        assert read_fasta(str(fasta_file)) == "ATCGATCG"

    def test_empty_fasta_raises_error(self, tmp_path: Path) -> None:
        """It should raise ValueError if no sequence data is found."""
        fasta_file = tmp_path / "empty.fasta"
        fasta_file.write_text(">Header only\n")
        with pytest.raises(ValueError, match="No sequence data found"):
            read_fasta(str(fasta_file))

    def test_missing_file_raises_error(self) -> None:
        """It should raise FileNotFoundError for missing paths."""
        with pytest.raises(FileNotFoundError):
            read_fasta("does_not_exist.fasta")


class TestWriteFastq:
    """Verifies FASTQ output options."""

    @pytest.fixture
    def sample_records(self) -> list[list[str]]:
        return [
            ["@READ1:1234", "ATCG", "+", "IIII"],
            ["@READ2:5678", "GCTA", "+", "::::"],
        ]

    def test_write_plain_file(
        self, tmp_path: Path, sample_records: list[list[str]]
    ) -> None:
        """It should write plain text records to a file."""
        out_file = tmp_path / "out.fastq"
        write_fastq(str(out_file), sample_records)

        content = out_file.read_text()
        expected = "@READ1:1234\nATCG\n+\nIIII\n@READ2:5678\nGCTA\n+\n::::\n"
        assert content == expected

    def test_write_stdout(self, sample_records: list[list[str]]) -> None:
        """It should write to sys.stdout when use_stdout=True."""
        # We patch sys.stdout with a StringIO object to capture output
        captured = StringIO()
        with patch("sys.stdout", captured):
            write_fastq(None, sample_records, use_stdout=True)

        expected = "@READ1:1234\nATCG\n+\nIIII\n@READ2:5678\nGCTA\n+\n::::\n"
        assert captured.getvalue() == expected

    def test_path_required_if_not_stdout(self, sample_records: list[list[str]]) -> None:
        """It should raise ValueError if path is None and not using stdout."""
        with pytest.raises(ValueError, match="path must be provided"):
            write_fastq(None, sample_records, use_stdout=False)
