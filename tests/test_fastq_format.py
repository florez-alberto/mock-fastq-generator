"""Tests verifying the standard 4-line FASTQ format and contents."""

from __future__ import annotations

import re

from mock_fastq_generator.core import assemble_sequences

TEMPLATE = "GAAAAAATTAACCAGAGCCTG"


class TestFastqFormat:
    """Verifies that generated records match the Sanger FASTQ format."""

    def test_record_has_four_elements(self) -> None:
        """Every FASTQ record must be a list of exactly 4 strings."""
        records = assemble_sequences(
            template_sequence=TEMPLATE, total_length=100, number_of_sequences=5
        )
        for record in records:
            assert isinstance(record, list)
            assert len(record) == 4
            assert all(isinstance(elem, str) for elem in record)

    def test_line1_header_format(self) -> None:
        """Line 1 must start with '@'."""
        records = assemble_sequences(
            template_sequence=TEMPLATE, total_length=100, number_of_sequences=10
        )
        for record in records:
            header = record[0]
            assert header.startswith("@"), f"Header missing '@': {header}"
            # Expecting @XXXX:1234
            assert re.match(r"^@[A-Z0-9]{4}:\d{4}$", header), (
                f"Invalid header format: {header}"
            )

    def test_line2_dna_alphabet(self) -> None:
        """Line 2 must contain only A, C, G, T characters."""
        records = assemble_sequences(
            template_sequence=TEMPLATE, total_length=150, number_of_sequences=10
        )
        for record in records:
            seq = record[1]
            assert set(seq).issubset({"A", "C", "G", "T"}), (
                f"Invalid characters in seq: {seq}"
            )

    def test_line3_separator(self) -> None:
        """Line 3 must be precisely '+'."""
        records = assemble_sequences(
            template_sequence=TEMPLATE, total_length=100, number_of_sequences=5
        )
        for record in records:
            assert record[2] == "+", f"Invalid separator: {record[2]}"

    def test_line4_quality_ascii(self) -> None:
        """Line 4 must be ASCII characters matching seq length."""
        records = assemble_sequences(
            template_sequence=TEMPLATE, total_length=100, number_of_sequences=5
        )
        for record in records:
            seq = record[1]
            qual = record[3]
            assert len(qual) == len(seq)
            for ch in qual:
                # 33 to 126 is the valid printable ASCII range for Phred (Sanger)
                assert 33 <= ord(ch) <= 126, f"Invalid Phred char {ch} (ord={ord(ch)})"
