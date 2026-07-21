"""Tests that generated sequence lengths match the requested total_length."""

from __future__ import annotations

import pytest

from mock_fastq_generator.core import assemble_sequences

# A short template that fits comfortably in the default total_length.
TEMPLATE = "GAAAAAATTAACCAGAGCCTG"


class TestSequenceLengths:
    """Every generated record must have sequence length == total_length,
    and the quality string must have the same length as the sequence."""

    @pytest.mark.parametrize("total_length", [100, 250, 500, 1000])
    def test_sequence_length_matches_total_length(self, total_length: int) -> None:
        """The DNA sequence in every record must have len == total_length."""
        records = assemble_sequences(
            template_sequence=TEMPLATE,
            left_margin=10,
            upstream_sequence="GCCGGCCATGGCG",
            total_length=total_length,
            number_of_sequences=5,
            center=total_length // 2,
            min_val=40,
            max_val=73,
            std_dev=50,
            noise_level=0.1,
        )
        for i, record in enumerate(records):
            seq = record[1]
            assert len(seq) == total_length, (
                f"Record {i}: expected len={total_length}, got {len(seq)}"
            )

    @pytest.mark.parametrize("total_length", [100, 250, 500, 1000])
    def test_quality_length_matches_sequence_length(self, total_length: int) -> None:
        """The quality string must be the same length as the DNA sequence."""
        records = assemble_sequences(
            template_sequence=TEMPLATE,
            left_margin=10,
            upstream_sequence="GCCGGCCATGGCG",
            total_length=total_length,
            number_of_sequences=5,
            center=total_length // 2,
            min_val=40,
            max_val=73,
            std_dev=50,
            noise_level=0.1,
        )
        for i, record in enumerate(records):
            seq, qual = record[1], record[3]
            assert len(qual) == len(seq), (
                f"Record {i}: seq len={len(seq)}, qual len={len(qual)}"
            )

    def test_all_records_same_length(self) -> None:
        """Within a single generation call, all records have equal length."""
        records = assemble_sequences(
            template_sequence=TEMPLATE,
            left_margin=15,
            upstream_sequence="ATCG",
            total_length=300,
            number_of_sequences=20,
            center=100,
            min_val=40,
            max_val=73,
            std_dev=50,
            noise_level=0.1,
        )
        lengths = {len(r[1]) for r in records}
        assert len(lengths) == 1, f"Expected uniform length, got {lengths}"

    def test_number_of_records(self) -> None:
        """The number of output records must match number_of_sequences."""
        for n in [1, 5, 50]:
            records = assemble_sequences(
                template_sequence=TEMPLATE,
                left_margin=10,
                upstream_sequence="GCCGGCCATGGCG",
                total_length=200,
                number_of_sequences=n,
                center=100,
                min_val=40,
                max_val=73,
                std_dev=50,
                noise_level=0.1,
            )
            assert len(records) == n, f"Expected {n} records, got {len(records)}"
