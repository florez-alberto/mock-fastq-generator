"""Tests for paired-end sequence generation."""

from __future__ import annotations

from mock_fastq_generator.core import assemble_paired_sequences, reverse_complement

TEMPLATE = "GAAAAAATTAACCAGAGCCTG"


class TestPairedEnd:
    """Verifies R1/R2 generation, reverse complementation, and headers."""

    def test_number_of_records(self) -> None:
        """Output should be a tuple of two equal-length lists (R1 and R2)."""
        r1, r2 = assemble_paired_sequences(
            template_sequence=TEMPLATE, total_length=150, number_of_sequences=5
        )
        assert len(r1) == 5
        assert len(r2) == 5

    def test_headers_suffix(self) -> None:
        """R1 headers must end with /1, R2 headers with /2."""
        r1, r2 = assemble_paired_sequences(
            template_sequence=TEMPLATE, total_length=150, number_of_sequences=5
        )
        for rec1, rec2 in zip(r1, r2):
            h1 = rec1[0]
            h2 = rec2[0]
            assert h1.endswith("/1"), f"R1 header doesn't end with /1: {h1}"
            assert h2.endswith("/2"), f"R2 header doesn't end with /2: {h2}"
            # Base headers must match
            assert h1[:-2] == h2[:-2], f"Headers mismatch: {h1} vs {h2}"

    def test_r2_is_reverse_complement_of_r1(self) -> None:
        """R2 sequence must be the exact reverse complement of R1 sequence."""
        # We test with zero noise to ensure exact match.
        r1, r2 = assemble_paired_sequences(
            template_sequence=TEMPLATE,
            total_length=150,
            number_of_sequences=10,
            noise_level=0.0,
        )
        for rec1, rec2 in zip(r1, r2):
            seq1 = rec1[1]
            seq2 = rec2[1]
            assert seq2 == reverse_complement(seq1), (
                f"R2 is not reverse complement of R1.\nR1: {seq1}\nR2: {seq2}"
            )

    def test_quality_scores_are_independent(self) -> None:
        """R1 and R2 should have independently generated noise in quality scores."""
        # High noise ensures they won't match randomly
        r1, r2 = assemble_paired_sequences(
            template_sequence=TEMPLATE,
            total_length=100,
            number_of_sequences=5,
            noise_level=0.5,
        )
        for rec1, rec2 in zip(r1, r2):
            qual1 = rec1[3]
            qual2 = rec2[3]
            # While the base curve is the same, noise means they should differ
            assert qual1 != qual2, "R1 and R2 quality strings are identical"


class TestReverseComplement:
    """Verifies the reverse_complement utility function."""

    def test_basic_rc(self) -> None:
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("A") == "T"
        assert reverse_complement("G") == "C"

    def test_case_preservation(self) -> None:
        """It should preserve the casing of the input."""
        assert reverse_complement("atcg") == "cgat"
        assert reverse_complement("aTcG") == "CgAt"

    def test_palindromes(self) -> None:
        assert reverse_complement("GAATTC") == "GAATTC"  # EcoRI site

    def test_empty(self) -> None:
        assert reverse_complement("") == ""
