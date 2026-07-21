"""Tests that Phred quality scores stay strictly within [min_val, max_val]."""

from __future__ import annotations

import pytest

from mock_fastq_generator.core import generate_phred_scores

ASCII_BASE = 33


def _decode_phred(quality_string: str, ascii_base: int = ASCII_BASE) -> list[int]:
    """Convert a Phred ASCII string back to integer scores."""
    return [ord(ch) - ascii_base for ch in quality_string]


class TestPhredBoundaries:
    """Phred scores must be clamped to [min_val, max_val]."""

    @pytest.mark.parametrize(
        "min_val, max_val",
        [
            (40, 73),
            (33, 93),
            (20, 40),
            (50, 50),  # edge case: min == max
            (60, 73),
        ],
    )
    def test_scores_within_bounds(self, min_val: int, max_val: int) -> None:
        """Every decoded Phred score must lie in [min_val, max_val]."""
        quality = generate_phred_scores(
            num_points=500,
            center=150,
            min_val=min_val,
            max_val=max_val,
            std_dev=75,
            noise_level=0.1,
        )
        scores = [ord(ch) for ch in quality]
        assert all(
            min_val <= s <= max_val for s in scores
        ), (
            f"Found out-of-range ASCII codes: "
            f"{[s for s in scores if not (min_val <= s <= max_val)]}"
        )

    def test_all_characters_valid_ascii(self) -> None:
        """Every character must map to an ASCII code in [33, 126]."""
        quality = generate_phred_scores(
            num_points=1000,
            center=250,
            min_val=40,
            max_val=73,
            std_dev=100,
            noise_level=0.2,
        )
        for ch in quality:
            code = ord(ch)
            assert 33 <= code <= 126, f"Invalid ASCII code {code} for char '{ch}'"

    def test_quality_string_length_matches_request(self) -> None:
        """The returned string length must equal num_points."""
        for n in [1, 10, 100, 500, 1000]:
            quality = generate_phred_scores(
                num_points=n,
                center=n // 2,
                min_val=40,
                max_val=73,
                std_dev=50,
                noise_level=0.05,
            )
            assert len(quality) == n, f"Expected {n}, got {len(quality)}"

    def test_min_equals_max_produces_uniform(self) -> None:
        """When min_val == max_val, all scores should be identical."""
        quality = generate_phred_scores(
            num_points=200,
            center=100,
            min_val=50,
            max_val=50,
            std_dev=30,
            noise_level=0.5,
        )
        scores = [ord(ch) for ch in quality]
        assert all(s == 50 for s in scores), "All ASCII codes must equal 50"

    def test_noise_level_zero(self) -> None:
        """With zero noise, scores should still respect boundaries."""
        quality = generate_phred_scores(
            num_points=300,
            center=150,
            min_val=40,
            max_val=73,
            std_dev=75,
            noise_level=0.0,
        )
        scores = [ord(ch) for ch in quality]
        assert all(40 <= s <= 73 for s in scores)
