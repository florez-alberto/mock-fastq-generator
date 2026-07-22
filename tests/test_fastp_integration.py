"""Integration tests verifying generated FASTQ compatibility with third-party tools."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from mock_fastq_generator.core import assemble_sequences
from mock_fastq_generator.io import write_fastq

# Skip the entire module if fastp is not installed
pytestmark = pytest.mark.skipif(
    shutil.which("fastp") is None, reason="fastp is not installed on the system PATH"
)

TEMPLATE = "GAAAAAATTAACCAGAGCCTGG"


class TestFastpIntegration:
    """Verifies that fastp correctly parses our files and catches low quality reads."""

    def test_fastp_low_quality_filtering(self, tmp_path: Path) -> None:
        """fastp should filter out reads intentionally generated with low quality."""
        # 1. Generate FASTQ with intentionally low quality scores
        # min_val=33 is Phred 0, max_val=40 is Phred 7.
        # fastp by default requires >= Phred 15 (Q15).
        records = assemble_sequences(
            template_sequence=TEMPLATE,
            total_length=150,
            number_of_sequences=100,
            min_val=33,
            max_val=40,
        )

        fastq_path = tmp_path / "low_qual.fastq"
        write_fastq(str(fastq_path), records)

        json_report_path = tmp_path / "fastp.json"
        html_report_path = tmp_path / "fastp.html"

        # 2. Run fastp via subprocess
        cmd = [
            "fastp",
            "-i",
            str(fastq_path),
            "-j",
            str(json_report_path),
            "-h",
            str(html_report_path),
        ]

        # Run fastp and capture output (fastp writes logs to stderr)
        result = subprocess.run(cmd, capture_output=True, text=True)

        assert result.returncode == 0, f"fastp failed: {result.stderr}"
        assert json_report_path.exists(), "fastp did not generate the JSON report"

        # 3. Verify fastp caught the low quality
        report = json.loads(json_report_path.read_text())

        # Total reads generated was 100
        total_reads = report["summary"]["before_filtering"]["total_reads"]
        assert total_reads == 100

        # Because we strictly bounded Phred to [0, 7], all reads should fail Q15
        passed_reads = report["summary"]["after_filtering"]["total_reads"]
        failed_due_to_quality = report["filtering_result"]["low_quality_reads"]

        assert passed_reads == 0, "fastp surprisingly passed some low-quality reads"
        assert failed_due_to_quality == 100, (
            "fastp did not attribute failures to low quality"
        )

    def test_fastp_high_quality_passes(self, tmp_path: Path) -> None:
        """fastp should retain reads intentionally generated with high quality."""
        # Generate FASTQ with high quality scores
        # min_val=60 is Phred 27, max_val=73 is Phred 40.
        records = assemble_sequences(
            template_sequence=TEMPLATE,
            total_length=150,
            number_of_sequences=100,
            min_val=60,
            max_val=73,
        )

        fastq_path = tmp_path / "high_qual.fastq"
        write_fastq(str(fastq_path), records)

        json_report_path = tmp_path / "fastp.json"

        cmd = ["fastp", "-i", str(fastq_path), "-j", str(json_report_path), "-Q"]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0

        report = json.loads(json_report_path.read_text())
        passed_reads = report["summary"]["after_filtering"]["total_reads"]

        # All reads should easily pass the Q15 threshold
        assert passed_reads == 100, (
            f"fastp surprisingly filtered {100 - passed_reads} high-quality reads"
        )

    def test_fastp_sigmoidal_trimming(self, tmp_path: Path) -> None:
        """fastp sliding-window trimming truncates sigmoidal drop reads at ~100 bp."""
        records = assemble_sequences(
            template_sequence=TEMPLATE,
            total_length=150,
            number_of_sequences=100,
            decay_model="sigmoidal",
            decay_rate=1.0,
            center=100,
            min_val=38,  # Phred 5
            max_val=73,  # Phred 40
            noise_level=0.5,
        )

        fastq_path = tmp_path / "sigmoidal.fastq"
        write_fastq(str(fastq_path), records)

        json_report_path = tmp_path / "fastp.json"
        trimmed_path = tmp_path / "trimmed.fastq"

        cmd = [
            "fastp",
            "-i",
            str(fastq_path),
            "-o",
            str(trimmed_path),
            "-j",
            str(json_report_path),
            "--cut_tail",
            "--cut_tail_window_size",
            "4",
            "--cut_tail_mean_quality",
            "15",
            "-Q",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"fastp failed: {result.stderr}"

        report = json.loads(json_report_path.read_text())
        after_mean_len = report["summary"]["after_filtering"]["read1_mean_length"]

        # Mean trimmed length should be within 100 +/- 1 bp
        assert 99.0 <= after_mean_len <= 101.0, (
            f"Expected mean trimmed length ~100 bp, got {after_mean_len}"
        )
