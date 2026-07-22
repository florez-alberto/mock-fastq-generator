"""Generate publication-ready 3-panel figure showcasing all three decay models.
Emulates realistic experimental FASTQ quality profiles with raw noisy curves and no artificial mean overlays.
"""

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from mock_fastq_generator.core import assemble_sequences

def main():
    # Large, crisp publication typography
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
        'font.size': 14,
        'axes.labelsize': 15,
        'axes.titlesize': 16,
        'xtick.labelsize': 13,
        'ytick.labelsize': 13,
        'legend.fontsize': 12,
        'figure.titlesize': 17
    })

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.8), sharey=True)

    template_seq = "GCCGGCCATGGCG" + "A" * 15 + "GCTAGCTAGCTA" + "T" * 20 + "CGATCGATCGAT"
    
    # -------------------------------------------------------------
    # Panel A: Gaussian Decay Model (Realistic Noise, No Mean Line)
    # -------------------------------------------------------------
    res_gauss = assemble_sequences(
        template_sequence=template_seq,
        upstream_sequence="GCCGGCCATGGCG",
        total_length=500,
        number_of_sequences=1,
        center=100,
        min_val=40,  # Phred 7
        max_val=73,  # Phred 40
        std_dev=300,
        noise_level=0.55,  # Strong realistic experimental noise
        decay_model="gaussian",
    )
    scores_gauss = [ord(c) - 33 for c in res_gauss[0][3]]
    x_coords = np.arange(len(scores_gauss))

    ax_a = axes[0]
    ax_a.plot(x_coords, scores_gauss, color='#1f77b4', linewidth=1.2, alpha=0.88, label='Experimental Read Profile')
    ax_a.axhline(y=40, color='#d62728', linestyle=':', linewidth=1.5, alpha=0.6, label='Q40 Upper Limit')
    ax_a.axhline(y=7, color='#2ca02c', linestyle=':', linewidth=1.5, alpha=0.6, label='Q7 Lower Floor')
    ax_a.set_title('A) Parametric Gaussian Decay', fontweight='bold', pad=12)
    ax_a.set_xlabel('Base Position (Index)')
    ax_a.set_ylabel('Phred Quality Score (Q)')
    ax_a.grid(True, linestyle='--', alpha=0.35)
    ax_a.set_ylim(0, 45)
    ax_a.legend(loc='lower left', framealpha=0.9)

    # -------------------------------------------------------------
    # Panel B: Exponential 3' Read Degradation (Progressive non-linear drop, No Mean Line)
    # -------------------------------------------------------------
    res_exp = assemble_sequences(
        template_sequence=template_seq,
        upstream_sequence="GCCGGCCATGGCG",
        total_length=500,
        number_of_sequences=1,
        center=100,
        min_val=40,
        max_val=73,
        noise_level=0.50,     # Experimental noise
        decay_model="exponential",
        decay_rate=0.0008,    # Progressive decay that accelerates near 3' end
    )
    scores_exp = [ord(c) - 33 for c in res_exp[0][3]]

    ax_b = axes[1]
    ax_b.plot(x_coords, scores_exp, color='#ff7f0e', linewidth=1.2, alpha=0.88, label='Experimental Read Profile')
    ax_b.set_title('B) Exponential 3\' Read Degradation', fontweight='bold', pad=12)
    ax_b.set_xlabel('Base Position (Index)')
    ax_b.grid(True, linestyle='--', alpha=0.35)
    ax_b.legend(loc='lower left', framealpha=0.9)

    # -------------------------------------------------------------
    # Panel C: Sigmoidal Drop + Binning + Homopolymer Penalty
    # -------------------------------------------------------------
    homopolymer_template = "GCCGGCCATGGCGAAAAAAATTTTTTTCCCCCCGGGGGGGGATCGATCGATCG"
    res_sig = assemble_sequences(
        template_sequence=homopolymer_template,
        upstream_sequence="GCCGGCCATGGCG",
        total_length=500,
        number_of_sequences=1,
        center=260,
        min_val=40,
        max_val=73,
        noise_level=0.15,
        decay_model="sigmoidal",
        decay_rate=0.020,
        binned_quality=True,
        homopolymer_penalty=True,
    )
    scores_sig = [ord(c) - 33 for c in res_sig[0][3]]

    ax_c = axes[2]
    ax_c.plot(x_coords, scores_sig, color='#2ca02c', linewidth=1.8, drawstyle='steps-mid', label='Binned + Homopolymer Penalty')

    # Discrete NovaSeq bins
    bins = [2, 12, 23, 37]
    for b in bins:
        ax_c.axhline(y=b, color='gray', linestyle='-.', alpha=0.5, linewidth=1.3)
    
    ax_c.set_title('C) Sigmoidal + Binning & Homopolymers', fontweight='bold', pad=12)
    ax_c.set_xlabel('Base Position (Index)')
    ax_c.grid(True, linestyle='--', alpha=0.35)
    
    ax_c.text(495, 37.5, 'Q37', color='#222222', fontsize=12, ha='right', va='bottom', fontweight='bold')
    ax_c.text(495, 23.5, 'Q23', color='#222222', fontsize=12, ha='right', va='bottom', fontweight='bold')
    ax_c.text(495, 12.5, 'Q12', color='#222222', fontsize=12, ha='right', va='bottom', fontweight='bold')
    ax_c.text(495, 2.5, 'Q2', color='#222222', fontsize=12, ha='right', va='bottom', fontweight='bold')
    ax_c.legend(loc='lower left', framealpha=0.9)

    plt.tight_layout()
    
    out_dir = os.path.join(os.path.dirname(__file__), "..", "img")
    os.makedirs(out_dir, exist_ok=True)
    
    three_panel_path = os.path.join(out_dir, "decay_models_three_panel.png")
    smooth_example_path = os.path.join(out_dir, "smooth_curve_example.png")
    
    fig.savefig(three_panel_path, dpi=300, bbox_inches='tight')
    fig.savefig(smooth_example_path, dpi=300, bbox_inches='tight')
    print(f"Saved experimental 3-panel figure to: {three_panel_path} and {smooth_example_path}")

if __name__ == '__main__':
    main()
