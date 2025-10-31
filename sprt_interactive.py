"""
Interactive SPRT boundary explorer.

Use sliders to adjust α, β, p0, p1 and see how the boundaries change in real-time.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from sprt_animate import compute_log_boundaries, prob_boundaries_curves


def interactive_sprt_boundaries(
    alpha_init: float = 0.05,
    beta_init: float = 0.1,
    fetal_fraction_init: float = 0.02,
    max_n: int = 500,
):
    """Create an interactive plot with sliders to explore SPRT boundaries."""
    
    # Fixed p0 for RHDO (balanced haplotype representation)
    p0_fixed = 0.5
    # Calculate p1 from fetal fraction: p1 = 0.5 + fetal_fraction
    p1_init = p0_fixed + fetal_fraction_init
    
    # Initialize figure with single subplot
    fig = plt.figure(figsize=(12, 7))
    
    # Create main plot area
    ax_prob = plt.subplot(1, 1, 1)
    
    # Compute initial boundaries
    lower_a, upper_b = compute_log_boundaries(alpha_init, beta_init)
    
    # Create n values for plotting curves
    n_curve = np.linspace(1, max_n, num=300)
    t_low, t_high = prob_boundaries_curves(n_curve, lower_a, upper_b, p0_fixed, p1_init)
    
    # Plot probability view
    line_high_prob, = ax_prob.plot(n_curve, t_high, color="#5cb85c", lw=2, label="Upper threshold")
    line_low_prob, = ax_prob.plot(n_curve, t_low, color="#a64ac9", lw=2, label="Lower threshold")
    # Add horizontal line at 0.5 to represent balanced state
    line_balanced = ax_prob.axhline(0.5, color="gray", linestyle="--", lw=1.5, alpha=0.7, label="Balanced (p₀ = 0.5)")
    ax_prob.set_xlim(0, max_n)
    ax_prob.set_ylim(0.0, 1.0)
    ax_prob.set_xlabel("Total cumulative reads (n)", fontsize=11, labelpad=10)
    ax_prob.set_ylabel("Fraction for target haplotype", fontsize=11)
    ax_prob.grid(True, alpha=0.25)
    ax_prob.legend(loc="best")
    ax_prob.set_title("SPRT Probability Thresholds", fontsize=12, fontweight='bold', pad=10)
    # Adjust tick label padding
    ax_prob.tick_params(axis='x', pad=8)
    
    # Add info text (positioned at bottom right)
    info_text = ax_prob.text(0.98, 0.02, "", transform=ax_prob.transAxes, 
                            fontsize=9, verticalalignment='bottom', horizontalalignment='right',
                            bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.9, ec="#999"))
    
    def update_info(alpha, beta, p0, p1, lower_a, upper_b):
        """Update the info text with current parameter values."""
        A = (1.0 - beta) / alpha
        B = beta / (1.0 - alpha)
        fetal_fraction = p1 - p0
        info = (
            f"α = {alpha:.3f},  β = {beta:.3f}\n"
            f"Fetal fraction = {fetal_fraction:.4f}\n"
            f"p₀ = {p0:.3f} (H₀: balanced),  p₁ = {p1:.4f} (H₁: imbalanced)\n"
            f"A = (1-β)/α = {A:.3f},  B = β/(1-α) = {B:.3f}\n"
            f"a = ln(B) = {lower_a:.3f},  b = ln(A) = {upper_b:.3f}"
        )
        info_text.set_text(info)
    
    update_info(alpha_init, beta_init, p0_fixed, p1_init, lower_a, upper_b)
    
    # Adjust subplots to make room for sliders - use larger bottom margin
    # The bottom parameter controls where the subplot ends, giving space for labels and sliders
    plt.subplots_adjust(bottom=0.50, left=0.1, right=0.95, top=0.92)
    
    # Create sliders stacked vertically underneath each other
    slider_width = 0.6
    slider_height = 0.03
    slider_spacing = 0.04
    slider_x = 0.2
    # Position sliders well below the subplot and x-axis labels
    # Subplot ends at bottom=0.50, labels extend below, so start sliders around 0.35
    start_y = 0.35  # Start sliders with clearance below x-axis labels
    
    ax_alpha = plt.axes([slider_x, start_y, slider_width, slider_height])
    ax_beta = plt.axes([slider_x, start_y - slider_spacing, slider_width, slider_height])
    ax_fetal = plt.axes([slider_x, start_y - 2*slider_spacing, slider_width, slider_height])
    
    slider_alpha = Slider(ax_alpha, 'α (Type I error)', 0.01, 0.30, valinit=alpha_init, 
                         valfmt='%.3f', valstep=0.01)
    slider_beta = Slider(ax_beta, 'β (Type II error)', 0.01, 0.30, valinit=beta_init, 
                        valfmt='%.3f', valstep=0.01)
    slider_fetal = Slider(ax_fetal, 'Fetal fraction (imbalance)', 0.02, 0.20, valinit=fetal_fraction_init, 
                         valfmt='%.4f', valstep=0.005)
    
    def update(val):
        """Update plots when sliders change."""
        alpha = slider_alpha.val
        beta = slider_beta.val
        fetal_fraction = slider_fetal.val
        
        # Calculate p0 and p1 from fetal fraction
        p0 = p0_fixed  # Always 0.5 for RHDO (balanced)
        p1 = p0_fixed + fetal_fraction  # p1 increases with fetal fraction
        
        # Recompute boundaries
        lower_a, upper_b = compute_log_boundaries(alpha, beta)
        
        # Update probability curves
        t_low, t_high = prob_boundaries_curves(n_curve, lower_a, upper_b, p0, p1)
        line_low_prob.set_ydata(t_low)
        line_high_prob.set_ydata(t_high)
        
        # Update info text
        update_info(alpha, beta, p0, p1, lower_a, upper_b)
        
        fig.canvas.draw_idle()
    
    # Connect sliders to update function
    slider_alpha.on_changed(update)
    slider_beta.on_changed(update)
    slider_fetal.on_changed(update)
    
    plt.suptitle("Interactive SPRT Boundary Explorer", fontsize=14, fontweight='bold', y=0.995)
    
    plt.show()
    
    return fig


if __name__ == "__main__":
    interactive_sprt_boundaries(
        alpha_init=0.05,
        beta_init=0.1,
        fetal_fraction_init=0.02,
        max_n=500,
    )

