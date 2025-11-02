import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from typing import Tuple


def compute_log_boundaries(alpha: float, beta: float) -> tuple[float, float]:
    """Return (lower_a, upper_b) boundaries in log-space for SPRT.

    Lower boundary is a = ln(B) where B = β / (1 - α)
    Upper boundary is b = ln(A) where A = (1 - β) / α
    """
    if not (0 < alpha < 1 and 0 < beta < 1):
        raise ValueError("alpha and beta must be in (0, 1)")
    A = (1.0 - beta) / alpha
    B = beta / (1.0 - alpha)
    return math.log(B), math.log(A)


def bernoulli_llr_sample(x: int, p0: float, p1: float) -> float:
    """Per-sample log-likelihood ratio for Bernoulli.
    Returns ln( P1(x) / P0(x) ).
    """
    if not (0.0 < p0 < 1.0 and 0.0 < p1 < 1.0):
        raise ValueError("p0 and p1 must be in (0, 1)")
    if x not in (0, 1):
        raise ValueError("Bernoulli sample x must be 0 or 1")
    log_p1_over_p0_when_one = math.log(p1) - math.log(p0)
    log_q1_over_q0_when_zero = math.log(1.0 - p1) - math.log(1.0 - p0)
    return log_p1_over_p0_when_one if x == 1 else log_q1_over_q0_when_zero


def binomial_llr_count(k: int, n: int, p0: float, p1: float) -> float:
    """Log-likelihood ratio for Binomial(k; n, p).

    LLR = ln( C(n,k) p1^k (1-p1)^(n-k) / (C(n,k) p0^k (1-p0)^(n-k) ) )
        = k ln(p1/p0) + (n-k) ln((1-p1)/(1-p0))
    """
    if n < 0 or k < 0 or k > n:
        raise ValueError("Require 0 <= k <= n and n >= 0")
    if not (0.0 < p0 < 1.0 and 0.0 < p1 < 1.0):
        raise ValueError("p0 and p1 must be in (0, 1)")
    if n == 0:
        return 0.0
    return k * (math.log(p1) - math.log(p0)) + (n - k) * (math.log(1.0 - p1) - math.log(1.0 - p0))


def simulate_bernoulli(p_true: float, n: int, rng: np.random.Generator) -> np.ndarray:
    if not (0.0 < p_true < 1.0):
        raise ValueError("p_true must be in (0, 1)")
    return rng.binomial(1, p_true, size=n)


def simulate_rhdo_counts(num_loci: int, p_true: float, rng: np.random.Generator, mean_n: int, fixed_n: int | None) -> Tuple[np.ndarray, np.ndarray]:
    """Simulate per-locus binomial counts for RHDO.

    Returns (k_counts, n_totals) arrays of length num_loci.
    If fixed_n is provided, each locus uses that coverage; otherwise n ~ Poisson(mean_n) clipped to >= 1.
    """
    if not (0.0 < p_true < 1.0):
        raise ValueError("p_true must be in (0, 1)")
    if fixed_n is not None:
        if fixed_n <= 0:
            raise ValueError("fixed_n must be positive")
        n_totals = np.full(num_loci, fixed_n, dtype=int)
    else:
        n_totals = rng.poisson(lam=max(1, int(mean_n)), size=num_loci).astype(int)
        n_totals = np.maximum(n_totals, 1)
    k_counts = rng.binomial(n_totals, p_true).astype(int)
    return k_counts, n_totals


def build_llr_trajectory(samples: np.ndarray, p0: float, p1: float) -> np.ndarray:
    llr_increments = np.array([bernoulli_llr_sample(int(x), p0, p1) for x in samples], dtype=float)
    return np.concatenate(([0.0], np.cumsum(llr_increments)))


def build_llr_trajectory_counts(k_counts: np.ndarray, n_totals: np.ndarray, p0: float, p1: float) -> np.ndarray:
    llr_increments = np.array([binomial_llr_count(int(k), int(n), p0, p1) for k, n in zip(k_counts, n_totals)], dtype=float)
    return np.concatenate(([0.0], np.cumsum(llr_increments)))


def prob_boundaries_curves(n_values: np.ndarray, a: float, b: float, p0: float, p1: float) -> Tuple[np.ndarray, np.ndarray]:
    """Return lower/upper probability threshold curves t_low(n), t_high(n) from LLR bounds.

    For LLR(n, t) = n[(A-B)t + B], where A=ln(p1/p0), B=ln((1-p1)/(1-p0)),
    solve for t at LLR=a or LLR=b: t = (LLR/n - B) / (A - B).
    """
    A = math.log(p1) - math.log(p0)
    B = math.log(1.0 - p1) - math.log(1.0 - p0)
    denom = (A - B)
    n_safe = np.maximum(n_values.astype(float), 1.0)
    t_low = (a / n_safe - B) / denom
    t_high = (b / n_safe - B) / denom
    return np.clip(t_low, 0.0, 1.0), np.clip(t_high, 0.0, 1.0)


def animate_sprt(
    p0: float,
    p1: float,
    p_true: float,
    alpha: float,
    beta: float,
    max_n: int,
    seed: int | None,
    save_path: str | None,
    fps: int,
    dpi: int,
    mode: str,
    rhdo_counts_csv: str | None,
    rhdo_mean_n: int,
    rhdo_fixed_n: int | None,
    view: str = "llr",
    prob_ylim_min: float = 0.3,
    prob_ylim_max: float = 0.7,
):
    # Simulation / data preparation
    rng = np.random.default_rng(seed)
    lower_a, upper_b = compute_log_boundaries(alpha, beta)

    if mode == "bernoulli":
        samples = simulate_bernoulli(p_true, max_n, rng)
        llr_path = build_llr_trajectory(samples, p0, p1)
        k_cum = np.cumsum(samples)
        n_cum = np.arange(1, len(samples) + 1)
        frac_cum = k_cum / n_cum
        title_base = "SPRT (Bernoulli)"
    elif mode == "rhdo":
        if rhdo_counts_csv:
            try:
                data = np.genfromtxt(rhdo_counts_csv, delimiter=",", dtype=float)
            except Exception as exc:
                raise RuntimeError(f"Failed to read CSV '{rhdo_counts_csv}': {exc}")
            if data.ndim == 1:
                if data.size != 2:
                    raise ValueError("CSV must have 2 columns: n,k (per row)")
                data = data.reshape(1, 2)
            if data.shape[1] != 2:
                raise ValueError("CSV must have exactly 2 columns: n,k (per row)")
            n_totals = np.asarray(data[:, 0], dtype=int)
            k_counts = np.asarray(data[:, 1], dtype=int)
            if max_n is not None and len(n_totals) > max_n:
                n_totals = n_totals[:max_n]
                k_counts = k_counts[:max_n]
        else:
            k_counts, n_totals = simulate_rhdo_counts(max_n, p_true, rng, mean_n=rhdo_mean_n, fixed_n=rhdo_fixed_n)
        llr_path = build_llr_trajectory_counts(k_counts, n_totals, p0, p1)
        k_cum = np.cumsum(k_counts)
        n_cum = np.cumsum(n_totals)
        frac_cum = k_cum / np.maximum(n_cum, 1)
        title_base = "SPRT (RHDO)"
    else:
        raise ValueError("Unsupported mode. Use 'bernoulli' or 'rhdo'.")

    # Decision detection on LLR path
    decision_index = None
    decision_label = None
    for i in range(1, len(llr_path)):
        if llr_path[i] >= upper_b:
            decision_index = i
            decision_label = "Reject H0 (accept H1)"
            break
        if llr_path[i] <= lower_a:
            decision_index = i
            decision_label = "Reject H1 (accept H0)"
            break

    # Compute a unified ending frame (inclusive index), then pass frames=end+1
    if view == "llr":
        end_frame = (decision_index if decision_index is not None else len(llr_path) - 1)
    elif view == "prob":
        end_frame = (decision_index if decision_index is not None else len(n_cum))
    else:  # both
        end_llr = (decision_index if decision_index is not None else len(llr_path) - 1)
        end_prob = (decision_index if decision_index is not None else len(n_cum))
        end_frame = min(end_llr, end_prob)

    plt.close('all')

    if view == "prob":
        fig, ax = plt.subplots(figsize=(8, 5))
        fig.suptitle(f"{title_base} — Probability view")
        x_vals = n_cum
        n_curve = np.linspace(1, float(x_vals[-1] if len(x_vals) > 0 else max_n), num=300)
        t_low, t_high = prob_boundaries_curves(n_curve, lower_a, upper_b, p0, p1)
        ax.plot(n_curve, t_high, color="#5cb85c", lw=2, label="Upper threshold")
        ax.plot(n_curve, t_low, color="#a64ac9", lw=2, label="Lower threshold")
        scatter, = ax.plot([], [], "^", color="#2c7be5", label="Cumulative fraction")
        ax.set_xlim(0, float(n_curve[-1]) if len(n_curve) else max_n)
        ax.set_ylim(prob_ylim_min, prob_ylim_max)
        ax.set_xlabel("Total no. of cumulative reads")
        ax.set_ylabel("Fraction of reads for target haplotype")
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best")
        decision_text = ax.text(0.02, 0.95, "", transform=ax.transAxes, va="top", ha="left",
                                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec="#999"))
        def init():
            scatter.set_data([], [])
            decision_text.set_text("")
            return scatter, decision_text
        def update(frame: int):
            xs = x_vals[: frame + 1] if frame >= 0 else []
            ys = frac_cum[: frame + 1] if frame >= 0 else []
            scatter.set_data(xs, ys)
            if decision_index is not None and frame >= decision_index:
                decision_text.set_text(f"Decision at n={decision_index}:\n{decision_label}")
            else:
                decision_text.set_text("Running...")
            return scatter, decision_text
        total_frames = max(1, int(end_frame) + 1)
        anim = FuncAnimation(fig, update, frames=total_frames, init_func=init, interval=int(1000 / fps), blit=True, repeat=False)
    elif view == "llr":
        fig, ax = plt.subplots(figsize=(8, 5))
        fig.suptitle(f"{title_base} — LLR view")
        x_data = np.arange(len(llr_path))
        line_llr, = ax.plot([], [], color="#2c7be5", lw=2, label="Cumulative LLR")
        marker_current, = ax.plot([], [], "o", color="#2c7be5")
        lower_line = ax.axhline(lower_a, color="#d9534f", lw=1.5, ls="--", label="Lower boundary: a = ln(β/(1-α))")
        upper_line = ax.axhline(upper_b, color="#5cb85c", lw=1.5, ls="--", label="Upper boundary: b = ln((1-β)/α)")
        ax.set_xlim(0, min(len(x_data) - 1, end_frame))
        min_y = min(lower_a, float(np.min(llr_path)))
        max_y = max(upper_b, float(np.max(llr_path)))
        pad = 0.2 * (max_y - min_y + 1e-6)
        ax.set_ylim(min_y - pad, max_y + pad)
        ax.fill_between([0, end_frame], upper_b, ax.get_ylim()[1], color="#5cb85c", alpha=0.1, step="pre")
        ax.fill_between([0, end_frame], ax.get_ylim()[0], lower_a, color="#d9534f", alpha=0.1, step="pre")
        ax.set_xlabel("Index")
        ax.set_ylabel("Cumulative log-likelihood ratio")
        ax.grid(True, alpha=0.25)
        ax.legend(handles=[line_llr, upper_line, lower_line], loc="best")
        decision_text = ax.text(0.02, 0.95, "", transform=ax.transAxes, va="top", ha="left", bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec="#999"))
        def init():
            line_llr.set_data([], [])
            marker_current.set_data([], [])
            decision_text.set_text("")
            return line_llr, marker_current, decision_text
        def update(frame: int):
            x = x_data[: frame + 1]
            y = llr_path[: frame + 1]
            line_llr.set_data(x, y)
            if len(x) > 0:
                marker_current.set_data([x[-1]], [y[-1]])
            else:
                marker_current.set_data([], [])
            if decision_index is not None and frame >= decision_index:
                decision_text.set_text(f"Decision at n={decision_index}:\n{decision_label}")
            else:
                decision_text.set_text("Running...")
            return line_llr, marker_current, decision_text
        total_frames = max(1, int(end_frame) + 1)
        anim = FuncAnimation(fig, update, frames=total_frames, init_func=init, interval=int(1000 / fps), blit=True, repeat=False)
    elif view == "both":
        fig, (ax_prob, ax_llr) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
        fig.suptitle(f"{title_base} — Probability and LLR views")
        # Probability panel
        x_vals = n_cum
        n_curve = np.linspace(1, float(x_vals[-1] if len(x_vals) > 0 else max_n), num=300)
        t_low, t_high = prob_boundaries_curves(n_curve, lower_a, upper_b, p0, p1)
        ax_prob.plot(n_curve, t_high, color="#5cb85c", lw=2, label="Upper threshold")
        ax_prob.plot(n_curve, t_low, color="#a64ac9", lw=2, label="Lower threshold")
        scatter, = ax_prob.plot([], [], "^", color="#2c7be5", label="Cumulative fraction")
        ax_prob.set_xlim(0, float(n_curve[-1]) if len(n_curve) else max_n)
        ax_prob.set_ylim(prob_ylim_min, prob_ylim_max)
        ax_prob.set_xlabel("Total cumulative reads")
        ax_prob.set_ylabel("Fraction for target haplotype")
        ax_prob.grid(True, alpha=0.25)
        ax_prob.legend(loc="best")
        # LLR panel
        x_data = np.arange(len(llr_path))
        line_llr, = ax_llr.plot([], [], color="#2c7be5", lw=2, label="Cumulative LLR")
        marker_current, = ax_llr.plot([], [], "o", color="#2c7be5")
        lower_line = ax_llr.axhline(lower_a, color="#d9534f", lw=1.5, ls="--", label="Lower boundary: a = ln(β/(1-α))")
        upper_line = ax_llr.axhline(upper_b, color="#5cb85c", lw=1.5, ls="--", label="Upper boundary: b = ln((1-β)/α)")
        ax_llr.set_xlim(0, end_frame)
        min_y = min(lower_a, float(np.min(llr_path)))
        max_y = max(upper_b, float(np.max(llr_path)))
        pad = 0.2 * (max_y - min_y + 1e-6)
        ax_llr.set_ylim(min_y - pad, max_y + pad)
        ax_llr.fill_between([0, end_frame], upper_b, ax_llr.get_ylim()[1], color="#5cb85c", alpha=0.1, step="pre")
        ax_llr.fill_between([0, end_frame], ax_llr.get_ylim()[0], lower_a, color="#d9534f", alpha=0.1, step="pre")
        ax_llr.set_xlabel("Index")
        ax_llr.set_ylabel("Cumulative LLR")
        ax_llr.grid(True, alpha=0.25)
        ax_llr.legend(handles=[line_llr, upper_line, lower_line], loc="best")
        decision_text = ax_llr.text(0.02, 0.95, "", transform=ax_llr.transAxes, va="top", ha="left", bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec="#999"))
        def init():
            scatter.set_data([], [])
            line_llr.set_data([], [])
            marker_current.set_data([], [])
            decision_text.set_text("")
            return scatter, line_llr, marker_current, decision_text
        def update(frame: int):
            xs = x_vals[: frame + 1] if frame >= 0 else []
            ys = frac_cum[: frame + 1] if frame >= 0 else []
            scatter.set_data(xs, ys)
            x = np.arange(len(llr_path))[: frame + 1]
            y = llr_path[: frame + 1]
            line_llr.set_data(x, y)
            if len(x) > 0:
                marker_current.set_data([x[-1]], [y[-1]])
            else:
                marker_current.set_data([], [])
            if decision_index is not None and frame >= decision_index:
                decision_text.set_text(f"Decision at n={decision_index}:\n{decision_label}")
            else:
                decision_text.set_text("Running...")
            return scatter, line_llr, marker_current, decision_text
        total_frames = max(1, int(end_frame) + 1)
        anim = FuncAnimation(fig, update, frames=total_frames, init_func=init, interval=int(1000 / fps), blit=True, repeat=False)
    else:
        raise ValueError("Invalid view. Use 'llr', 'prob', or 'both'.")

    if save_path:
        ext = save_path.lower().rsplit('.', 1)
        ext = ext[1] if len(ext) == 2 else ""
        if ext in ("gif",):
            anim.save(save_path, writer="pillow", dpi=dpi, fps=fps)
        elif ext in ("mp4", "m4v"):
            anim.save(save_path, writer="ffmpeg", dpi=dpi, fps=fps)
        else:
            print(f"Unrecognized extension '.{ext}'. Use .gif or .mp4. Showing interactive window instead.")
            plt.show()
    else:
        plt.show()

    return anim


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Animate SPRT for Bernoulli or RHDO using matplotlib")
    parser.add_argument("--mode", type=str, choices=["bernoulli", "rhdo"], default="bernoulli", help="Which experiment model to animate")
    parser.add_argument("--view", type=str, choices=["llr", "prob", "both"], default="llr", help="Visualization: LLR view, probability view, or both")
    parser.add_argument("--p0", type=float, default=0.5, help="Null hypothesis probability p0")
    parser.add_argument("--p1", type=float, default=0.7, help="Alternative hypothesis probability p1")
    parser.add_argument("--p_true", type=float, default=0.6, help="Data-generating probability")
    parser.add_argument("--alpha", type=float, default=0.05, help="Type I error target")
    parser.add_argument("--beta", type=float, default=0.1, help="Type II error target")
    parser.add_argument("--max_n", type=int, default=200, help="Maximum number of samples / loci / frames")
    parser.add_argument("--seed", type=int, default=7, help="RNG seed")
    parser.add_argument("--save_path", type=str, default=None, help="Optional output path (.gif or .mp4)")
    parser.add_argument("--fps", type=int, default=20, help="Frames per second for animation")
    parser.add_argument("--dpi", type=int, default=150, help="DPI for saved animation frames")
    parser.add_argument("--rhdo_counts_csv", type=str, default=None, help="Path to CSV with two columns per row: n,k")
    parser.add_argument("--rhdo_mean_n", type=int, default=100, help="Mean coverage per locus when simulating (Poisson)")
    parser.add_argument("--rhdo_fixed_n", type=int, default=None, help="If set, use fixed coverage per locus instead of Poisson")
    parser.add_argument("--prob_ylim_min", type=float, default=0.3, help="Probability view Y-axis min (default 0.3)")
    parser.add_argument("--prob_ylim_max", type=float, default=0.7, help="Probability view Y-axis max (default 0.7)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    animate_sprt(
        p0=args.p0,
        p1=args.p1,
        p_true=args.p_true,
        alpha=args.alpha,
        beta=args.beta,
        max_n=args.max_n,
        seed=args.seed,
        save_path=args.save_path,
        fps=args.fps,
        dpi=args.dpi,
        mode=args.mode,
        rhdo_counts_csv=args.rhdo_counts_csv,
        rhdo_mean_n=args.rhdo_mean_n,
        rhdo_fixed_n=args.rhdo_fixed_n,
        view=args.view,
        prob_ylim_min=args.prob_ylim_min,
        prob_ylim_max=args.prob_ylim_max,
    )
