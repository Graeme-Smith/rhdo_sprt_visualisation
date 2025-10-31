## SPRT Animation

This project animates the Sequential Probability Ratio Test (SPRT) using matplotlib for RHDO (Relative Haplotype Dosage) with per-locus binomial counts.

You can visualize either the cumulative log-likelihood ratio (LLR view), the cumulative fraction with curved SPRT decision boundaries (probability view), or both side-by-side.

### Install

```bash
python -m venv .venv
. .venv/Scripts/activate  # Windows PowerShell: . .venv/Scripts/Activate.ps1
pip install -r requirements.txt
```

### RHDO mode (simulate counts)

Simulate per-locus binomial counts, e.g. if haplotype A is slightly overrepresented (`p_true > 0.5`).

```bash
python sprt_animate.py --mode rhdo --view both \
  --p0 0.5 --p1 0.52 --p_true 0.515 \
  --alpha 0.05 --beta 0.1 --max_n 200 \
  --rhdo_mean_n 80 \
  --prob_ylim_min 0.3 --prob_ylim_max 0.7
```

- `p0`, `p1`: hypothesized allele/reads proportion under H0 and H1 (commonly 0.5 vs >0.5)
- `p_true`: simulation truth (e.g., 0.515)
- `rhdo_mean_n`: average per-locus coverage (Poisson). Use `--rhdo_fixed_n` to set fixed coverage instead.
- `--view prob` renders cumulative fraction vs. cumulative reads with curved boundaries.
- `--view llr` renders the cumulative LLR vs. constant log boundaries.
- `--view both` renders probability and LLR views side-by-side in one figure.
- `--prob_ylim_min/--prob_ylim_max`: set probability-view Y range (defaults 0.3–0.7).

### RHDO mode (use CSV counts)

Provide your own per-locus counts as a two-column CSV with `n,k` per row (n = total reads, k = reads supporting target haplotype):

```csv
# n,k
105,61
97,49
120,74
...
```

Run:

```bash
python sprt_animate.py --mode rhdo --view both \
  --p0 0.5 --p1 0.52 --alpha 0.05 --beta 0.1 \
  --rhdo_counts_csv my_counts.csv \
  --prob_ylim_min 0.3 --prob_ylim_max 0.7
```

### Save animation (GIF or MP4)

```bash
# GIF (requires Pillow)
python sprt_animate.py --mode rhdo --view both --save_path out.gif --dpi 150

# MP4 (requires ffmpeg on PATH)
python sprt_animate.py --mode rhdo --view both --save_path out.mp4 --fps 30 --dpi 150
```

### Notes
- Boundaries: A = (1-β)/α, B = β/(1-α); plotted as a = ln B (lower), b = ln A (upper).
- Probability-view thresholds for each cumulative n come from solving the LLR equality for the success proportion.
- RHDO binomial increment per locus: `LLR(k,n) = k ln(p1/p0) + (n-k) ln((1-p1)/(1-p0))`.
