## SPRT Animation

This project animates the Sequential Probability Ratio Test (SPRT) using matplotlib for RHDO (Relative Haplotype Dosage) with per-locus binomial counts.

You can visualize either the cumulative log-likelihood ratio (LLR view), the cumulative fraction with curved SPRT decision boundaries (probability view), or both side-by-side.

### Install

```bash
python -m venv .venv
. .venv/Scripts/activate  # Windows PowerShell: . .venv/Scripts/Activate.ps1
pip install -r requirements.txt
```

**Note:** For the interactive boundary explorer in Jupyter notebooks, `ipympl` is required and is included in the requirements. If using Jupyter, you may need to restart the kernel after installation.

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

### SNP Classification Animation

Visualize SNPs along genomic coordinates being classified into RHDO SNP types:

**Simulate SNPs:**
```bash
python snp_classifier_animate.py --chr_start 1000000 --chr_end 2000000 --num_snps 50 --chr_name chr1
```

**Load from CSV:**
```bash
python snp_classifier_animate.py --snps_csv my_snps.csv --chr_name chr1
```

CSV format (with header):
```csv
position, maternal_hap1, maternal_hap2, paternal_hap1, paternal_hap2
1000123, A, A, B, B
1000456, A, B, A, A
...
```

**Save animation:**
```bash
# GIF
python snp_classifier_animate.py --num_snps 50 --save_path snp_classification.gif --fps 5

# MP4
python snp_classifier_animate.py --num_snps 50 --save_path snp_classification.mp4 --fps 5
```

The animation shows:
- SNPs appearing sequentially along genomic coordinates
- Color-coded classification by SNP type (Type 1-5)
- Real-time counts of each SNP type
- Legend and info panel showing classification status

### Interactive boundary explorer

Explore how changing parameters (α, β, p₀, p₁) affects the SPRT boundaries using interactive sliders:

**Standalone script:**
```bash
python sprt_interactive.py
```

**In Jupyter notebook:**
```python
%matplotlib widget
from sprt_interactive import interactive_sprt_boundaries
interactive_sprt_boundaries()
```

This opens an interactive window with:
- **Probability view**: Shows how the curved thresholds change with n (total reads)
- **Sliders**: Adjust α (Type I error), β (Type II error), p₀ (null hypothesis), and p₁ (alternative hypothesis) in real-time
- **Info panel**: Displays current parameter values and computed boundaries

**Note:** For Jupyter notebooks, `ipympl>=0.9.0` is required (included in requirements.txt). After installing, restart your Jupyter kernel.

### RHDO SNP Types

For RHDO analysis, different SNP types provide different information:

| Type | Subtype | Maternal Genotype<br/>HapI | Maternal Genotype<br/>HapII | Paternal Genotype<br/>HapI | Paternal Genotype<br/>HapII | Description | Useful Information |
|------|---------|---------------------------|----------------------------|---------------------------|----------------------------|-------------|-------------------|
| 1 | A | A | A | B | B | Parents homozygous for different genotypes | Fetal fraction (non-specific detection of paternal allele) |
| 1 | B | B | B | A | A | Parents homozygous for different genotypes | Fetal fraction (non-specific detection of paternal allele) |
| 2 | A | A | A | A | A | Parents homozygous for same genotype | Sequencing error rate (quality control) |
| 2 | B | B | B | B | B | Parents homozygous for same genotype | Sequencing error rate (quality control) |
| 3 | A | A | A | A | B | Mother homozygous and father heterozygous | Paternal inheritance (specific detection of paternal allele), fetal fraction |
| 3 | B | B | A | A | B | Mother homozygous and father heterozygous | Paternal inheritance (specific detection of paternal allele), fetal fraction |
| 3 | C | B | B | A | B | Mother homozygous and father heterozygous | Paternal inheritance (specific detection of paternal allele), fetal fraction |
| 3 | D | B | B | B | A | Mother homozygous and father heterozygous | Paternal inheritance (specific detection of paternal allele), fetal fraction |
| 4a | 1 | A | B | A | A | Mother heterozygous and father homozygous | Maternal inheritance |
| 4a | 2 | B | A | B | B | Mother heterozygous and father homozygous | Maternal inheritance |
| 4b | 1 | B | A | A | A | Mother heterozygous and father homozygous | Maternal inheritance |
| 4b | 2 | A | B | B | B | Mother heterozygous and father homozygous | Maternal inheritance |
| 5 | A | A | A | B | A | Parents heterozygous for same genotype | Consanguineous |
| 5 | B | B | B | A | B | Parents heterozygous for same genotype | Consanguineous |

**Type 1 and Type 3** SNPs are most useful for fetal fraction detection, as they allow detection of paternal alleles in maternal plasma. The SPRT test in this package is designed to detect imbalances in haplotype representation using these informative SNPs.

### Notes
- Boundaries: A = (1-β)/α, B = β/(1-α); plotted as a = ln B (lower), b = ln A (upper).
- Probability-view thresholds for each cumulative n come from solving the LLR equality for the success proportion.
- RHDO binomial increment per locus: `LLR(k,n) = k ln(p1/p0) + (n-k) ln((1-p1)/(1-p0))`.
