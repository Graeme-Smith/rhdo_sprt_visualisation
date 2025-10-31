"""
Animate SNP classification along genomic coordinates.

Visualizes SNPs appearing along a genomic region and their classification
into RHDO SNP types (Type 1-5) based on maternal and paternal genotypes.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from typing import Tuple, Optional
from dataclasses import dataclass


@dataclass
class SNP:
    """Represents a SNP with genotype information."""
    position: int
    maternal_hap1: str
    maternal_hap2: str
    paternal_hap1: str
    paternal_hap2: str
    snp_type: Optional[int] = None
    subtype: Optional[str] = None
    

def classify_snp(maternal_hap1: str, maternal_hap2: str, 
                 paternal_hap1: str, paternal_hap2: str) -> Tuple[int, Optional[str]]:
    """
    Classify SNP into RHDO types based on maternal and paternal genotypes.
    
    Returns:
        (type, subtype) where type is 1-5 and subtype is letter/number identifier
    """
    mat_alleles = {maternal_hap1, maternal_hap2}
    pat_alleles = {paternal_hap1, paternal_hap2}
    
    # Type 1: Parents homozygous for different genotypes
    if len(mat_alleles) == 1 and len(pat_alleles) == 1 and mat_alleles != pat_alleles:
        return (1, 'A' if maternal_hap1 == 'A' else 'B')
    
    # Type 2: Parents homozygous for same genotype
    if len(mat_alleles) == 1 and len(pat_alleles) == 1 and mat_alleles == pat_alleles:
        return (2, 'A' if maternal_hap1 == 'A' else 'B')
    
    # Type 3: Mother homozygous, father heterozygous
    if len(mat_alleles) == 1 and len(pat_alleles) == 2:
        if maternal_hap1 == 'A' and maternal_hap2 == 'A':
            if paternal_hap1 == 'A' and paternal_hap2 == 'B':
                return (3, 'A')
            if paternal_hap1 == 'B' and paternal_hap2 == 'A':
                return (3, 'A')
        if maternal_hap1 == 'B' and maternal_hap2 == 'A':
            if paternal_hap1 == 'A' and paternal_hap2 == 'B':
                return (3, 'B')
        if maternal_hap1 == 'B' and maternal_hap2 == 'B':
            if 'A' in pat_alleles and 'B' in pat_alleles:
                return (3, 'C' if paternal_hap1 == 'A' else 'D')
        return (3, 'A')  # Default for type 3
    
    # Type 4a: Mother heterozygous (A/B), father homozygous
    if len(mat_alleles) == 2 and len(pat_alleles) == 1:
        if 'A' in mat_alleles and 'B' in mat_alleles:
            if 'A' in pat_alleles:
                return (4, 'a1' if maternal_hap1 == 'A' else 'a2')
            else:
                return (4, 'b1' if maternal_hap1 == 'B' else 'b2')
    
    # Type 5: Both parents heterozygous for same genotype
    if len(mat_alleles) == 2 and len(pat_alleles) == 2 and mat_alleles == pat_alleles:
        return (5, 'A' if 'A' in mat_alleles else 'B')
    
    return (0, None)  # Unclassified


def simulate_snps(chr_start: int, chr_end: int, num_snps: int, 
                  rng: np.random.Generator) -> list[SNP]:
    """
    Simulate SNPs along a genomic region.
    
    Args:
        chr_start: Start position of genomic region
        chr_end: End position of genomic region
        num_snps: Number of SNPs to simulate
        rng: Random number generator
        
    Returns:
        List of SNP objects with random positions and genotypes
    """
    positions = sorted(rng.integers(chr_start, chr_end, size=num_snps))
    snps = []
    
    # Common allele types for simulation
    alleles = ['A', 'B']
    
    for pos in positions:
        # Randomly assign genotypes (simplified - in practice would use real data)
        mat_h1 = rng.choice(alleles)
        mat_h2 = rng.choice(alleles)
        pat_h1 = rng.choice(alleles)
        pat_h2 = rng.choice(alleles)
        
        snp_type, subtype = classify_snp(mat_h1, mat_h2, pat_h1, pat_h2)
        
        snp = SNP(
            position=pos,
            maternal_hap1=mat_h1,
            maternal_hap2=mat_h2,
            paternal_hap1=pat_h1,
            paternal_hap2=pat_h2,
            snp_type=snp_type,
            subtype=subtype
        )
        snps.append(snp)
    
    return snps


def load_snps_from_csv(csv_path: str) -> list[SNP]:
    """
    Load SNPs from CSV file.
    
    Expected CSV format:
        position, maternal_hap1, maternal_hap2, paternal_hap1, paternal_hap2
        or
        position, mat_hap1, mat_hap2, pat_hap1, pat_hap2
    
    Args:
        csv_path: Path to CSV file
        
    Returns:
        List of SNP objects
    """
    import csv
    
    snps = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Try different column name variations
            pos = int(row.get('position', row.get('pos', 0)))
            mat_h1 = row.get('maternal_hap1', row.get('mat_hap1', row.get('maternal_h1', '')))
            mat_h2 = row.get('maternal_hap2', row.get('mat_hap2', row.get('maternal_h2', '')))
            pat_h1 = row.get('paternal_hap1', row.get('pat_hap1', row.get('paternal_h1', '')))
            pat_h2 = row.get('paternal_hap2', row.get('pat_hap2', row.get('paternal_h2', '')))
            
            snp_type, subtype = classify_snp(mat_h1, mat_h2, pat_h1, pat_h2)
            
            snp = SNP(
                position=pos,
                maternal_hap1=mat_h1,
                maternal_hap2=mat_h2,
                paternal_hap1=pat_h1,
                paternal_hap2=pat_h2,
                snp_type=snp_type,
                subtype=subtype
            )
            snps.append(snp)
    
    return snps


def animate_snp_classification(
    snps: list[SNP],
    chr_name: str = "chr1",
    save_path: Optional[str] = None,
    fps: int = 5,
    dpi: int = 100,
):
    """
    Animate SNP classification along genomic coordinates.
    
    Args:
        snps: List of SNP objects to display
        chr_name: Chromosome name for labeling
        save_path: Optional path to save animation (.gif or .mp4)
        fps: Frames per second
        dpi: Resolution for saved animation
    """
    if not snps:
        raise ValueError("No SNPs provided")
    
    # Color mapping for SNP types
    type_colors = {
        0: 'gray',      # Unclassified
        1: '#5cb85c',   # Green - Type 1 (fetal fraction)
        2: '#999999',   # Gray - Type 2 (QC only)
        3: '#2c7be5',   # Blue - Type 3 (paternal inheritance)
        4: '#f0ad4e',   # Orange - Type 4 (maternal inheritance)
        5: '#d9534f',   # Red - Type 5 (consanguineous)
    }
    
    type_labels = {
        0: 'Unclassified',
        1: 'Type 1: Fetal fraction',
        2: 'Type 2: QC only',
        3: 'Type 3: Paternal inheritance',
        4: 'Type 4: Maternal inheritance',
        5: 'Type 5: Consanguineous',
    }
    
    fig, ax = plt.subplots(figsize=(14, 6))
    
    positions = [snp.position for snp in snps]
    chr_start = min(positions)
    chr_end = max(positions)
    
    # Initialize empty scatter plot
    scatter = ax.scatter([], [], s=100, alpha=0.7, edgecolors='black', linewidths=1)
    
    ax.set_xlim(chr_start - (chr_end - chr_start) * 0.05, 
                chr_end + (chr_end - chr_start) * 0.05)
    ax.set_ylim(-0.5, 5.5)
    ax.set_xlabel(f"Genomic position ({chr_name})", fontsize=12)
    ax.set_ylabel("SNP Type", fontsize=12)
    ax.set_yticks([0, 1, 2, 3, 4, 5])
    ax.set_yticklabels(['Unclassified', 'Type 1', 'Type 2', 'Type 3', 'Type 4', 'Type 5'])
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_title(f"SNP Classification along {chr_name}", fontsize=14, fontweight='bold')
    
    # Create legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=type_colors[i], label=type_labels[i]) 
                      for i in range(1, 6)]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
    
    # Info text
    info_text = ax.text(0.02, 0.98, "", transform=ax.transAxes, 
                       fontsize=10, verticalalignment='top',
                       bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.9, ec="#999"))
    
    # Type counts
    type_counts = {i: 0 for i in range(6)}
    
    def init():
        scatter.set_offsets(np.empty((0, 2)))
        info_text.set_text("")
        return scatter, info_text
    
    def update(frame: int):
        nonlocal type_counts
        if frame >= len(snps):
            frame = len(snps) - 1
        
        # Get SNPs up to current frame
        visible_snps = snps[:frame + 1]
        
        # Collect positions and types
        x_pos = [snp.position for snp in visible_snps]
        y_pos = [snp.snp_type if snp.snp_type is not None else 0 for snp in visible_snps]
        colors = [type_colors.get(snp.snp_type, 'gray') for snp in visible_snps]
        
        # Update scatter plot
        scatter.set_offsets(np.column_stack([x_pos, y_pos]))
        scatter.set_color(colors)
        
        # Update counts
        type_counts = {i: 0 for i in range(6)}
        for snp in visible_snps:
            if snp.snp_type is not None:
                type_counts[snp.snp_type] = type_counts.get(snp.snp_type, 0) + 1
        
        # Update info
        info = f"SNPs processed: {frame + 1}/{len(snps)}\n"
        info += f"Type 1: {type_counts[1]}  Type 3: {type_counts[3]}\n"
        info += f"Type 2: {type_counts[2]}  Type 4: {type_counts[4]}  Type 5: {type_counts[5]}"
        info_text.set_text(info)
        
        return scatter, info_text
    
    anim = FuncAnimation(fig, update, frames=len(snps), init_func=init, 
                        interval=int(1000 / fps), blit=True, repeat=False)
    
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


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Animate SNP classification along genomic coordinates")
    parser.add_argument("--snps_csv", type=str, default=None, help="Path to CSV file with SNP data")
    parser.add_argument("--chr_start", type=int, default=1000000, help="Start position for simulation (default: 1000000)")
    parser.add_argument("--chr_end", type=int, default=2000000, help="End position for simulation (default: 2000000)")
    parser.add_argument("--num_snps", type=int, default=50, help="Number of SNPs to simulate (default: 50)")
    parser.add_argument("--chr_name", type=str, default="chr1", help="Chromosome name (default: chr1)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    parser.add_argument("--save_path", type=str, default=None, help="Output path (.gif or .mp4)")
    parser.add_argument("--fps", type=int, default=5, help="Frames per second (default: 5)")
    parser.add_argument("--dpi", type=int, default=100, help="DPI for saved animation (default: 100)")
    
    args = parser.parse_args()
    
    if args.snps_csv:
        snps = load_snps_from_csv(args.snps_csv)
    else:
        rng = np.random.default_rng(args.seed)
        snps = simulate_snps(args.chr_start, args.chr_end, args.num_snps, rng)
    
    animate_snp_classification(
        snps=snps,
        chr_name=args.chr_name,
        save_path=args.save_path,
        fps=args.fps,
        dpi=args.dpi,
    )

