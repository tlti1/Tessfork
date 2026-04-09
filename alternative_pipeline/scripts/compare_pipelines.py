import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ── 1. Load paper's expression data (TopHat/Cufflinks FPKM) ──────────────────
print("Loading paper's expression data...")
paper_df = pd.read_csv(
    "GSE75140_hOrg.fetal.master.data.frame.txt.gz",
    sep="\t",
    index_col=0
)
# Remove quotes from index
paper_df.index = paper_df.index.str.replace('"', '')
paper_df.columns = paper_df.columns.str.replace('"', '')
print(f"Paper data shape: {paper_df.shape}")  # cells x genes

# Filter to organoid cells only (hOrg in cell name)
organoid_mask = paper_df.index.str.contains('hOrg')
paper_org = paper_df[organoid_mask]
print(f"Paper organoid cells: {paper_org.shape[0]}")

# ── 2. Load our expression data (HISAT2/featureCounts) ───────────────────────
print("\nLoading our expression data...")
our_df = pd.read_csv(
    "all_cells_counts.txt",
    sep="\t",
    comment="#",
    index_col=0
)
our_df = our_df.drop(columns=["Chr", "Start", "End", "Strand", "Length", "gene_id"])
our_df.columns = [c.split("/")[-1].replace("_sorted.bam", "") for c in our_df.columns]

# Load metadata to filter to organoid cells
meta = pd.read_csv("SraRunTable.csv", index_col="Run")
organoid_srr = meta[meta["tissue"].isin([
    "Dissociated whole cerebral organoid",
    "Microdissected cortical-like ventricle from cerebral organoid"
])].index.tolist()
organoid_srr = [s for s in organoid_srr if s in our_df.columns]
our_org = our_df[organoid_srr].T  # cells x genes

# Log normalise our counts to make comparable to paper's log2 FPKM
our_norm = our_org.div(our_org.sum(axis=1), axis=0) * 1e4
our_log = np.log1p(our_norm)
print(f"Our organoid cells: {our_log.shape[0]}")

# ── 3. Find common genes ──────────────────────────────────────────────────────
common_genes = list(set(paper_org.columns) & set(our_log.columns))
print(f"\nGenes in paper: {len(paper_org.columns)}")
print(f"Genes in our data: {len(our_log.columns)}")
print(f"Common genes: {len(common_genes)}")

# ── 4. Compare top 10 expressed genes ────────────────────────────────────────
print("\n" + "="*60)
print("TOP 10 EXPRESSED GENES COMPARISON")
print("="*60)

# Paper top 10 (mean expression across organoid cells)
paper_mean = paper_org[common_genes].mean().sort_values(ascending=False)
our_mean = our_log[common_genes].mean().sort_values(ascending=False)

paper_top10 = paper_mean.head(10)
our_top10 = our_mean.head(10)

print("\nPaper (TopHat/Cufflinks log2 FPKM):")
for gene, val in paper_top10.items():
    print(f"  {gene}: {val:.3f}")

print("\nOurs (HISAT2/featureCounts log-normalised counts):")
for gene, val in our_top10.items():
    print(f"  {gene}: {val:.3f}")

# Overlap
paper_top_genes = set(paper_top10.index)
our_top_genes = set(our_top10.index)
overlap = paper_top_genes & our_top_genes
print(f"\nOverlap in top 10: {len(overlap)} genes")
print(f"Shared genes: {overlap}")
print(f"Only in paper: {paper_top_genes - our_top_genes}")
print(f"Only in ours: {our_top_genes - paper_top_genes}")

# ── 5. Compare top 50 expressed genes ────────────────────────────────────────
paper_top50 = set(paper_mean.head(50).index)
our_top50 = set(our_mean.head(50).index)
overlap50 = paper_top50 & our_top50
print(f"\nOverlap in top 50: {len(overlap50)} genes ({len(overlap50)/50*100:.1f}%)")

# ── 6. Check key marker genes in both ────────────────────────────────────────
print("\n" + "="*60)
print("KEY MARKER GENE EXPRESSION COMPARISON")
print("="*60)

marker_genes = {
    'Cortical NPC': ['FOXG1', 'NFIA', 'NFIB'],
    'Cortical neuron': ['NEUROD6', 'BCL11B', 'DCX'],
    'Non-cortical NPC': ['OTX2', 'FABP7'],
    'Hem/RSPO': ['WNT2B', 'RSPO2', 'RSPO3'],
    'Mesenchyme': ['COL3A1', 'LUM', 'DCN'],
    'Cycling': ['MKI67', 'TOP2A'],
}

print(f"\n{'Gene':<12} {'Cell type':<20} {'Paper mean':<15} {'Our mean':<15} {'Direction match'}")
print("-"*70)

for cell_type, genes in marker_genes.items():
    for gene in genes:
        if gene not in common_genes:
            continue
        p_mean = paper_org[gene].mean()
        o_mean = float(our_log[gene].mean())
        # Both should be > 0 for expressed genes
        match = "✓" if (p_mean > 0) == (o_mean > 0) else "✗"
        print(f"{gene:<12} {cell_type:<20} {p_mean:<15.3f} {o_mean:<15.3f} {match}")

# ── 7. Correlation between pipelines ─────────────────────────────────────────
print("\n" + "="*60)
print("CORRELATION BETWEEN PIPELINES")
print("="*60)

# Mean expression per gene
paper_gene_means = paper_org[common_genes].mean()
our_gene_means = our_log[common_genes].mean()

correlation = paper_gene_means.corr(our_gene_means, method='spearman')
print(f"\nSpearman correlation of mean gene expression: {correlation:.3f}")

# ── 8. Plot comparison figures ────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot 1: Top 10 genes bar chart comparison
ax = axes[0]
top_genes_union = list(paper_top_genes | our_top_genes)
x = np.arange(len(top_genes_union))
width = 0.35

paper_vals = [paper_mean.get(g, 0) for g in top_genes_union]
our_vals = [our_mean.get(g, 0) for g in top_genes_union]

# Normalise to 0-1 for visual comparison
paper_vals_norm = np.array(paper_vals) / max(paper_vals) if max(paper_vals) > 0 else paper_vals
our_vals_norm = np.array(our_vals) / max(our_vals) if max(our_vals) > 0 else our_vals

bars1 = ax.bar(x - width/2, paper_vals_norm, width,
               label='Paper (TopHat/Cufflinks)', color='steelblue', alpha=0.8)
bars2 = ax.bar(x + width/2, our_vals_norm, width,
               label='Ours (HISAT2/featureCounts)', color='coral', alpha=0.8)

ax.set_xlabel('Gene')
ax.set_ylabel('Normalised mean expression')
ax.set_title('Top 10 expressed genes\n(normalised for comparison)')
ax.set_xticks(x)
ax.set_xticklabels(top_genes_union, rotation=45, ha='right', fontsize=8)
ax.legend(fontsize=8)

# Plot 2: Scatter plot of mean gene expression
ax = axes[1]
ax.scatter(paper_gene_means, our_gene_means, alpha=0.3, s=5, color='grey')

# Highlight marker genes
colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown']
for (cell_type, genes), color in zip(marker_genes.items(), colors):
    for gene in genes:
        if gene in common_genes:
            ax.scatter(paper_gene_means[gene], our_gene_means[gene],
                      s=60, color=color, zorder=5, label=gene)
            ax.annotate(gene,
                       (paper_gene_means[gene], our_gene_means[gene]),
                       fontsize=7, xytext=(3, 3),
                       textcoords='offset points')

ax.set_xlabel('Paper mean expression (log2 FPKM)')
ax.set_ylabel('Our mean expression (log-normalised counts)')
ax.set_title(f'Gene expression correlation\nSpearman r = {correlation:.3f}')

# Plot 3: Top 50 gene overlap
ax = axes[2]
only_paper = len(paper_top50 - our_top50)
only_ours = len(our_top50 - paper_top50)
shared = len(overlap50)

bars = ax.bar(['Only in paper', 'Shared', 'Only in ours'],
              [only_paper, shared, only_ours],
              color=['steelblue', 'green', 'coral'],
              alpha=0.8, edgecolor='black')

for bar, val in zip(bars, [only_paper, shared, only_ours]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            str(val), ha='center', va='bottom', fontweight='bold')

ax.set_ylabel('Number of genes')
ax.set_title('Top 50 expressed genes\noverlap between pipelines')
ax.set_ylim(0, 55)

plt.tight_layout()
fig.savefig("figures/pipeline_comparison.png", dpi=150, bbox_inches='tight')
plt.close()
print("\nSaved: figures/pipeline_comparison.png")

print("\nDone!")
