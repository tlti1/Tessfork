
#Scanpy Analysis of Salmon-Quantified Organoid Data
#Alternative pipeline to paper's TopHat/Cufflinks + Seurat/Monocle approach

#This script analyzes 16 fetal neocortex samples (12 weeks post-conception)
#quantified using Salmon instead of the paper's Cufflinks pipeline.

#Last edited: Tess 27/02/2026 10:14am

#NOTE: Claude AI was used for code formatting, debugging, and troubleshooting technical issues. 
#All biological analysis and interpretation performed by researcher.

# ═══════════════════════════════════════════════════════════════════════════
# IMPORT LIBRARIES
# ═══════════════════════════════════════════════════════════════════════════

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

print("Libraries loaded successfully")

# ═══════════════════════════════════════════════════════════════════════════
# LOAD SALMON DATA
# ═══════════════════════════════════════════════════════════════════════════

df = pd.read_csv('salmon_combined_expression.txt', sep='\t', index_col=0)
print(f"Loaded data: {df.shape[0]} cells × {df.shape[1]} genes")

# Convert to numeric and handle missing values
df = df.apply(pd.to_numeric, errors='coerce')
df = df.fillna(0)

# Create AnnData object (Scanpy's data format)
adata = sc.AnnData(df)
print(f"Created AnnData object: {adata}")

# ═══════════════════════════════════════════════════════════════════════════
# ADD METADATA
# ═══════════════════════════════════════════════════════════════════════════

# All samples are from fetal neocortex 12 weeks post-conception, chip 1
adata.obs['experiment'] = 'fetal_12wpc_c1'
adata.obs['tissue'] = 'fetal_neocortex'
adata.obs['stage'] = '12_weeks'

print("\nSample metadata:")
print(adata.obs['experiment'].value_counts())

# ═══════════════════════════════════════════════════════════════════════════
# PREPROCESSING
# ═══════════════════════════════════════════════════════════════════════════

# Filter genes expressed in fewer than 2 cells
print("Filtering lowly expressed genes...")
sc.pp.filter_genes(adata, min_cells=2)
print(f"After filtering: {adata.n_vars} genes remaining")

# Identify highly variable genes (top 2000)
print("Finding highly variable genes...")
sc.pp.highly_variable_genes(adata, 
                             n_top_genes=2000,
                             flavor='cell_ranger')
print(f"Selected {adata.var['highly_variable'].sum()} highly variable genes")

# Scale data (mean=0, variance=1) to prevent highly expressed genes from dominating
print("Scaling data...")
sc.pp.scale(adata, max_value=10)

# ═══════════════════════════════════════════════════════════════════════════
# DIMENSIONALITY REDUCTION
# ═══════════════════════════════════════════════════════════════════════════

# Run PCA on highly variable genes
print("Running PCA...")
sc.tl.pca(adata, 
          svd_solver='arpack',
          use_highly_variable=True)

# Save variance explained plot
sc.pl.pca_variance_ratio(adata, n_pcs=15, save='_salmon_variance.png', show=False)
print("PCA complete - saved variance plot")

# ═══════════════════════════════════════════════════════════════════════════
# CLUSTERING
# ═══════════════════════════════════════════════════════════════════════════

# Build k-nearest neighbor graph
print("Building neighbor graph...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)

# Leiden clustering - adjusted for smaller dataset
# Note: With only 16 cells, we expect fewer clusters than the paper's 11
print("Running Leiden clustering...")
sc.tl.leiden(adata, resolution=0.5)

# Renumber clusters starting from 1 instead of 0
adata.obs['leiden'] = (adata.obs['leiden'].astype(int) + 1).astype(str)

n_clusters = adata.obs['leiden'].nunique()
print(f"Found {n_clusters} clusters")
print(adata.obs['leiden'].value_counts().sort_index())

# ═══════════════════════════════════════════════════════════════════════════
# t-SNE VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════

# Run t-SNE with same perplexity as paper
# Note: Perplexity=5 is recommended for small datasets
sc.tl.tsne(adata,
           perplexity=5,
           use_rep='X_pca',
           random_state=42,
           n_jobs=1)

print("t-SNE complete")

# ═══════════════════════════════════════════════════════════════════════════
# CREATE VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════

# Extract data for plotting
tsne_coords = adata.obsm['X_tsne']
clusters = adata.obs['leiden'].astype(str)
unique_clusters = sorted(clusters.unique(), key=int)

# Create color map for clusters
colours = plt.cm.tab20(np.linspace(0, 1, len(unique_clusters)))
colour_map = dict(zip(unique_clusters, colours))

# Create figure
fig, ax = plt.subplots(figsize=(10, 8))

# Draw smooth background regions using Gaussian KDE
print("Drawing cluster backgrounds...")
for cluster in unique_clusters:
    mask = clusters == cluster
    if mask.sum() > 2:  # Need at least 3 points for KDE
        cluster_points = tsne_coords[mask]
        
        # Create grid for contour
        x_min, x_max = cluster_points[:, 0].min() - 5, cluster_points[:, 0].max() + 5
        y_min, y_max = cluster_points[:, 1].min() - 5, cluster_points[:, 1].max() + 5
        
        xx, yy = np.meshgrid(
            np.linspace(x_min, x_max, 100),
            np.linspace(y_min, y_max, 100)
        )
        
        # Calculate density
        try:
            positions = np.vstack([xx.ravel(), yy.ravel()])
            kernel = gaussian_kde(cluster_points.T)
            f = np.reshape(kernel(positions).T, xx.shape)
            
            # Draw filled contour
            ax.contourf(xx, yy, f, levels=[f.max()*0.4, f.max()], 
                       colors=[colour_map[cluster]], 
                       alpha=0.3, 
                       zorder=0)
        except:
            pass

# Plot individual points
print("Plotting cells...")
for cluster in unique_clusters:
    mask = clusters == cluster
    if mask.sum() > 0:
        ax.scatter(
            tsne_coords[mask, 0],
            tsne_coords[mask, 1],
            c=[colour_map[cluster]],
            marker='o',
            s=100,
            alpha=0.9,
            edgecolors='black',
            linewidths=0.5,
            label=f'Cluster {cluster}',
            zorder=2
        )

# Add cluster labels at centroids
for cluster in unique_clusters:
    mask = clusters == cluster
    if mask.sum() > 0:
        centroid_x = tsne_coords[mask, 0].mean()
        centroid_y = tsne_coords[mask, 1].mean()
        ax.text(centroid_x, centroid_y, f'c{cluster}',
                fontsize=12, fontweight='bold',
                ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', 
                         facecolor='white', 
                         alpha=0.8,
                         edgecolor='none'),
                zorder=3)

# Labels and title
ax.set_xlabel('tSNE 1', fontsize=12)
ax.set_ylabel('tSNE 2', fontsize=12)
ax.set_title('Salmon Pipeline: t-SNE Clustering\nFetal Neocortex 12wpc (n=16 cells)', 
             fontsize=14, fontweight='bold', pad=20)

# Legend
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), 
          frameon=True, fontsize=9)

plt.tight_layout()
plt.savefig('salmon_tsne_clustering.png', dpi=300, bbox_inches='tight')
plt.close()

print("Saved: salmon_tsne_clustering.png")

# ═══════════════════════════════════════════════════════════════════════════
# FIND MARKER GENES
#this part probably wont apply to the 16 cell alt pipeline recreation 
# ═══════════════════════════════════════════════════════════════════════════

print("\n=== Finding marker genes ===")

# Perform differential expression analysis
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# Create marker gene plot
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, 
                        save='_salmon_markers.png', show=False)
print("Saved: figures/rank_genes_groups_leiden_salmon_markers.png")

# Print top marker genes for each cluster
print("\nTop 5 marker genes per cluster:")
result = adata.uns['rank_genes_groups']
for cluster in unique_clusters:
    print(f"\nCluster {cluster}:")
    genes = result['names'][cluster][:5]
    for i, gene in enumerate(genes, 1):
        print(f"  {i}. {gene}")

# ═══════════════════════════════════════════════════════════════════════════
# SAVE RESULTS
# ═══════════════════════════════════════════════════════════════════════════

print("\n=== Saving results ===")

# Save cluster assignments
cluster_df = pd.DataFrame({
    'cell_id': adata.obs_names,
    'cluster': adata.obs['leiden'],
    'tissue': adata.obs['tissue'],
    'stage': adata.obs['stage']
})
cluster_df.to_csv('salmon_cluster_assignments.csv', index=False)
print("Saved: salmon_cluster_assignments.csv")

# Save AnnData object for future use
adata.write('salmon_analysis.h5ad')
print("Saved: salmon_analysis.h5ad")

# ═══════════════════════════════════════════════════════════════════════════
# SUMMARY STATISTICS
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print(f"Total cells analyzed: {adata.n_obs}")
print(f"Total genes analyzed: {adata.n_vars}")
print(f"Clusters identified: {n_clusters}")
print(f"\nCluster sizes:")
for cluster in unique_clusters:
    count = (adata.obs['leiden'] == cluster).sum()
    print(f"  Cluster {cluster}: {count} cells")

print("\nOutput files:")
print("  - salmon_tsne_clustering.png")
print("  - salmon_cluster_assignments.csv")
print("  - salmon_analysis.h5ad")
print("  - figures/pca_variance_ratio_salmon_variance.png")
print("  - figures/rank_genes_groups_leiden_salmon_markers.png")
print("="*70)

