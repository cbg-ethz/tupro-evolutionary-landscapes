import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.optimize import nnls
import seaborn as sns
import os
import sys

from utils_plots import plotMatrix


def getWMatrixOrCompute(path, X, Z):
  if os.path.isfile(path):
    print(path, "found ...")
    W = pd.read_csv(path, index_col=0)
  else:
    print("Computing W ...")
    W = recoverWMatrix(X, Z)
    W.to_csv(path)
  return W


def plotMatrices(X_df, Z_df, W_df, heatmap_filename_prefix, thresholds=None):
  threshold=0
  if thresholds and "X" in thresholds:
    threshold = thresholds["X"]
  plotMatrix(X_df, heatmap_filename_prefix + "_X", filtering_threshold=threshold)

  threshold=0
  if thresholds and "Z" in thresholds:
    threshold = thresholds["Z"]
  plotMatrix(Z_df, heatmap_filename_prefix + "_Z", filtering_threshold=threshold)

  threshold=0
  if thresholds and "W" in thresholds:
    threshold = thresholds["W"]
  plotMatrix(W_df, heatmap_filename_prefix + "_W", filtering_threshold=threshold)


def recoverWMatrix(X_df, Z_df):

  X = X_df.to_numpy()
  Z = Z_df.to_numpy()

  r = Z.shape[1]  # number of components
  n = X.shape[1]  # number of samples
  W = np.zeros((r, n))

  for i in range(n):
    W[:, i], _ = nnls(Z, X[:, i])

  X_reconstructed = Z @ W
  reconstruction_error = np.linalg.norm(X - X_reconstructed, ord='fro')
  original_norm = np.linalg.norm(X, ord='fro')
  print("Reconstructing W with error:", reconstruction_error / original_norm)

  W_df = pd.DataFrame(
      W,
      index = Z_df.columns,   # rows = components (same as Z columns)
      columns = X_df.columns) # columns = same as X columns 

  return W_df


def process_scDEF(h5ad_file, neutral_clones, out_dir):

  adata = ad.read_h5ad(h5ad_file)
  adata = adata[~adata.obs['batch'].isin(["OHAMUME-T", "OHACAHY-T"])].copy()

  adata.obs['subclone'] = adata.obs['batch'].astype(str).str.split('-').str[0] + '_' + adata.obs['scatrex_node'].astype(str)
  adata.obs['subclone'] = adata.obs['subclone'].astype('category')
  adata.obs['is_neutral'] = adata.obs['subclone'].isin(neutral_clones)

  # Plot raw expression profiles.
  '''
  sc.pp.neighbors(adata)
  sc.tl.umap(adata)
  #sc.pl.umap(adata, color="is_neutral", legend_loc='right margin', save="_allCells_neutral_or_not.png")
  sc.pl.umap(adata, color="batch", legend_loc='right margin', save="_allCells_raw.png")
  sys.exit()
  '''

  palette = sns.color_palette("husl", n_colors=192)
  category_to_color = dict(zip(adata.obs['subclone'].cat.categories, palette))
  samples = list(adata.obs.index)

  layers = ["factors", "hfactors", "hhfactors", "hhh"]
  subclone_signatures = {}
  nmf_matrices = {}
  for layer, obsm_name in enumerate(layers):
    '''
    sc.pp.neighbors(adata, use_rep="X_" + obsm_name)
    sc.tl.umap(adata)
    sc.pl.umap(
        adata, 
        color='subclone', 
        #title='Cell factors layer ' + str(layer),
        size=20, # point size
        alpha=0.5,
        frameon=False, 
        palette=category_to_color, 
        #legend_loc='right',
        save="_allCells_factors_layer" + str(layer) + ".png")
    '''

    df_factors = pd.DataFrame(adata.obsm["X_" + obsm_name], adata.obs.index)
    nmf_matrices["Z_" + str(layer)] = df_factors.copy()
    df_factors['subclone'] = adata.obs['subclone']
    subclone_signatures[layer] = df_factors.groupby('subclone', observed=False).mean()

  X = df = pd.DataFrame(
      adata.X,
      index=adata.obs_names,
      columns=adata.var_names
  )

  # X = Z_0 * W_0
  W_0 = getWMatrixOrCompute("W_0.csv", X, nmf_matrices["Z_0"])
  plotMatrices(X, nmf_matrices["Z_0"], W_0, os.path.join(out_dir, "nmf_layer0"), {"X":4}) 

  # Z_0 = Z_1 * W_1
  W_1 = getWMatrixOrCompute("W_1.csv", nmf_matrices["Z_0"], nmf_matrices["Z_1"])
  plotMatrices(nmf_matrices["Z_0"], nmf_matrices["Z_1"], W_1, os.path.join(out_dir, "nmf_layer1"))   

  # Z_1 = Z_2 * W_2
  W_2 = getWMatrixOrCompute("W_2.csv", nmf_matrices["Z_1"], nmf_matrices["Z_2"])
  plotMatrices(nmf_matrices["Z_1"], nmf_matrices["Z_2"], W_2, os.path.join(out_dir, "nmf_layer2"))         

  # Z_2 = Z_3 * W_3
  W_3 = getWMatrixOrCompute("W_3.csv", nmf_matrices["Z_2"], nmf_matrices["Z_3"])
  plotMatrices(nmf_matrices["Z_2"], nmf_matrices["Z_3"], W_3, os.path.join(out_dir, "nmf_layer3"))

  # Z_3 = Z_4 * W_4
  #W_4 = getWMatrixOrCompute("W_4.csv", nmf_matrices["Z_3"], nmf_matrices["Z_4"])
  #plotMatrices(nmf_matrices["Z_3"], nmf_matrices["Z_4"], W_4, os.path.join(out_dir, "nmf_layer4"))

  # Top factors per cluster
  #  top_factors_by_cluster = mean_factor_by_cluster.apply(
  #      lambda row: row.sort_values(ascending=False).index[:2], axis=1
  #  )

  # Top genes for each factor
  # factor_i = 0
  # top_gene_indices = W[factor_i].argsort()[::-1][:10]
  # top_genes = gene_names[top_gene_indices]

  # Top genes for those dominant factors
  '''
  for cluster_id in top_factors_by_cluster.index:
    print(f"\nCluster {cluster_id}:")

    for factor in top_factors_by_cluster.loc[cluster_id]:
        factor = int(factor)  # ensure integer index
        top_genes_idx = W[factor].argsort()[::-1][:10]  # Top 10 genes
        top_genes = [gene_names[i] for i in top_genes_idx]
        print(f"  Factor {factor} → Top genes: {top_genes}")  
  '''

  return subclone_signatures
    
    

