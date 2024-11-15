import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import gcf
import seaborn as sns

def plot_hierarchical_clustering(
        filename_base,
        df,
        show_hierarchical_clustering = True,
        precomputed_clustering = None,
        clusters = None,
        metric="euclidean",
        method="ward",
        sample_label_colors=None,
        color_codes=None,
        cmap="coolwarm",
        annot=True,
        vmin = 0,
        vmax = 1,
        font_scale=0.6):

    sns.set(font="Arial", font_scale=font_scale)
    plot = sns.clustermap(
        df,
        row_cluster=show_hierarchical_clustering,
        col_cluster=show_hierarchical_clustering,
        row_colors=sample_label_colors,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        annot=annot,
        fmt=".2f",
        annot_kws={"size": 5},
        row_linkage=precomputed_clustering,
        col_linkage=precomputed_clustering,
        xticklabels=True,
        yticklabels=True)

    plot.cax.set_visible(False)
    figure = plt.gcf()
    figure.set_size_inches(15, 9)
    plot.savefig(filename_base + ".png", format='png', dpi=900)

    if clusters:
      hline_index = []
      cnt = 0
      for idx in range(len(clusters) - 1):
        cnt = cnt + len(clusters[idx])
        hline_index.append(cnt)
      ax_heatmap = plot.ax_heatmap
      ax_heatmap.hlines(hline_index, *ax_heatmap.get_xlim(), color='whitesmoke', linewidth=0.5)

    if color_codes: 
      sns.set(font_scale=0.7)
      for label, color_code in color_codes.items():
        plot.ax_col_dendrogram.bar(0, 0, color="white", label=label, linewidth=0)
        for key, color in color_code.items():
          plot.ax_col_dendrogram.bar(0, 0, color=color, label=key, linewidth=0)
        plot.ax_col_dendrogram.bar(0, 0, color="white", label="", linewidth=0)
  
      ncol = 4 
      legend_box_position = (0.5, 1.15)
      l = plot.ax_col_dendrogram.legend(title="", loc="center", ncol=ncol, bbox_to_anchor=legend_box_position,
          bbox_transform=gcf().transFigure, facecolor='white', framealpha=1)

    plot.cax.set_visible(False)
    figure = plt.gcf()
    figure.set_size_inches(15, 9)
    plot.savefig("_".join([filename_base, "clusters"]) + ".png", format='png', dpi=900)

def generateHeatmaps(
    df_jaccard_distances,
    precomputed_clustering, 
    cluster_list, 
    sample_label_colors,
    color_codes,
    out_dir_prefix,
    plot_histogram=False,
    plot_cluster_heatmaps=False):

  run_clustering = True
  if type(precomputed_clustering) == type(None):
    run_clustering = False
  plot_hierarchical_clustering(
          out_dir_prefix + "_heatmap",
          df_jaccard_distances,
          show_hierarchical_clustering = run_clustering,
          precomputed_clustering = precomputed_clustering,
          clusters=cluster_list,
          sample_label_colors=sample_label_colors,
          color_codes=color_codes,
          cmap="YlOrBr_r",
          annot=False,
          font_scale=0.1)

  if plot_histogram:
    # Get the values from below the diagonal and compute the histogram.
    df_jaccard_distances_below_diag = df_jaccard_distances.mask(np.triu(np.ones(df_jaccard_distances.shape, dtype=np.bool_)))
    values = df_jaccard_distances_below_diag.values.tolist()
    values = [item for sublist in values for item in sublist] # flatten list of lists
    values = [x for x in values if str(x) != 'nan'] # remove nan value (last element)

    plt.figure()
    sns.set(font_scale=1)
    histogram = sns.histplot(data = values, bins=10)
    histogram.set_yscale('log')
    histogram.set(ylabel='Counts', xlabel='Jaccard distances')
    histogram_fig = histogram.get_figure()
    histogram_filename = "_".join([out_dir_prefix, "histo"]) + ".png"
    histogram_fig.savefig(histogram_filename, format='png', dpi=300)

  if plot_cluster_heatmaps:
    # Plot the annotated heatmap for sub-matrices of size 40x40, with overlap 5x5.
    num_rows = df_jaccard_distances.shape[0]
    cnt = -1
    for cluster in cluster_list:
      cnt = cnt + 1
      submatrix = df_jaccard_distances.loc[cluster, cluster]
      # TODO: assert submatrix.max().max() <= distance_threshold
      if submatrix.shape[0] == 1 and submatrix.values[0] == 0:
        continue
      plot_hierarchical_clustering(
            "_".join([out_dir_prefix, "cluster_" + str(cnt)]),
            submatrix,
            cmap="YlOrBr_r",
            annot=True,
            font_scale=1)
      # TODO assert submatrix.max().max() <= distance_threshold

def generateStackedHistogram(df, output_path, x_column="size", hue_column="threshold"):
  plt.figure()
  histogram = sns.histplot(df, x=x_column, hue=hue_column, alpha=1, discrete=True) 
  histogram.set_yscale('log')
  histogram.set(ylabel='Counts', xlabel='cluster sizes')
  max_value = df[x_column].max() + 1
  histogram.set_xticks(list(range(0, max_value, 1)))

  histogram.set_yscale('log')
  histogram.set(ylabel='Counts', xlabel='cluster sizes')
  histogram_fig = histogram.get_figure()
  histogram_fig.savefig(output_path, format='png', dpi=300)

