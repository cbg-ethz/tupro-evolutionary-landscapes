import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import gcf
import seaborn as sns

from parse_metadata_tupro_ovarian import *

def plotHierarchicalClustering(
        filename,
        df,
        do_clustering = True,
        precomputed_clustering = None,
        clusters = None,
        metric="euclidean",
        method="ward",
        cmap="coolwarm",
        annot=False,
        font_scale=0.6,
        metadata=None):

    print("Generating heatmap:", filename + ".pdf ...")

    if metadata:
      sample_label_colors, color_codes = parse_metadata_tupro_ovarian(df.index, metadata)
    else:
      sample_label_colors = None

    # Add color for sample name
    sample_names = df.index.str.split("_").str[0].unique().tolist()
    num_samples = len(sample_names)
    cmap_samples = plt.cm.get_cmap("tab20", num_samples)  
    palette = [cmap_samples(i) for i in range(num_samples)]
    sample_color_dict = dict(zip(sample_names, palette))
    sample_label_colors.index.name = "Sample"
    sample_label_colors['Sample'] = sample_label_colors.index.str.split('_').str[0].map(sample_color_dict)

    sns.set(font="Arial", font_scale=font_scale)
    plot = sns.clustermap(
        df,
        row_cluster=do_clustering,
        col_cluster=do_clustering,
        row_colors=sample_label_colors,
        cmap=cmap,
        annot=annot,
        fmt=".2f",
        annot_kws={"size": 5},
        row_linkage=precomputed_clustering,
        col_linkage=precomputed_clustering,
        xticklabels=True,
        yticklabels=True
    )

    # Adjust axis label text size.
    ax = plot.ax_heatmap
    bbox = ax.get_window_extent()
    width_in = bbox.width 
    font_size = width_in / len(df.index) - 0.5
    ax.tick_params(axis='both', labelsize=font_size)    

    if clusters:
      hline_index = []
      cnt = 0
      for idx in range(len(clusters) - 1):
        cnt = cnt + len(clusters[idx])
        hline_index.append(cnt)
      ax_heatmap = plot.ax_heatmap
      ax_heatmap.hlines(hline_index, *ax_heatmap.get_xlim(), color='whitesmoke', linewidth=1)

    if metadata: 
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
    plot.savefig(filename + ".png", format='png', dpi=300, bbox_inches='tight')

def plotMatrix(df, file_path, filtering_threshold=0):

  print("Generating heatmap of size", df.shape, ":", file_path)

  print(df.min().min(), df.max().max())
  flat_vals = df.values.flatten()

  log_df = np.log1p(df)

  plt.clf()
  plt.hist(log_df.values.flatten(), color='steelblue', edgecolor='black')
  plt.xlabel("Value")
  plt.ylabel("Frequency")
  plt.title("Histogram of All Values in DataFrame")
  plt.savefig(file_path + "_histo.png", format='png', dpi=150, bbox_inches='tight')

  rows_to_keep = (log_df > filtering_threshold).any(axis=1)
  cols_to_keep = (log_df > filtering_threshold).any(axis=0)
  filtered_df = log_df.loc[rows_to_keep, cols_to_keep]

  plt.clf()
  sns.heatmap(
      filtered_df, 
      #mask=mask,
      cmap="viridis",
      vmin=0, 
  )
  plt.savefig(file_path + ".png", format='png', dpi=150, bbox_inches='tight')  
  

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
  histogram_fig.savefig(output_path + ".pdf", format='pdf', dpi=300, bbox_inches='tight')


def plotPieChartsClusterHeterogeneity(clusters, clone_sizes, file_path):

  samples = list(set([item.split("_")[0] for sublist in clusters for item in sublist]))
  num_samples = len(samples)
  cols = 5 #math.ceil(math.sqrt(len(samples)))  
  rows = math.ceil(num_samples / cols)

  fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
  axes = axes.flatten()
  colors = sns.color_palette("tab10")
  for i, sample in enumerate(samples):
    values = []
    for cluster in clusters:
      subclones = [subclone for subclone in cluster if subclone.startswith(sample)]
      if len(subclones):
        values.append(sum(clone_sizes[k] if k in clone_sizes else 0  for k in subclones))
      else:
        values.append(0)
    values = [x for x in values if x != 0] # remove the 0s
    axes[i].pie(values, colors=colors[:len(values)])#, labels=values.index)
    axes[i].text(0, -1.6, sample, fontsize=50, ha='center')

  # Hide any unused subplots
  for j in range(i + 1, len(axes)):
    axes[j].axis('off')
  plt.tight_layout()
  plt.savefig(file_path + ".pdf", format='pdf', dpi=300, bbox_inches='tight')

