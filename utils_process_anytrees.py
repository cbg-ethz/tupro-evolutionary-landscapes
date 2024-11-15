import pandas as pd
import numpy as np
import copy
import os
from anytree import PreOrderIter
from scipy.spatial import distance
from scipy.cluster import hierarchy

from utils_clustering import getSimilarityClusters
from utils_heatmaps import generateHeatmaps
from utils_heatmaps import generateStackedHistogram
from utils_oncotree2vec_input import saveOncotree2vecInput 
from utils_oncotreevis import writeOncotreeVISInput
from utils import create_dir
from utils_anytree import isRoot
from parse_metadata_tupro_ovarian import *

def jaccard_distance(set_1, set_2, is_malignant=True):
  # By default the jaccard distance between two empty sets is 0. If the nodes are malignant, then we want to return distance 1.
  if len(set_1) == 0 and len(set_2) == 0 and is_malignant:
    return 1
  intersection = len(set_1.intersection(set_2))
  union = len(set_1.union(set_2))
  if union > 0:
    iou = intersection / union
  else: # union size is 0, i.e., no genes are affected whatssoever (i.e., neutral clone/event node)
    iou = 1 # similarity one betweeb neutral clones / event nodes.
  return 1-iou

def isRoot(node):
  return node.parent is None

def getAffectedGenes(anytrees):
  genes = set()
  for key, tree in anytrees.items():
    for node in PreOrderIter(tree):
      if hasattr(node, 'gene_cn_events'):
        genes.update(node.gene_cn_events.keys())
  return list(genes)

def geneCNToString(gene, cn, neutral_state):
  gene = gene.split(".")[0]
  if cn > 0:
    return gene
  else:
    if cn == -neutral_state:
      return gene + "_lost"
    else:
      return gene + "_del"

def getSubcloneDistances(anytrees_list, marker_genes):
  def check_distance_dataframe_integrity(df):
    # Assert the matrix is symmetric.
    assert np.all(np.abs(df - df.T) < 0.0001)
    # Assert diagonal values are 0.
    assert (np.diag(df) == [0] * len(df.index)).all()

  # Populate malignant_node_gene_map.
  malignant_node_gene_map = {}
  for idx, anytrees in enumerate(anytrees_list):
    for sample_name, anytree in anytrees.items():
      for node in PreOrderIter(anytree):
        if node.is_neutral:
          continue
        assert hasattr(node,"gene_cn_events")
        node_key = "_".join([sample_name, str(node.node_id)])
        if node_key not in malignant_node_gene_map:
          malignant_node_gene_map[node_key] = {}
        gene_set = set([geneCNToString(gene, node.gene_cn_events[gene], node.cn_neutral_state) for gene in node.gene_cn_events])
        malignant_node_gene_map[node_key][idx] = gene_set

  # Add marker genes.
  marker_genes_deleted = [gene + "_del" for gene in marker_genes]
  marker_genes_lost = [gene + "_lost" for gene in marker_genes]
  marker_genes = marker_genes + marker_genes_deleted + marker_genes_lost
  idx = len(anytrees_list) # next index
  for node, gene_set_map in malignant_node_gene_map.items():
    gene_set_map[idx] = gene_set_map[idx-1].intersection(set(marker_genes))

  # Combine distances between gene sets.
  nodes = list(malignant_node_gene_map.keys())
  dataset_idxs = list(range(idx+1))
  dataframes = {}
  for idx in dataset_idxs:
    df_jaccard_distances = pd.DataFrame(0, columns=nodes, index=nodes).astype(float)
    for node_1 in nodes:
      for node_2 in nodes:
        if node_1 < node_2:
          dist = jaccard_distance(malignant_node_gene_map[node_1][idx], malignant_node_gene_map[node_2][idx])
          df_jaccard_distances[node_1][node_2] = dist
          df_jaccard_distances[node_2][node_1] = dist
    check_distance_dataframe_integrity(df_jaccard_distances)
    dataframes[idx] = df_jaccard_distances
  return dataframes

'''
def getSubcloneClusters(df_distances, distance_threshold, plot_histogram=None, plotting_dir_prefix=None): 
  # Cluster the clones.
  clustering = hierarchy.linkage(
    distance.pdist(df_distances), metric="euclidean", method="ward")
  nodes = df_distances.index
  return getSimilarityClusters(clustering, nodes, df_distances, distance_threshold)

  if plotting_dir_prefix:
    # Plot the big pair-wise clone distance heatmap.
    sample_label_colors, color_codes = parse_metadata_tupro_ovarian(df_distances.index, metadata)
    generate_heatmaps(df_distances, clustering, clone_clusters, sample_label_colors, color_codes, plotting_dir_prefix,
        plot_histogram=True, plot_cluster_heatmaps=False)
'''

def setMatchingLabels(anytrees, clone_clusters):
  malignant_node_matching_labels = {}
  clone_cluster_id = 2 # start the labeling of the malignant clones from 2 (0 is reserved for the root and 1 for the neutral clones).
  for cluster in clone_clusters:
    for clone in cluster: 
      malignant_node_matching_labels[clone] = clone_cluster_id
    clone_cluster_id += 1

  for sample_name, tree in anytrees.items():
    for node in PreOrderIter(tree):
      node_key = "_".join([sample_name, str(node.node_id)])
      if isRoot(node): 
        node.matching_label = 0 # root       
      elif node_key in malignant_node_matching_labels:
        assert malignant_node_matching_labels[node_key] >= 2
        node.matching_label = malignant_node_matching_labels[node_key]
      else:
        node.matching_label = 1 # neutral or aux nodes

def processTrees(anytrees_list, dna_marker_genes, highlighted_genes, metadata, out_dir_prefix, df_distances_csv=None):
 
  print("Processing trees", out_dir_prefix, "...")

  create_dir(out_dir_prefix)
  plotting_dir = os.path.join(out_dir_prefix, "plots")
  create_dir(plotting_dir)
  filename_prefix = os.path.basename(out_dir_prefix)

  if df_distances_csv:
    df_combined_distances = pd.read_csv(df_distances_csv, index_col=0)  
  else:
    dataframes = getSubcloneDistances(anytrees_list, dna_marker_genes)
    df_combined_distances = pd.concat(list(dataframes.values())).min(level=0) #groupby(level=0).min()
    df_combined_distances.to_csv(os.path.join(out_dir_prefix, "_".join([filename_prefix, "df_distances.csv"])))  
  
  reference_anytrees = anytrees_list[-1]

  clustering = hierarchy.linkage(distance.pdist(df_combined_distances), metric="euclidean", method="ward")
  all_cluster_sizes = []
  for distance_threshold in np.linspace(0, 0.9, num=10):
    distance_threshold = round(distance_threshold,1)
    clone_clusters = getSimilarityClusters(clustering, df_combined_distances, distance_threshold)
    clone_cluster_sizes = [len(cluster) for cluster in clone_clusters]
    for item in clone_cluster_sizes:
      all_cluster_sizes.append({"threshold":distance_threshold, "size":item})

    anytrees = copy.deepcopy(reference_anytrees)
    setMatchingLabels(anytrees, clone_clusters)
    oncotree2vec_dir = os.path.join(out_dir_prefix, "_".join([filename_prefix, "oncotree2vec_input", str(distance_threshold)]))
    saveOncotree2vecInput(anytrees, oncotree2vec_dir, duplicate_trees=True)

    oncotreeVIS_path = os.path.join(out_dir_prefix, filename_prefix + "_".join(["oncotreeVIS", str(distance_threshold)]) + ".json")
    oncotreeVIS_path_oncotree2vec = os.path.join(oncotree2vec_dir, "trees.json")
    writeOncotreeVISInput(anytrees, metadata, [oncotreeVIS_path, oncotreeVIS_path_oncotree2vec], "", highlighted_genes)

    if distance_threshold == 0.8:
      sample_label_colors, color_codes = parse_metadata_tupro_ovarian(df_combined_distances.index, metadata)
      generateHeatmaps(df_combined_distances, clustering, clone_clusters, sample_label_colors, color_codes, plotting_dir + "/" + filename_prefix + "_subclone_dist",
          plot_histogram=True, plot_cluster_heatmaps=True)

  generateStackedHistogram(pd.DataFrame.from_dict(all_cluster_sizes), plotting_dir + "/" + filename_prefix + "_histo_stacked_subclone_sizes.png")
  return 
  




