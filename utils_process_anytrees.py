import pandas as pd
import numpy as np
import copy
import os
import sys
from anytree import PreOrderIter, RenderTree
from scipy.spatial import distance
from scipy.cluster import hierarchy

from utils_clustering import getSimilarityClusters
from utils_anytree import removeTreeAttributes
from utils_oncotree2vec import saveOncotree2vecInput
from utils_plots import plotHierarchicalClustering 
from utils_plots import generateStackedHistogram
from utils_plots import plotPieChartsClusterHeterogeneity

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
          df_jaccard_distances.loc[node_1, node_2] = dist
          df_jaccard_distances.loc[node_2, node_1] = dist
    check_distance_dataframe_integrity(df_jaccard_distances)
    dataframes[idx] = df_jaccard_distances
  return dataframes


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

def processCNTrees(anytrees_list, marker_genes, df_distances_csv, out_dir, metadata):
 
  print("Processing CN trees...")

  # Combine gene sets. 
  if os.path.isfile(df_distances_csv):
    df_combined_distances = pd.read_csv(df_distances_csv, index_col=0)  
    print("Found", df_distances_csv)
  else:
    dataframes = getSubcloneDistances(anytrees_list, marker_genes)
    df_combined_distances = pd.concat(dataframes.values()).groupby(level=0).min()
    df_combined_distances.to_csv(df_distances_csv)
    print("Saved ", df_distances_csv)

  # Process trees in oncotree2vec input format. 
  cn_dir = os.path.join(out_dir, "cn")
  os.makedirs(cn_dir, exist_ok=True)
  oncotree2vec_dir = os.path.join(cn_dir, "oncotree2vec_cnTrees")
  os.makedirs(oncotree2vec_dir, exist_ok=True)

  empty_trees = {}
  for key, tree in anytrees_list[-1].items():
    empty_trees[key] = removeTreeAttributes(tree)

  clustering = hierarchy.linkage(distance.pdist(df_combined_distances), metric="euclidean", method="ward")
  all_cluster_sizes = []
  for distance_threshold in np.linspace(0, 0.9, num=10):
    distance_threshold = round(distance_threshold,1)
    clone_clusters = getSimilarityClusters(clustering, df_combined_distances, distance_threshold)

    anytrees = copy.deepcopy(empty_trees)
    setMatchingLabels(anytrees, clone_clusters)
    if "ODAFUGU" in anytrees:
      test_tree = anytrees["ODAFUGU"]
      assert test_tree.root.matching_label == 0
      print(PreOrderIter(test_tree))
      assert next((n for n in PreOrderIter(test_tree) if getattr(n, "node_id", None) == -100), None).matching_label == 1
      assert next((n for n in PreOrderIter(test_tree) if getattr(n, "node_id", None) == 103), None).matching_label == 1
      assert next((n for n in PreOrderIter(test_tree) if getattr(n, "node_id", None) == 109), None).matching_label == 1
      assert next((n for n in PreOrderIter(test_tree) if getattr(n, "node_id", None) == 116), None).matching_label == 1
      assert next((n for n in PreOrderIter(test_tree) if getattr(n, "node_id", None) == 96), None).matching_label > 1
      assert next((n for n in PreOrderIter(test_tree) if getattr(n, "node_id", None) == 114), None).matching_label > 1
    oncotree2vec_subdir = os.path.join(oncotree2vec_dir, str(distance_threshold))
    os.makedirs(oncotree2vec_subdir)
    saveOncotree2vecInput(anytrees, oncotree2vec_subdir, duplicate_trees=True)

    clone_cluster_sizes = [len(cluster) for cluster in clone_clusters]
    for item in clone_cluster_sizes:
      all_cluster_sizes.append({"threshold":distance_threshold, "size":item})

  figure_dir = os.path.join(cn_dir, "figures")
  os.makedirs(figure_dir, exist_ok=True)
  generateStackedHistogram(pd.DataFrame.from_dict(all_cluster_sizes), os.path.join(figure_dir, "subclone_sizes"))


def processRNATrees(scatrex_trees, rna_signatures, scatrex_cn_nodeid_map, out_dir, metadata):

  print("Building RNA signature trees...")

  rna_dir = os.path.join(out_dir, "rna")
  os.makedirs(rna_dir, exist_ok=True) 
  oncotree2vec_dir = os.path.join(rna_dir, "oncotree2vec_rnaTrees")
  os.makedirs(oncotree2vec_dir, exist_ok=True)
  figures_dir = os.path.join(rna_dir, "figures")
  os.makedirs(figures_dir, exist_ok=True)

  empty_trees = {}
  clone_sizes = {}
  for key, tree in scatrex_trees.items():
    for node in PreOrderIter(tree):
      node_name = "_".join([key, node.node_id])
      clone_sizes[node_name] = node.num_cells 
    empty_trees[key] = removeTreeAttributes(tree) 

  metric = "correlation"
  distance_matrices = {}
  for key, df in rna_signatures.items():
    pairwise_distances = distance.pdist(df, metric=metric)
    df_distances = pd.DataFrame(distance.squareform(pairwise_distances))
    df_distances.index = df.index
    df_distances.columns = df.index
    distance_matrices[key] = df_distances

    clustering = hierarchy.linkage(pairwise_distances, metric="euclidean", method="ward")
    clone_clusters = getSimilarityClusters(clustering, df_distances, 0.75)

    plotHierarchicalClustering(
        df = df_distances,
        do_clustering = True,
        precomputed_clustering = clustering,
        clusters = clone_clusters,
        cmap="YlOrBr_r",
        metadata=metadata,
        filename=os.path.join(figures_dir, "subclone_distances_" + str(key)),      
    )

    plotPieChartsClusterHeterogeneity(
      clone_clusters, 
      clone_sizes,
      os.path.join(figures_dir, "pie_chart_" + str(key)), 
    )

    anytrees = copy.deepcopy(empty_trees)
    setMatchingLabels(anytrees, clone_clusters)
    if key == 1 and "OBABACY" in anytrees:
      test_tree = anytrees["OBABACY"]
      assert test_tree.root.matching_label == 0
      assert next((n for n in PreOrderIter(test_tree) if getattr(n, "node_id", None) == "B"), None).matching_label > 1
      assert next((n for n in PreOrderIter(test_tree) if getattr(n, "node_id", None) == "C"), None).matching_label == 1 # missing from scDEF
      label_node_1 = next((n for n in PreOrderIter(test_tree) if getattr(n, "node_id", None) == "D"), None).matching_label
      label_node_2 = next((n for n in PreOrderIter(test_tree) if getattr(n, "node_id", None) == "E"), None).matching_label
      label_node_3 = next((n for n in PreOrderIter(test_tree) if getattr(n, "node_id", None) == "F"), None).matching_label
      assert label_node_1 > 1
      assert label_node_2 > 1
      assert label_node_3 > 1
      assert len({label_node_1, label_node_2, label_node_3}) == 3 # values are not equal to each other

    oncotree2vec_subdir = os.path.join(oncotree2vec_dir, str(key))
    os.makedirs(oncotree2vec_subdir)
    # Set SCATrEx the node ids to the corresponding SCICoNE ones, because ocnotree2vec wants numerical ids.
    for key, tree in anytrees.items():
      for node in PreOrderIter(tree):
        node.node_id = scatrex_cn_nodeid_map[key][node.node_id]
    saveOncotree2vecInput(anytrees, oncotree2vec_subdir, duplicate_trees=True)

