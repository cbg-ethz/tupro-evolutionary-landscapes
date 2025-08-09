# Usage: python main.py tupro-ovarian 

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import json
import re
import os,sys
import time
import sys
import statistics

from anytree.importer import JsonImporter
from anytree import RenderTree

from utils import *
from utils_data import *
from utils_clone_tree import getCloneTrees
from utils_oncotreevis import writeOncotreeVISInput
from utils_process_anytrees import processCNTrees
from utils_scatrex_trees import getCorrespondingSCATrExNode
from utils_process_anytrees import processRNATrees
from utils_process_scdef import process_scDEF

datasets = {
  "tupro-melanoma": {
      "scicone_trees_priority": "data/samples_melanoma_v1.15_priority_genes_26032024.js",
      "scicone_trees_all": "data/samples_melanoma_v1.15_all_genes_26032024.js",
      "scatrex_trees": "data/1754242433_scatrex_melanoma_18june2024_-2:2.json",
      "neutral_clones": neutral_clones_melanoma,
      "metadata": "data/metadata_melanoma.csv",
      "samples_to_remove": samples_to_remove_melanoma,
      "dna_marker_genes": dna_marker_genes_melanoma
  },
  "tupro-ovarian": {
      "scicone_trees_priority": "data/samples_ovarian_v1.15_priority_genes_26032024.js",
      "scicone_trees_all": "data/samples_ovarian_v1.15_all_genes_26032024.js",
      "scicone_df_distances_event": "data/tupro-ovarian_event_df_distances.csv",
      "scicone_df_distances_clone": "data/tupro-ovarian_clone_df_distances.csv",
      "scatrex_trees": "data/out/1754242433_scatrex_ovarian_14nov2024b_tota_exp_1:1.json",
      "scDEF": "data/scDEF-data/hgsoc_scatrex_scdef.h5ad",
      "neutral_clones": neutral_clones_ovarian,
      "metadata": "data/metadata_ovarian.csv",
      "samples_to_remove": samples_to_remove_ovarian,
      "dna_marker_genes": dna_marker_genes_ovarian
  },
}

def removeGenes(anytrees, gene_list):
  for sample, anytree in anytrees.items():
    for node in PreOrderIter(anytree):
      if hasattr(node, 'gene_cn_events'):
        genes_to_remove = set(node.gene_cn_events.keys()).intersection(gene_list)
        for gene in genes_to_remove:
          del node.gene_cn_events[gene]
  return anytrees

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("dataset", choices=['tupro-melanoma', 'tupro-ovarian'])
  parser.add_argument("-o", "--output_dir", default="out")
  args = parser.parse_args()

  dataset_name = args.dataset
  if dataset_name not in datasets:
    sys.exit('Dataset not found')

  os.makedirs(args.output_dir, exist_ok=True)

  timestamp = str(int(time.time()))
  filename_prefix = "_".join([timestamp, args.dataset])
  out_dir = os.path.join(args.output_dir, filename_prefix)

  # CN trees.
  data = datasets[dataset_name]
  dna_marker_genes = data["dna_marker_genes"]
  neutral_clones = data["neutral_clones"]
  samples_to_remove = data["samples_to_remove"]
  metadata = readMetadata(data["metadata"])
  gene_chr_map = getGeneToChrMap("data/gene_chr_map.csv", "index")

  scicone_event_trees_priority_genes = readJsonTrees(data["scicone_trees_priority"], metadata, samples_to_remove, neutral_clones)
  scicone_event_trees_all_genes = readJsonTrees(data["scicone_trees_all"], metadata, samples_to_remove, neutral_clones)
  if "ovarian" in dataset_name:
    genes_chr_y = [gene for gene in gene_chr_map if gene_chr_map[gene]=="Y"]
    scicone_event_trees_priority_genes = removeGenes(scicone_event_trees_priority_genes, genes_chr_y)
    scicone_event_trees_all_genes = removeGenes(scicone_event_trees_all_genes, genes_chr_y)
  scicone_clone_trees_priority_genes = getCloneTrees(scicone_event_trees_priority_genes)
  scicone_clone_trees_all_genes = getCloneTrees(scicone_event_trees_all_genes)

  # SCATrEx trees.
  scatrex_clone_trees = readJsonTrees(data["scatrex_trees"], metadata, samples_to_remove)
  # Remove the root.
  for key, tree in scatrex_clone_trees.items():
    assert len(tree.children) == 1
    next_node = copy.deepcopy(tree.children[0])
    next_node.parent = None
    scatrex_clone_trees[key] = next_node
  # Remove broken samples
  del scatrex_clone_trees["OHACAHY"]
  del scatrex_clone_trees["OHAMUME"]

  printDatasetStatistics(scicone_event_trees_all_genes, "SCICoNE trees")
  printDatasetStatistics(scatrex_clone_trees, "SCATrEx trees")
  print()
  #printMetadataStatistics(metadata, scicone_event_trees_all_genes.keys())

  # Combine datasets for oncotreeVIS input. 
  '''
  combined_datasets = {}
  for key, tree in scatrex_clone_trees.items():
    combined_datasets[key + "_C_SCATrEx"] = tree
  for key, tree in scicone_clone_trees_all_genes.items():
    combined_datasets[key + "_B_SCICoNE_clone"] = tree
  for key, tree in scicone_event_trees_all_genes.items():
    combined_datasets[key] = tree + "_A_SCICoNE_event"] = tree
  combined_datasets = dict(sorted(combined_datasets.items()))
  oncotreeVIS_dir = os.path.join(out_dir, "oncotreeVIS", )
  os.makedirs(oncotreeVIS_dir)
  filename = "_".join(["scicone", "and", "scatrex", "trees.json"])
  writeOncotreeVISInput(
      combined_datasets, 
      metadata, 
      os.path.join(oncotreeVIS_dir, filename),
  )
  '''

  # Create oncotree2vec input for CN trees.
  '''
  processCNTrees(
      anytrees_list = [scicone_event_trees_priority_genes, scicone_event_trees_all_genes], 
      marker_genes = dna_marker_genes,
      df_distances_csv = "df_cn_subclone_distances_" + args.dataset + ".csv",
      out_dir = out_dir,
      metadata = metadata,
  )
  '''

  # Find the correspondance between the SCICoNE and SCATrEx nodes.
  # Add node num_cells to SCATrEx trees.
  scatrex_scicone_node_correspondance = {} 
  for sample_name, scicone_tree in scicone_clone_trees_all_genes.items():
    scatrex_scicone_node_correspondance[sample_name] = {}
    if sample_name in scatrex_clone_trees:
      # TODO: assert matchTreeNodes(scicone_tree, scatrex_clone_trees[sample_name])
      for node_scicone in PreOrderIter(scicone_tree):
        flag = False
        node_scatrex = getCorrespondingSCATrExNode(node_scicone, scatrex_clone_trees[sample_name])
        node_scatrex.num_cells = node_scicone.num_cells
        node_scatrex_label = node_scatrex.node_id
        assert node_scatrex is not None
        node_scicone.node_id_scatrex = node_scatrex_label
        scatrex_scicone_node_correspondance[sample_name][node_scatrex_label] = node_scicone.node_id

  # Get scDEF subclone RNA profiles for different levels in the hierarchy.
  neutral_clones_scatrex = []
  for sample_name, root in scicone_clone_trees_all_genes.items():
    if sample_name in scatrex_clone_trees:
      neutral_clones_scatrex.extend("_".join([sample_name, node.node_id_scatrex]) for node in PreOrderIter(root) if node.is_neutral)  

  scdef_dir = os.path.join(out_dir, "scDEF")
  os.makedirs(scdef_dir, exist_ok=True)
  signatures_different_resolutions = process_scDEF(data["scDEF"], neutral_clones_scatrex, scdef_dir)

  processRNATrees(
    scatrex_clone_trees, 
    signatures_different_resolutions, 
    scatrex_scicone_node_correspondance,
    out_dir,
    metadata,
  )






