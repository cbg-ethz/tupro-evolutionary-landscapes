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

from utils_data import *
from utils_clone_tree import getCloneTrees
from utils_clone_tree import getChrTrees
from utils_process_anytrees import getAffectedGenes
from utils_process_anytrees import processTrees
from utils_clustering import getSimilarityClusters
#from utils_scatrex_trees import parse_scatrex_trees
from utils_scatrex_trees import matchSCATrExClones
from utils_oncotreevis import writeOncotreeVISInput
from utils import *

#from utils_tupro_input import *
#from utils_visualization import *
#from utils_graph2vec_input import *

_OUT_DIR = "out"

datasets = {
  "tupro-melanoma": {
      "scicone_trees_priority": "data/samples_melanoma_v1.15_priority_genes_26032024.js",
      "scicone_trees_all": "data/samples_melanoma_v1.15_all_genes_26032024.js",
      "scatrex_trees": "data/1729089878_scatrex_melanoma_18june2024_-2:2.json",
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
      #"scatrex_trees": "data/out/1731483421_scatrex_ovarian_12nov2024_tota_exp_1:1.json",
      #"scatrex_trees": "data/out/1731589422_scatrex_ovarian_14nov2024b_tota_exp_1:1.json",
      "scatrex_trees": "data/out/1731648452_scatrex_ovarian_14nov2024b_tota_exp_1:1.json",
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
  args = parser.parse_args()

  dataset_name = args.dataset

  if dataset_name not in datasets:
    sys.exit('Dataset not found')

  timestamp = str(int(time.time()))
  filename_prefix = "_".join([timestamp, args.dataset])
  _OUT_DIR = os.path.join(_OUT_DIR, filename_prefix)

  # CN trees.
  data = datasets["tupro-ovarian"]
  dna_marker_genes = data["dna_marker_genes"]
  neutral_clones = data["neutral_clones"]
  samples_to_remove = data["samples_to_remove"]
  metadata = readMetadata(data["metadata"])
  gene_chr_map = getGeneToChrMap("data/gene_chr_map.csv", "index")

  scicone_trees_event_priority_genes = readJsonTrees(data["scicone_trees_priority"], metadata, samples_to_remove, neutral_clones)
  scicone_trees_event_all_genes = readJsonTrees(data["scicone_trees_all"], metadata, samples_to_remove, neutral_clones)
  genes_chr_y = [gene for gene in gene_chr_map if gene_chr_map[gene]=="Y"]
  scicone_trees_event_priority_genes = removeGenes(scicone_trees_event_priority_genes, genes_chr_y)
  scicone_trees_event_all_genes = removeGenes(scicone_trees_event_all_genes, genes_chr_y)

  highlighted_genes = {}
  for gene in getAffectedGenes(scicone_trees_event_priority_genes):
    highlighted_genes[gene] = "underline"
  for gene in dna_marker_genes:
    highlighted_genes[gene] = "bold"
  # Process trees. The results are in the last tree from the list.
  #processTrees([scicone_trees_event_priority_genes, scicone_trees_event_all_genes], dna_marker_genes, highlighted_genes, metadata, 
  #    "_".join([_OUT_DIR, "event"]), df_distances_csv=data["scicone_df_distances_event"])

  scicone_trees_clone_priority_genes = getCloneTrees(scicone_trees_event_priority_genes)
  scicone_trees_clone_all_genes = getCloneTrees(scicone_trees_event_all_genes)
  # Process trees. The results are in the last tree from the list.
  #processTrees([scicone_trees_clone_priority_genes, scicone_trees_clone_all_genes], dna_marker_genes, highlighted_genes, metadata, 
  #    "_".join([_OUT_DIR, "clone"]), df_distances_csv=data["scicone_df_distances_clone"])

  # SCATrEx.
  scatrex_trees_clone = readJsonTrees(data["scatrex_trees"], metadata, samples_to_remove)
  # Remove the root.
  for key, tree in scatrex_trees_clone.items():
    next_node = copy.deepcopy(tree.children[0])
    next_node.parent = None
    scatrex_trees_clone[key] = next_node

  #for node in PreOrderIter(scatrex_trees_clone["AFACYKY"]):
  #  if hasattr(node, "gene_cn_events") and "BCAT1" in node.gene_cn_events:
  #    print(node.node_id, node.gene_cn_events["BCAT1"]) 

  scicone_trees_clone_all_genes = matchSCATrExClones(scicone_trees_clone_all_genes, scatrex_trees_clone)

  '''
  for key, tree in scicone_trees_clone_all_genes.items():
    mean_ploidy = 0
    mean_exp = 0
    mean_xi = 0
    for node in PreOrderIter(tree):
      if hasattr(node, "gene_cn_events") and len(node.gene_cn_events.values()):
        mean_ploidy = mean_ploidy + statistics.mean(node.gene_cn_events.values())
        mean_exp = mean_exp + statistics.mean(node.absolute_exp_profile.values())
        mean_xi = mean_xi + statistics.mean(node.non_cn_exp_profile.values())
    num_clones = len(list(PreOrderIter(tree)))
    print(key, round(mean_ploidy/num_clones, 1), round(mean_exp/num_clones, 1), round(mean_xi/num_clones, 1))
  '''

  ########
  # Combine datasets for oncotreeVIS input. 
  '''
  combined_datasets = {}
  for key, tree in scatrex_trees_clone.items():
    combined_datasets[key + "_C_SCATrEx"] = tree
  for key, tree in scicone_trees_clone_all_genes.items():
    combined_datasets[key + "_B_SCICoNE_clone"] = tree
  for key, tree in scicone_trees_event_all_genes.items():
    combined_datasets[key + "_A_SCICoNE_event"] = tree
  combined_datasets = dict(sorted(combined_datasets.items()))
  oncotreeVIS_path = "_".join([_OUT_DIR, "scicone", "scatrex", "trees", "oncotreeVIS.json"])
  writeOncotreeVISInput(combined_datasets, metadata, [oncotreeVIS_path], dataset_name)  
  '''
  
  ###### TODO: assert that the same gene in not deleted and amplified in the same clone


  # Remove neutral clones.
  #removeNodes()


  printDatasetStatistics(scicone_trees_event_all_genes, "SCICoNE trees")
  printDatasetStatistics(scatrex_trees_clone, "SCATrEx trees")

  
