import json
import os
import copy
import pandas as pd
from anytree import PreOrderIter
from anytree.importer import JsonImporter
from utils_anytree import isRoot
import sys

def isNaN(num):
  return num != num

# Read gene-chromosome correspondance from file.
def getGeneToChrMap(csv, key):
  if os.path.isfile(csv):
    gene_chr_map = {}
    dict = pd.read_csv(csv).to_dict(key)
    for key, pair in dict.items(): 
      gene_chr_map[pair["gene"]] = pair["chr"]
    return gene_chr_map

# Read metadata dict from file.
def readMetadata(metadata_file):
  metadata = {}
  if os.path.isfile(metadata_file):
    metadata = pd.read_csv(metadata_file).to_dict("index")
    keys = copy.deepcopy(list(metadata.keys()))
    for key in keys:
      sample_name_key = ""
      if "sample_name" in metadata[key]:
        sample_name_key = "sample_name"
      elif "Patient_ID" in metadata[key]:
        sample_name_key = "Patient_ID"
      else:
        break
      sample_name = metadata[key][sample_name_key]
      metadata[sample_name] = metadata[key]
      for kkey in metadata[key]:
        if isNaN(metadata[key][kkey]):
          metadata[sample_name][kkey] = ""
      del metadata[sample_name][sample_name_key]
      del metadata[key]
  return metadata

# Load anytree data structure from json file.
def readJsonTrees(path, metadata, samples_to_remove=[], neutral_clones=[]):
  importer = JsonImporter()
  with open(path) as json_file:
    sample_map = json.load(json_file)
    for key in samples_to_remove:
      sample_map.pop(key, None)
    anytrees = {}
    for key, json_tree in sample_map.items():
      if "event_tree" in json_tree:
        json_tree = json_tree["event_tree"]
      anytrees[key] = importer.import_(json.dumps(json_tree))
      for node in PreOrderIter(anytrees[key]):
        # CN neutral state.
        node.cn_neutral_state = 4 if metadata[key]["has_wgd"]=="yes" else 2
        # Is neutral subclone.
        if isRoot(node) or node.node_id == -100:
          node.is_neutral = True
          continue
        else: 
          if hasattr(node,"node_label") and '_'.join([key, str(node.node_label)]) in neutral_clones:
            node.is_neutral = True
            continue
        node.is_neutral = False
    return anytrees


# Remove the genes from the gene list from all the nodes in the input trees. 
def removeGenes(anytrees, gene_list):
  for sample, anytree in anytrees.items():
    for node in PreOrderIter(anytree):
      if hasattr(node, 'gene_cn_events'):
        genes_to_remove = set(node.gene_cn_events.keys()).intersection(gene_list)
        for gene in genes_to_remove:
          del node.gene_cn_events[gene]
  return anytrees


# Get statistics.
def printDatasetStatistics(anytrees, dataset_name):
  num_clones = 0
  sum_degree = 0
  for key, root in anytrees.items():
    num_clones = num_clones + len(list(PreOrderIter(root)))
    for node in PreOrderIter(root):
      sum_degree = sum_degree + len(node.children)
      if not isRoot(node):
        sum_degree = sum_degree + 1
  print()
  print("Overview", dataset_name)
  print("Num trees:", len(anytrees.keys()))
  print("Num clones:", num_clones)
  print("Avg node degree:", sum_degree / num_clones)

def printMetadataStatistics(metadata, samples):

  metadata_subset = {k: metadata[k] for k in samples if k in metadata}
  df = pd.DataFrame.from_dict(metadata_subset, orient="index")
  for col in df.columns:
    print(f"\nValue counts for '{col}':")
    print(df[col].value_counts())  
  print()




