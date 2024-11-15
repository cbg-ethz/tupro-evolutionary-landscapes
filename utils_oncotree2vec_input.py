import os
import copy
import json
import pandas as pd
import sys
import shutil

import networkx as nx
from anytree.importer import JsonImporter
from anytree import RenderTree, PreOrderIter, PostOrderIter
from utils import create_dir
from utils_anytree import isRoot

def generateTreeFiles(edges, features, sample_index, dir_path):

    create_dir(dir_path)

    input_map = {}
    input_map["edges"] = edges
    input_map["features"] = features

    json_path = os.path.join(dir_path, str(sample_index) + ".json")
    file = open(json_path, "w")
    file.write(str(input_map).replace('\'','\"'))
    file.close()

    gexf_path = os.path.join(dir_path, str(sample_index) + ".gexf")
    data = json.load(open(json_path))
    graph = nx.from_edgelist(data["edges"])

    features = data["features"]
    features = {int(k): v for k, v in features.items()}
    #features = {k: v for k, v in features.items()}
    node_attributes = {}
    for node in features:
      node_attributes[node] = {"Label" : features[node]}
    nx.set_node_attributes(graph, node_attributes)
    nx.write_gexf(graph, os.path.splitext(gexf_path)[0] + ".gexf")

def anytreeToOncotree2vecFile(root, dir_path, sample_index):
  edges = []
  features = {}
  for node in PreOrderIter(root):
    if not isRoot(node): # if not root
      edges.append([node.parent.node_id, node.node_id])
    features[str(node.node_id)] = node.matching_label
  generateTreeFiles(edges, features, sample_index, dir_path)

def saveOncotree2vecInput(anytrees, dir_path, duplicate_trees=False):

  print("Writing oncotre2vec input at", dir_path, "...")

  index_sample_map = {}
  paths = []
  sample_index = 0
  for sample_name, root in anytrees.items():
    anytreeToOncotree2vecFile(root, dir_path, sample_index)
    index_sample_map[sample_index] = sample_name
    sample_index = sample_index + 1

  # Add duplicates.
  if duplicate_trees:
    indices = list(index_sample_map.keys())
    for index in indices:
      source_path = os.path.join(dir_path, str(index) + ".gexf") 
      target_path = os.path.join(dir_path, str(sample_index) + ".gexf")
      shutil.copyfile(source_path, target_path)
    
      sample_name = index_sample_map[index]
      index_sample_map[sample_index] = sample_name + "#2"
      sample_index = sample_index + 1

  index_path = os.path.join(dir_path, "filename_index.csv")
  (pd.DataFrame.from_dict(data=index_sample_map, orient='index').to_csv(index_path, header=False))   
