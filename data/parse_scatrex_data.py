import os 
import pandas as pd
import numpy as np
import json
from anytree import AnyNode, RenderTree, PreOrderIter
from anytree.exporter import JsonExporter
import statistics
import matplotlib.pyplot as plt
import seaborn as sns
import time
import sys

def plotDistribution(values, filename_prefix, filename_suffix):
  # Plot distribution of gene differential expressions.
  plt.figure()
  sns.set(font_scale=1)
  histogram = sns.histplot(data = values)
  plt.yscale('log')
  histogram.set(ylabel='Log counts', xlabel='Expression levels')
  histogram_fig = histogram.get_figure()
  histogram_filename = filename_prefix + "_expression-values_" + filename_suffix + ".png"
  histogram_fig.savefig(histogram_filename, format='png', dpi=300)

def parse_scatrex_trees(dir_path, metadata_file, root_label="Root"):
  has_wgd = {}
  if os.path.isfile(metadata_file):
    metadata = pd.read_csv(metadata_file).to_dict("index")
    for key, row in metadata.items():
      sample_name = row["sample_name"]
      has_wgd[sample_name] = True if row["has_wgd"]=="yes" else False

  anytrees_json = {}
  root_label = "Root"
  xi_values = []  
  total_exp_values = []
  neutral_gene_expression_interval = [1,1]
  exporter = JsonExporter(indent=2, sort_keys=False)
  for filename in os.listdir(dir_path):
    sample_name = filename.split('.')[0].split("-")[0]
    cn_neutral_state = 2
    if has_wgd[sample_name]:
      cn_neutral_state = 4

    with open(os.path.join(dir_path, filename)) as json_file:
      nodes = json.load(json_file)

    genes = nodes["genes"]
    node_id_object_map = {}
    # Add empty root.
    root = AnyNode(name=root_label, node_id=root_label)
    node_id_object_map[root_label] = root

    for key in nodes:
      if key == "genes": 
        continue
      node_id = key
      node = nodes[node_id]
      if node["parent"] == "NULL":
        parent_id = root_label
      else:
        parent_id = node["parent"]

      # CNA.
      absolute_gene_cn_events = dict(zip(genes, node["cnv"]))
      relative_gene_cn_events = {} 
      for key, value in absolute_gene_cn_events.items():
        relative_cn_value = value - cn_neutral_state
        if relative_cn_value:
          relative_gene_cn_events[key] = relative_cn_value

      # Gene expression.
      xi = dict(zip(genes, node["xi"]))
      xi_values.extend(list(xi.values()))

      if np.isnan(list(xi.values())).any():
        print(sample_name, node_id)

      total_exp = dict(zip(genes, node["total_exp"]))
      total_exp_values.extend(list(total_exp.values()))
      relative_expression = {}
      for key, exp_value in total_exp.items():
        if exp_value > neutral_gene_expression_interval[1]:
          relative_expression[key] = round(exp_value,2)
        elif exp_value < neutral_gene_expression_interval[0]:
          relative_expression[key] = -abs(round(exp_value,2))

      #gene_exp = [gene + "_up" for gene in xi if xi[gene] > neutral_gene_expression_interval[1]]
      #gene_exp.extend([gene + "_down" for gene in xi if xi[gene] < neutral_gene_expression_interval[0]])
      #gene_exp.sort()

      scicone_node_id = -1
      if "node_id" in node:
        scicone_node_id = node["node_id"]

      # Populate anytree node.
      anytree_node = AnyNode(name=node_id,
          node_id = node_id,
          scicone_node_id = scicone_node_id,
          parent_id = parent_id,
          mean_cn = statistics.mean(absolute_gene_cn_events.values()),
          gene_cn_events = relative_gene_cn_events,
          gene_exp_events = relative_expression, 
          absolute_exp_profile = total_exp,
          non_cn_exp_profile = xi)
      node_id_object_map[node_id] = anytree_node

    # Populate parent and node_depth.
    for node_id, node in node_id_object_map.items():
      if node_id == root_label: # if root
        continue
      node.parent = node_id_object_map[node.parent_id]
      del node.parent_id

    # Replace the node_ids with integers.
    idx = 2
    for node in PreOrderIter(root):
      if node.node_id == root_label: # if root
        node.node_id = 0
      else:
        node.node_id = idx
        idx = idx + 1
    anytrees_json[sample_name] = json.loads(exporter.export(root))

  timestamp = str(int(time.time()))
  filename_prefix = "out/" + timestamp + "_" + os.path.basename(dir_path)
  tree_filename = filename_prefix + "_tota_exp_" + ':'.join(str(x) for x in neutral_gene_expression_interval) + ".json"
  file = open(tree_filename, "w")
  file.write(json.dumps(anytrees_json))

  # Plot distribution of gene differential expressions.
  plotDistribution(xi_values, filename_prefix, "xi")
  plotDistribution(total_exp_values, filename_prefix, "total_exp")

if __name__ == "__main__":
  #parse_scatrex_trees("scatrex_melanoma_18june2024", "metadata_melanoma.csv")
  parse_scatrex_trees("scatrex_ovarian_14nov2024b", "metadata_ovarian.csv")
