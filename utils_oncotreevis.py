import json
from anytree.exporter import JsonExporter
from anytree import PreOrderIter, RenderTree
import copy
import sys
from utils_anytree import isRoot

def createOncotreeVISInput(original_anytrees, metadata):
 
  print("Writing oncotreVIS input...") 

  anytrees = copy.deepcopy(original_anytrees)
  tree_map = {}
  exporter = JsonExporter(indent=2, sort_keys=False)
  for sample_name in anytrees:
    tree = anytrees[sample_name]
    for node in PreOrderIter(tree):
        del node.name
        if hasattr(node, 'region_cn_state'):
          del node.region_cn_state
        if hasattr(node, 'gene_state'):
          del node.gene_state
        if hasattr(node, 'node_label'):
          del node.node_label
        if hasattr(node, 'node_depth'):
          del node.node_depth

        gene_events = {}
        if hasattr(node, 'gene_cn_events'):
          for gene in node.gene_cn_events:
            gene_events[gene] = {}
            gene_events[gene]["CNA"] = node.gene_cn_events[gene]
          del node.gene_cn_events

        if hasattr(node, 'gene_exp_events'):
          for gene in node.gene_exp_events:
            if gene not in gene_events:
              gene_events[gene] = {}
            relative_expression_value = node.gene_exp_events[gene]
            if relative_expression_value > 1:
              gene_events[gene]["+expressed"] = relative_expression_value
            elif relative_expression_value < 1:
              gene_events[gene]["-expressed"] = relative_expression_value
          del node.gene_exp_events

        if not isRoot(node) or len(gene_events): 
          node.gene_events = gene_events
        
    json_tree = json.loads(exporter.export(tree))
    tree_map[sample_name] = {}
    tree_map[sample_name]["tree"] = json_tree 
    sample_name_prefix = sample_name.split("_")[0]
    if sample_name_prefix in metadata:
      tree_map[sample_name]["metadata"] = metadata[sample_name_prefix]
    else:
      tree_map[sample_name]["metadata"] = {}
  return tree_map

def writeOncotreeVISInput(anytrees, metadata, out_path, highlighted_genes={}):
  trees_oncotreevis = createOncotreeVISInput(anytrees, metadata)
  json_oncotreevis = {}
  json_oncotreevis["trees"] = trees_oncotreevis
  if len(highlighted_genes):
    json_oncotreevis["highlighted_genes"] = highlighted_genes

  file = open(out_path, "w")
  file.write(json.dumps(json_oncotreevis))
  file.close()



