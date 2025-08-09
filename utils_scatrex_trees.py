import json
import pandas as pd
import sys

from anytree.importer import JsonImporter
from anytree.exporter import JsonExporter
from anytree import AnyNode, RenderTree, PreOrderIter

def nodesMatch(node_scicone, node_scatrex):
  if node_scicone.depth != node_scatrex.depth or len(node_scicone.children) != len(node_scatrex.children):
    return False

  '''
  gene_set_scicone = set(node_scicone.gene_state.keys())
  gene_set_scatrex = set(node_scatrex.gene_cn_events.keys())
  if len(gene_set_scatrex.intersection(gene_set_scicone)) != len(gene_set_scatrex):
    #print(node_scicone.node_id, node_scatrex.node_id, "num genes", 
    #    len(gene_set_scatrex.intersection(gene_set_scicone)), len(gene_set_scatrex)) 
    return False

  if not gene_set_scatrex.intersection(gene_set_scicone) == gene_set_scatrex:
    #print(node_scicone.node_id, node_scatrex.node_id, "gene set not the same")
    return False

  gene_list_scicone = list(gene_set_scatrex.intersection(gene_set_scicone))
  gene_list_scicone.sort()
  gene_list_scatrex = list(gene_set_scatrex)
  gene_list_scatrex.sort()
  print("gene set scicone", len(gene_set_scatrex.intersection(gene_set_scicone)), gene_list_scicone) 
  print("gene set scatrex", len(gene_set_scatrex), gene_list_scatrex)

  for gene in gene_set_scatrex:
    if node_scatrex.gene_cn_events[gene] != node_scicone.gene_state[gene]:
      #print(node_scicone.node_id, node_scatrex.node_id, "CNAs don't match")
      return False
  ''' 
  return True

def getCorrespondingSCATrExNode(node_scicone, scatrex_tree):

  # Look for the scicone_node_id in SCATrEx first.
  for node_scatrex in list(PreOrderIter(scatrex_tree)):
    if node_scatrex.scicone_node_id and str(node_scicone.node_id) == str(node_scatrex.scicone_node_id):
      return node_scatrex

  # No corresponding labeling is providing, look for the most similar node.
  for node_scatrex in list(PreOrderIter(scatrex_tree)):
    if nodesMatch(node_scicone, node_scatrex):
      return node_scatrex
  return None

def getCorrespondingSCATrExSubclone(subclone_scicone, correspondence_map):
  split = subclone_scicone.split('_')
  sample_name = split[0]
  scatrex_node_id = correspondence_map[split[1]]
  return "_".join(sample_name, scatrex_node_id)

