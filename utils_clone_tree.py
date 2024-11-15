import copy
from anytree import PreOrderIter, RenderTree
from utils_anytree import isRoot

def getCloneTrees(anytrees):
  for key, tree in anytrees.items():
    for node in PreOrderIter(tree):
      # Remove auxiliary nodes (eg nodes with negative ids).
      '''
      if int(node.node_id) < 0:
        # Link the direct children to the parent.
        for child in node.children:
          child.parent = node.parent
        node.parent = None
      '''
      # Populate subclone gene state.
      if isRoot(node): # root
        node.gene_state = {}
      else:
        assert hasattr(node, 'gene_cn_events')
        # Copy the gene states from the parent node.
        gene_state = copy.deepcopy(node.parent.gene_state)
        # Update the gene states according to the current CN events.
        for gene, cn_event in node.gene_cn_events.items():
          new_gene_state = cn_event
          if gene in gene_state:
            new_gene_state = gene_state[gene] + cn_event
          gene_state[gene] = new_gene_state
          if new_gene_state == 0:
            del gene_state[gene]
        node.gene_state = gene_state

  clone_trees = copy.deepcopy(anytrees)
  for key, tree in clone_trees.items():
    for node in PreOrderIter(tree):
      if hasattr(node, 'gene_cn_events'):
        node.gene_cn_events = node.gene_state

  return clone_trees

def getChrTrees(anytrees, gene_chr_map):
  chr_trees = copy.deepcopy(anytrees)
  for sample_name, tree in chr_trees.items():
    for node in PreOrderIter(tree):
      if hasattr(node, "gene_cn_events"):
        chrs = set()
        for gene in node.gene_cn_events:
          chrs.add("chr" + gene_chr_map[gene])
        node.chrs = chrs
        del node.gene_cn_events
      if hasattr(node, 'gene_state'):
        del node.gene_state
      if hasattr(node, 'region_cn_state'):
        del node.region_cn_state
    print(RenderTree(tree))
    print("-------------") 
  return chr_trees
 
