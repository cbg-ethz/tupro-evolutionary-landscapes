import numpy as np
from anytree import Node, PreOrderIter
from utils_anytree import isRoot

def getLabelSet(node_map, labels, value):
  if value < len(labels):
    return [labels[int(value)]]
  else:
    return getLeaves(node_map[value], labels)

def getLeaves(root, labels):
  num_labels = len(labels)
  output = []
  for node in PreOrderIter(root):
    if node.name < num_labels:
      output.append(labels[node.name])
  return output

# Get max over all elements except the diagonal.
def getMaxPairwiseDistance(df_distances, label_set):
  submatrix = df_distances.loc[label_set, label_set]
  submatrix = submatrix.mask(np.eye(len(submatrix.index), dtype = bool))
  return submatrix.max().max()

def mergeChildren(node, threshold):
  if hasattr(node, "weight"):
    if isRoot(node) and node.weight <= threshold: # if root
      return True # merge everything
    if node.weight <= threshold and node.parent.weight > threshold:
      return True
    else:
      return False
  else: # leaf
    if node.parent.weight > threshold:
      return True
  return False

## Get the clusters from hierarchical clustering based on pairwise similarity.
def getSimilarityClusters(hierarchy_linkage, df_distances, distance_threshold=1):
  labels = list(df_distances.index)
  num_labels = len(labels)
  #assert num_labels - 1 == len(hierarchy_linkage)

  node_map = {}
  root = None

  for i, pair in enumerate(hierarchy_linkage):
    # Create tree with distance edge weights.
    distance = -1
    if pair[0] < num_labels and pair[1] < num_labels:
      distance = df_distances[labels[int(pair[0])]].loc[labels[int(pair[1])]]
    else:
      label_set_0 = getLabelSet(node_map, labels, pair[0])
      label_set_1 = getLabelSet(node_map, labels, pair[1])
      distance = getMaxPairwiseDistance(df_distances, label_set_0 + label_set_1)

    parent = Node(i + num_labels, weight = distance)
    node_map[i + num_labels] = parent

    if pair[0] < num_labels:
      node_0 = Node(int(pair[0]), parent=parent, sample_name=labels[int(pair[0])])
    else:
      node_map[pair[0]].parent = parent

    if pair[1] < num_labels:
      node_1 = Node(int(pair[1]), parent=parent, sample_name=labels[int(pair[1])])
    else:
      node_map[pair[1]].parent = parent

    if i == num_labels - 2:
      root = parent

  clusters = []
  for node in PreOrderIter(root):
    # if root and weight below threshold, merge everything.
    if isRoot(node) and node.weight <= distance_threshold:
      return [getLeaves(node, labels)]
    if not isRoot(node): # if not root
      if mergeChildren(node, distance_threshold):
        cluster = getLeaves(node, labels)
        if len(cluster) != 1:
          assert getMaxPairwiseDistance(df_distances, cluster) <= distance_threshold
        clusters.append(cluster)

  return clusters
