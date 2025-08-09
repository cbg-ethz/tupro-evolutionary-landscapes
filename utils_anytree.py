import copy
from anytree import PreOrderIter, RenderTree


def isRoot(node):
  return node.parent is None


def removeTreeAttributes(root, excepted_fields=[]):
  anytree = copy.deepcopy(root)
  for node in PreOrderIter(anytree):
    # Remove all attributes except node_id.
    for key in list(vars(node).keys()):
      if not key.startswith("_") and key != "node_id" and key not in excepted_fields:
        delattr(node, key)
  return anytree

