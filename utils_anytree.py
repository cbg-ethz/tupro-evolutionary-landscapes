def isRoot(node):
  return node.parent is None

'''
def not_root(node_id):
  return node_id != 0 and node_id != -100

def get_node_by_id(root, id):
  for node in PreOrderIter(root):
    if node.node_id == int(id):
      return node
  return None

def get_tree_nodes(root):
  node_ids = []
  for node in PreOrderIter(root):
    node_ids.append(node.node_id)
  return node_ids

def get_node_depth(clone_id, root):
  node_id = get_id_from_clone_name(clone_id)
  node = get_node_by_id(root, node_id)
  return node.node_depth

def get_node_degree(node):
  return len(node.children) + 1

def get_path_to_root(root, node):
  nodes = []
  while node != root:
    nodes.append(node)
    node = node.parent
  return nodes

def filter_node_attributes(node):
  del node.name
  if hasattr(node, 'node_label'):
    del node.node_label
  if hasattr(node, 'num_cells'):
    del node.num_cells
  #if hasattr(node, 'gene_cn_events'):
  #  del node.gene_cn_events
  #if hasattr(node, 'gene_state'):
  #  del node.gene_state

def filter_tupro_node_attributes(root):
  for node in PreOrderIter(root):
    filter_node_attributes(node)

def get_clone_names_from_tree(sample, root, event_tree=False):
  ids = []
  for node in PreOrderIter(root):
    # TODO: write nicer code.
    if event_tree:
      if node.node_id and node.is_malignant:
        ids.append(sample + "_" + str(node.node_id))
    else:
      if node.node_id and node.node_label and node.is_malignant: #len(node.gene_state):
        ids.append(sample + "_" + str(node.node_id))
  return ids

def make_root(node):
  node.parent = None
  node.is_root
  return node

def tree_height(root):
  return max([node.node_depth for node in root.leaves])

def setNodeColor(node_list, node_id, color):
  for node in node_list:
    if node["id"] == node_id:
      node["color"] = color
      break

def AddXOffsets(tree):
  num_leaves = len(tree.leaves)
  max_depth = max([node.node_depth for node in tree.leaves])
  
  offset = 10
  iter = 0
  for node in tree.leaves:
    node.dx = iter * offset
    iter = iter + 1

  for node in PostOrderIter(tree):
    if node not in tree.leaves:
      sum_children_dx = 0
    for child in node.children:
      sum_children_dx = sum_children_dx + child.dx
      node.dx = sum_children_dx / len(node.children)
  return tree

#####################################
###### Helpers tree matching  #######
#####################################

class CloneMatch:
  # init method or constructor
  def __init__(self, clone_1, clone_2, similarity):
    self.clone_1 = clone_1
    self.clone_2 = clone_2
    self.similarity = similarity
 
  def print(self):
    print(self.clone_1, self.clone_2, self.similarity)

def EdgeToMap(source, target, similarity=-1, color="#999"):
  return {'source': source, 'target': target, 'similarity': similarity, 'color': color}

def NodeToMap(sample, node, label, color="white"):
  node = {'id': sample + "_" + str(node.node_id),
      'label': label, 'depth': node.node_depth, 'dx': node.dx, 'matching_label': node.matching_label,  'color': color}
  if hasattr(node, 'size_percent'):
    node["size_percent"] = node.size_percent
  return node

def AnytreeToMap(sample, root):
  nodes = []
  links = []
  root = AddXOffsets(root)
  for node in anytree.PreOrderIter(root):
    if node.node_id and node.parent: #check if not root
      links.append(EdgeToMap(sample + "_" + str(node.parent.node_id), sample + "_" + str(node.node_id)))

    if hasattr(node, 'node_label'):
      label = node.node_label
    else:
      label = ""

    nodes.append(NodeToMap(sample, node, label))
  return {'nodes': nodes, 'links': links}

def get_root_to_leaf_paths(root):
  # If leaf, return (it's the end of the path).
  if not root.children: 
    return [[root]]

  paths = []
  for node in root.children:
    descending_paths = get_root_to_leaf_paths(node)
    for idx in range(len(descending_paths)):
      descending_paths[idx] = [root] + descending_paths[idx]
    paths.append(descending_paths) 
  paths = [item for sublist in paths for item in sublist]
  return paths

# Assign clone labels. Anytrees should contain the attribute "labels" with the list of labels,
# and the root need to have the node name `root_label`.
# `clone_gene_map` and `distance_threshold` are only used when `use_multi_labels`.
def assign_clone_labels(anytrees, cancer_type, tree_type="", use_multi_labels=False, clone_gene_map=None, 
    distance_threshold=None, plot_heatmaps=None, output_dir=None, marker_genes_cn=None, marker_genes_exp=None):
  if not use_multi_labels: # each node has a single label
    gene_label_map = {}
    for sample_name, tree in anytrees.items():
      for node in PreOrderIter(tree):
        if not node.parent: # if root
          node.matching_label = 0
        else:
          if type(node.label) is str:
            node.label = [node.label]
          node_label = node.label[0]
          if node_label not in gene_label_map:
            gene_label_map[node_label] = len(gene_label_map) + 2 # Labels 0 and 1 are reserved for root and neutral nodes respectively.
          node.matching_label = gene_label_map[node_label]
    if plot_heatmaps:
      plot_label_histogram(anytrees, cancer_type, output_dir)

  else: # Relabel the nodes of the clone tree according to the Jaccard distance threshold.
    if marker_genes_cn:
      marker_genes_cn = [[gene+"_amp", gene+"_del"] for gene in marker_genes_cn]
      marker_genes_cn = [item for sublist in marker_genes_cn for item in sublist]      
    else: 
      marker_genes_cn = []

    if marker_genes_exp:
      marker_genes_exp = [[gene+"_up", gene+"_down"] for gene in marker_genes_exp]
      marker_genes_exp = [item for sublist in marker_genes_exp for item in sublist]
    else:
      marker_genes_exp = []

    clones = list(clone_gene_map.keys())

    # Save the distance matrix in an auxiliary file.
    aux_file = "df_jaccard_distances_" + cancer_type + "_" + tree_type + ".csv"
    if os.path.exists(aux_file):
      print("Loading distance matrix from ", aux_file)
      df_distances = pd.read_csv(aux_file, index_col=0)
    else:
      print(aux_file, "does not exist. Computing the distance matrix from the scratch.")

      def get_jaccard_similarity(set_1, set_2):
        set_1 = set(set_1)
        set_2 = set(set_2)
        intersection = len(set_1.intersection(set_2))
        union = len(set_1.union(set_2))
        assert union > 0
        return intersection / union

      df_distances = pd.DataFrame(0, columns=clones, index=clones).astype(float)
      for clone_1 in clones:
        list_1_genes = clone_gene_map[clone_1]
        if marker_genes_cn or marker_genes_exp:
          list_1_marker_genes = set(list_1_genes).intersection(set(marker_genes_cn + marker_genes_exp))

        for clone_2 in clones:
          if clone_1 == clone_2:
            df_distances[clone_1][clone_2] = 1
            continue
          list_2_genes = clone_gene_map[clone_2]
          jaccard_sim = get_jaccard_similarity(list_1_genes, list_2_genes)

          jaccard_sim_marker = 0
          if marker_genes_cn or marker_genes_exp:
            list_2_marker_genes = set(list_2_genes).intersection(set(marker_genes_cn + marker_genes_exp))
            if list_1_marker_genes.union(list_2_marker_genes):
              jaccard_sim_marker = get_jaccard_similarity(list_1_marker_genes, list_2_marker_genes)
           
          df_distances[clone_1][clone_2] = 1 - max(jaccard_sim, jaccard_sim_marker)
      df_distances.to_csv(aux_file)

    clustering = hierarchy.linkage(
        distance.pdist(df_distances), metric="euclidean", method="ward")
    clone_clusters = get_similarity_clusters(clustering, clones, df_distances, distance_threshold)

    clone_labels = {}
    for clone_cluster in clone_clusters:
      label = len(clone_labels) + 2 # Labels 0 and 1 are reserved for root and neutral nodes respectively.
      for clone in clone_cluster:
        clone_labels[clone] = label

    for sample_name, root in anytrees.items():
      for node in PreOrderIter(root):
        if not node.parent: # if root
          node.matching_label = 0
        else:
          clone_name = sample_name + "_" + str(node.name)
          if clone_name in clone_labels:
            node.matching_label = clone_labels[clone_name]
          else:
            node.matching_label = 1 # neutral node

    if plot_heatmaps:
      generate_heatmaps(df_distances, clustering, clone_clusters, 0.5, output_dir,
          cancer_type, "event", plot_cluster_heatmaps=False)
      #generate_heatmaps(clustering, clone_clusters, clones, df_distances, output_dir,
      #    cancer_type, tree_type, distance_threshold, plot_cluster_heatmaps=False)

# Propagate the labels down to the children. The tree should contain
# the attribute "label" with the list of labels.
def create_clone_trees(anytrees):
  clone_gene_map = {}
  anytrees_clone = {}
  for sample_name, tree in anytrees.items():
    tree_copy = copy.deepcopy(tree)
    for node in PreOrderIter(tree_copy):
      if not node.parent: # if root
        node.label = []
      else:
        node.label = list((set(node.parent.label)).union(set(node.label)))
        if node.label:
          clone_gene_map[sample_name + "_" + str(node.node_id)] = node.label
    anytrees_clone[sample_name] = tree_copy

  return anytrees_clone, clone_gene_map

def removeNodes(anytrees, function):
  for key, tree in anytrees.items():
    print()
    print(before)
    print(RenderTree(tree))
    for node in PreOrderIter(tree):
      if func(node):
        print("******************************")
        parent = node.parent
        for child in node.children:
          child.parent = parent
        node.parent=None
    print(after)
    print(RenderTree(tree))
    print("=========")
''' 
