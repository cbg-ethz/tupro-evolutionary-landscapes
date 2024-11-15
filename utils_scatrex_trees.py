import json
import pandas as pd
import sys

from anytree.importer import JsonImporter
from anytree.exporter import JsonExporter
from anytree import AnyNode, RenderTree, PreOrderIter

chr_y = ['RF00019', 'RF00156', 'RF00017', 'RF01518', 'RF00401', 'LL0YNC03-29C1.1', 'PLCXD1', 'GTPBP6', 'NCRNA00107', 'PPP2R3B', 'RP13-465B17.4', 'RP13-465B17.5', 'RP13-465B17.2', 'RP13-465B17.3', 'SHOX', 'RP11-309M23.1', 'RP11-309M23.2', 'CRLF2', 'CSF2RA', 'MIR3690', 'RNA5SP498', 'IL3RA', 'SLC25A6', 'NCRNA00106', 'NCRNA00105', 'ASMTL', 'P2RY8', 'AKAP17A', 'ASMT', 'RP13-297E16.4', 'RP13-297E16.5', 'DHRSX', 'RP13-858C7.1', 'ZBED1', 'MIR6089', 'CD99L1', 'NCRNA00102', 'CD99', 'SPRY3', 'AJ271735.3', 'AJ271736.1', 'VAMP7', 'AJ271736.8', 'AJ271736.3', 'IL9R', 'AJ271736.10', 'AJ271736.5', 'WASH6P', 'DDX11L16', 'RNU6-1334P', 'SRY', 'AC006040.4', 'AC006040.5', 'RPS4Y1', 'AC006040.7', 'AC006157.1', 'RP11-414C23.1', 'ZFY', 'AC006157.4', 'AC006157.5', 'AC006032.1', 'AC006152.1', 'AC019058.1', 'TGIF2LY', 'AC024038.1', 'AC012078.1', 'AC012078.2', 'AC012078.3', 'AC010084.1', 'RNU6-303P', 'AC019060.1', 'RP11-122L9.1', 'PCDH11Y', 'RNU2-57P', 'AC010685.1', 'AC010129.1', 'AC012667.1', 'SNX3P1Y', 'AC010874.1', 'TUSC2P1', 'AC010140.1', 'AC006335.6', 'TSPY2', 'FAM197Y9', 'AC006335.5', 'RP11-492C2.1', 'AC006335.9', 'AC006335.10', 'AC006335.11', 'AC006335.13', 'AC006335.14', 'TTTY1B', 'TTTY2B', 'TTTY21B', 'TTTY7', 'TTTY8B', 'AC010154.2', 'AC010144.1', 'AC013412.1', 'RP11-507A3.1', 'AC013412.2', 'AMELY', 'AC011297.1', 'TBL1Y', 'RP11-115H13.1', 'AC012068.2', 'PRKY', 'RN7SKP282', 'RNU6-941P', 'RNU6-521P', 'RP11-105L10.1', 'AC007274.1', 'AC007274.2', 'AC007274.4', 'AC007274.5', 'TTTY16', 'AC007275.1', 'TTTY12', 'RP11-109F19.1', 'AC007275.3', 'RP11-108F14.1', 'AC010678.1', 'AC010678.2', 'AC010902.1', 'AC016749.2', 'AC016749.1', 'AC051663.1', 'AC051663.3', 'AC025731.1', 'AC025731.2', 'AC025731.3', 'AC025731.4', 'AC025731.5', 'AC025731.6', 'AC025731.7', 'RP11-17E15.1', 'AC064829.1', 'TTTY18', 'TTTY19', 'TTTY11', 'AC007967.2', 'AC007967.3', 'RP11-373F14.1', 'AC007967.4', 'AC068719.1', 'AC009952.1', 'AC009952.2', 'AC009952.3', 'AC009952.4', 'RBMY1A3P', 'RP11-441G8.1', 'TTTY20', 'TSPY4', 'AC006158.5', 'TSPY8', 'FAM197Y7', 'AC006158.8', 'FAM197Y6', 'TSPY3', 'FAM197Y5', 'TSPY1', 'FAM197Y4', 'TSPY9P', 'FAM197Y3', 'AC006156.2', 'FAM197Y2', 'TSPY10', 'FAM197Y1', 'AC006156.5', 'AC025819.1', 'AC025819.2', 'AC025819.3', 'AC017019.1', 'TTTY8', 'TTTY7B', 'TTTY21', 'TTTY2', 'TTTY1', 'TTTY22', 'AC010891.2', 'AC010891.3', 'RP11-373F14.2', 'AC006986.1', 'AC006986.2', 'TTTY23', 'AC006986.4', 'AC006987.1', 'AC006987.2', 'AC006987.3', 'RP11-160K17.1', 'AC006987.4', 'RNA5SP518', 'AC006987.5', 'RNA5SP519', 'AC006987.6', 'AC006987.7', 'AC006987.8', 'AC010970.1', 'AC010970.2', 'RNA5-8SP6', 'RP1-85D24.4', 'RP1-85D24.5', 'RP1-85D24.3', 'RP1-85D24.2', 'RP1-85D24.1', 'RP11-131M6.1', 'MAFIP', 'RP11-131M6.3', 'RP11-131M6.4', 'DUX4L16', 'DUX4L17', 'DUX4L18', 'DUX4L19', 'PABPC1P5', 'RP11-886I11.6', 'RP11-295P22.1', 'RP11-295P22.3', 'RP11-295P22.5', 'RP11-295P22.4', 'RP11-295P22.2', 'AC011293.1', 'AC012502.1', 'AC011302.1', 'AC004772.1', 'AC004772.2', 'RN7SL702P', 'AC002992.1', 'AC002992.2', 'AC002992.3', 'GYG2P1', 'AC004617.1', 'AC004617.2', 'AC004810.1', 'TTTY15', 'USP9Y', 'AC002531.2', 'DDX3Y', 'AC004474.2', 'AC005820.1', 'AC005820.2', 'AC005820.3', 'UTY', 'AC010877.1', 'TMSB4Y', 'KALP', 'VCY', 'VCY1B', 'AC018677.2', 'AC010723.1', 'NLGN4Y', 'AC010979.2', 'AC010879.1', 'AC017032.1', 'AC006989.1', 'AC006989.2', 'AC006998.1', 'AC006382.1', 'RNU6-109P', 'RNU6-184P', 'AC006999.1', 'RP11-529I21.1', 'FAM41AY1', 'RP11-529I21.2', 'TUBB1P2', 'AC015978.3', 'RNA5SP520', 'RNA5SP521', 'AC068704.2', 'AC007742.1', 'RP11-357E16.1', 'AC007742.2', 'AC007742.3', 'AC007742.4', 'ACTG1P2', 'AC007742.7', 'AC007742.9', 'AC007742.10', 'RNU1-128P', 'AC007742.8', 'AC095381.1', 'AC095381.2', 'CDY2B', 'AC009976.2', 'AC009976.8', 'AC009976.3', 'AC009976.4', 'AC009976.5', 'AC009976.9', 'AC009976.7', 'CDY2A', 'AC095380.1', 'AC024183.1', 'AC024183.2', 'AC024183.9', 'RNU1-95P', 'AC024183.10', 'AC024183.3', 'AC024183.5', 'AC024183.6', 'AC024183.7', 'AC024183.8', 'TAF9P2', 'CLUHP2', 'AC007241.4', 'RNA5SP522', 'RNA5SP523', 'TUBB1P1', 'RP11-157F24.1', 'FAM41AY2', 'TCEB1P13', 'OFD1P4Y', 'AC068541.1', 'RNU1-48P', 'AC068541.2', 'AC022486.1', 'AC068541.3', 'AC022486.2', 'AC022486.3', 'AC022486.4', 'HSFY1', 'AC022486.8', 'TTTY9B', 'OFD1P5Y', 'AC022486.9', 'AC022486.10', 'OFD1P6Y', 'HSFY2', 'AC007379.4', 'AC007379.6', 'AC007379.7', 'AC007379.8', 'AC007379.9', 'AC007379.10', 'AC007379.11', 'RNU1-41P', 'OFDYP6', 'AC009235.2', 'NCRNA00185', 'AC009235.4', 'CD24L4', 'RNU6-255P', 'RP11-118E9.1', 'AC010133.1', 'BCORP1', 'TXLNG2P', 'RP11-576C2.1', 'RP11-424G14.1', 'KDM5D', 'RP11-356K22.1', 'AC009233.1', 'AC009233.2', 'AC009240.1', 'TTTY10', 'RP11-256K9.1', 'EIF1AY', 'AC009494.1', 'RPS4Y2', 'AC009494.4', 'AC009494.3', 'AC009494.5', 'AC009489.1', 'RP11-65G9.1', 'AC007876.1', 'AC009239.1', 'AC010086.1', 'PRORY', 'RBMY2EP', 'AC010086.5', 'AC010086.6', 'AC010141.1', 'AC010141.2', 'RBMY1H', 'AC010141.4', 'RBMY1B', 'RBMY1A1', 'AC010141.6', 'AC010141.8', 'TTTY13', 'RP11-178M5.1', 'AC021107.2', 'AC021107.3', 'AC021107.4', 'AC021107.5', 'OFD1P16Y', 'AC007322.1', 'RBMY1D', 'AC007322.3', 'RBMY1E', 'AC007322.5', 'RBMY2AP', 'AC007322.7', 'AC007322.8', 'AC007359.1', 'AC007359.2', 'CDY12P', 'TCEB1P15', 'PRY2', 'AC007359.6', 'TTTY6B', 'RBMY1F', 'TSPY23P', 'RBMY2UP', 'AC025227.2', 'TTTY5', 'TSPY22P', 'RBMY2FP', 'TTTY25P', 'TSPY21P', 'RBMY1J', 'AC008175.1', 'TTTY6', 'PRY', 'AC008175.5', 'AC008175.6', 'AC008175.7', 'AC008175.8', 'AC008175.9', 'RP11-427G18.1', 'RBMY2BP', 'AC016694.2', 'AC010080.1', 'TTTY17A', 'RP11-473E1.1', 'AC006366.1', 'AC006366.5', 'BPY2', 'AC006366.3', 'AC006366.4', 'AC010088.1', 'AC010088.2', 'DAZ1', 'DAZ2', 'AC006983.2', 'AC009947.1', 'AC009947.2', 'AC009947.3', 'AC009947.5', 'AC016752.1', 'AC016752.2', 'AC016752.3', 'AC016752.4', 'AC016752.5', 'PRYP3', 'AC016752.7', 'AC016752.8', 'AC016752.9', 'RNU1-97P', 'AC016752.10', 'AC016752.11', 'AC025246.1', 'AC073649.1', 'AC073649.2', 'AC073649.3', 'AC073893.1', 'AC073893.2', 'TTTY3B', 'AC073893.3', 'RNU1-86P', 'AC073893.4', 'AC073893.5', 'AC068601.1', 'AC068601.2', 'CDY1B', 'AC068601.4', 'AC023274.1', 'AC023274.2', 'AC023274.3', 'AC023274.4', 'AC023274.5', 'DNM1P24', 'AC023274.6', 'CSPG4LYP1', 'AC023274.7', 'RP11-533E23.1', 'AC012005.2', 'RN7SL818P', 'AC012005.3', 'AC012005.4', 'AC012005.6', 'AC012005.5', 'RP11-533E23.2', 'AC016698.1', 'AC016698.2', 'AC016698.3', 'AC016698.4', 'AC010153.2', 'BPY2B', 'AC010153.3', 'AC025735.3', 'AC025735.2', 'AC006338.2', 'DAZ3', 'DAZ4', 'RP11-539D10.1', 'AC006338.3', 'RP11-539D10.2', 'AC016728.1', 'BPY2C', 'TTTY4C', 'AC006386.1', 'RP11-566H16.1', 'AC006386.2', 'AC006386.3', 'RP11-566H16.2', 'AC006328.1', 'AC006328.3', 'AC006328.2', 'AC006328.4', 'AC006328.5', 'RN7SL725P', 'RP11-102O5.1', 'AC006328.7', 'AC006328.6', 'AC006328.8', 'AC006328.9', 'AC006328.10', 'AC007562.1', 'AC007562.2', 'AC007562.3', 'AC007562.4', 'AC007562.5', 'CDY1', 'RP11-251M8.1', 'AC010682.1', 'AC010682.2', 'AC010682.3', 'AC010682.4', 'RNU1-107P', 'TTTY3', 'AC010682.6', 'AC010682.5', 'AC010682.7', 'AC017005.1', 'AC017005.2', 'AC017005.3', 'AC007965.1', 'AC007965.2', 'AC007965.3', 'RNU1-40P', 'AC007965.4', 'AC007965.5', 'PRYP4', 'AC007965.7', 'AC007965.8', 'AC007965.9', 'AC007965.10', 'AC007965.11', 'AC006991.1', 'AC006991.2', 'AC024067.1', 'AC024067.2', 'RNU6-1314P', 'AC013734.1', 'AC013734.3', 'RP11-557B9.1', 'AC019099.1', 'AC019099.2', 'AC019099.4', 'AC019099.5', 'AC019099.6', 'AC025226.2']

def nodesMatch(sample_name, node_scicone, node_scatrex):
  if node_scicone.depth != node_scatrex.depth or len(node_scicone.children) != len(node_scatrex.children):
    return False

  gene_set_scicone = set(node_scicone.gene_state.keys())
  gene_set_scatrex = set(node_scatrex.gene_cn_events.keys())
  if len(gene_set_scatrex.intersection(gene_set_scicone)) != len(gene_set_scatrex):
    print(node_scicone.node_id, node_scatrex.node_id, "num genes", len(gene_set_scatrex.intersection(gene_set_scicone)), len(gene_set_scatrex)) #, gene_set_scatrex - gene_set_scicone)
    if sample_name == "OQEFIRU":
      print("scatrex in chr y:", gene_set_scatrex.intersection(set(chr_y)))
      print("scicone genes:", gene_set_scicone)
      print("scatrex genes:", gene_set_scatrex)
    return False

  if not gene_set_scatrex.intersection(gene_set_scicone) == gene_set_scatrex:
    print(node_scicone.node_id, node_scatrex.node_id, "gene set not the same")
    return False

  '''
  gene_list_scicone = list(gene_set_scatrex.intersection(gene_set_scicone))
  gene_list_scicone.sort()
  gene_list_scatrex = list(gene_set_scatrex)
  gene_list_scatrex.sort()
  print("gene set scicone", len(gene_set_scatrex.intersection(gene_set_scicone)), gene_list_scicone) 
  print("gene set scatrex", len(gene_set_scatrex), gene_list_scatrex)
  ''' 

  for gene in gene_set_scatrex:
    if node_scatrex.gene_cn_events[gene] != node_scicone.gene_state[gene]:
      print(node_scicone.node_id, node_scatrex.node_id, "CNAs don't match")
      return False
  return True

# Match SCICoNE and SCATrEx subclones.
def matchSCATrExClones(clone_trees_scicone, clone_trees_scatrex):
 
  for key, scicone_tree in clone_trees_scicone.items():
    print()
    if key in clone_trees_scatrex:
      scicone_nodes = [node for node in list(PreOrderIter(scicone_tree)) if node.parent and hasattr(node, 'gene_state')]
      scatrex_nodes = [node for node in list(PreOrderIter(clone_trees_scatrex[key])) if node.parent and hasattr(node, 'gene_cn_events')]
      not_matched_scicone_nodes = []
      for node_scicone in scicone_nodes:
        found_subclone_match = False
        for node_scatrex in scatrex_nodes:
          #node_key_1 = key + "_" + str(node_scicone.node_id)
          #node_key_2 = key + "_" + str(node_scatrex.node_id)
          if node_scatrex.scicone_node_id != -1 and str(node_scicone.node_id) == str(node_scatrex.scicone_node_id):
            if nodesMatch(key, node_scicone, node_scatrex):
              print("found match", key, node_scicone.node_id, node_scatrex.node_id)
              node_scicone.mean_cn = node_scatrex.mean_cn
              node_scicone.absolute_exp_profile = node_scatrex.absolute_exp_profile
              node_scicone.non_cn_exp_profile = node_scatrex.non_cn_exp_profile

              found_subclone_match = True
              scatrex_nodes.remove(node_scatrex)
              break    

        if not found_subclone_match:
          not_matched_scicone_nodes.append(node_scicone)
          print("no SCATrEx corresponding subclone found for clone", key, "node_id", node_scicone.node_id)
      '''
      assert len(not_matched_scicone_nodes) == len(scatrex_nodes) - 1
      if len(scatrex_nodes) == 2: # one not matched node left (+ the root)
        node_scicone = not_matched_scicone_nodes[-1]
        node_scatrex = scatrex_nodes[-1]
        node_scicone.mean_cn = node_scatrex.mean_cn
        node_scicone.absolute_exp_profile = node_scatrex.absolute_exp_profile
        node_scicone.non_cn_exp_profile = node_scatrex.non_cn_exp_profile

        node_id_1 = node_scicone.node_id
        node_id_2 = node_scatrex.node_id
        #print("proposed node match:", key, node_id_1, node_id_2)
        
        #genes_node_1 = set(node_scicone.gene_state.keys())
        #genes_node_2 = set(node_scatrex.gene_cn_events.keys())
        #print("!!! non-matchign node ids:", node_id_1, node_id_2, "; #genes in scatrex and not in scicone", len(genes_node_2 - genes_node_1)) 
      '''  
  return clone_trees_scicone
