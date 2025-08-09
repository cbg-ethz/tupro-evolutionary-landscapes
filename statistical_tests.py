from scipy.stats import fisher_exact, chi2_contingency
#from lifelines import KaplanMeierFitter
#from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt

#######
# Testing that: "“WGD is less common in early-divergent tumors.”
#                WGD
#               Yes   No  
# EarlyDiv  Yes  a=0     b=15    
#           No   c=20    d=12   # I exclude the 2 uncertain ones

contingency_table = [[0, 15],
                     [20, 12]]
print(fisher_exact(contingency_table, alternative='less'))

# Survival analisys for your two groups — tumors with vs. without early divergence.

