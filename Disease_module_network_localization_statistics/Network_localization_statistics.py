import networkx as nx
import pandas as pd
import numpy as np
import seaborn as sns 
import sys
from scipy.special import ndtr
DEG_list=[]
dt=[]
ddt=[]
filename = open(r'disease.txt')
for dt1 in filename:
    dt1 = str(dt1).replace('\n', '')
    dt.append(dt1)
    ddt.append(dt1)
DEG_list=ddt
focal_genes = ddt
Gint=nx.Graph()
Gint = nx.read_edgelist(r"PPI.txt", create_using=nx.DiGraph())
#print(Gint.edges())
focal_genes = list(np.intersect1d(ddt, Gint.nodes()))
g = nx.Graph(nx.subgraph(Gint, ddt))
#---------------Loalization-----------------------------
num_reps=1000
def get_degree_binning(g, bin_size, lengths = None):
    degree_to_nodes = {}
    if sys.version_info >= (3, 0):
        for node, degree in dict(g.degree()).items():
            if lengths is not None and node not in lengths:
                continue
            degree_to_nodes.setdefault(degree, []).append(node)
    else:
        for node, degree in dict(g.degree()).iteritems():
            if lengths is not None and node not in lengths:
                continue
            degree_to_nodes.setdefault(degree, []).append(node)        
    values = list(degree_to_nodes.keys())
    values.sort()
    bins = []
    i = 0
    while i < len(values):        
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:           
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])            
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1 
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins
bins = get_degree_binning(Gint, 10)
min_degree, max_degree, genes_binned = zip(*bins)
bin_df = pd.DataFrame({'min_degree':min_degree, 'max_degree':max_degree, 'genes_binned':genes_binned})
actual_degree_to_bin_df_idx = {}
for i in range(0, bin_df['max_degree'].max() + 1):
    idx_temp = bin_df[ (bin_df['min_degree'].lt(i + 1)) & (bin_df['max_degree'].gt(i - 1)) ].index.tolist()
    if len(idx_temp) > 0:
        actual_degree_to_bin_df_idx[i] = idx_temp[0]
focal_genes = list(np.intersect1d(focal_genes, Gint.nodes())) 
numedges_list = []
numedges_rand = []
LCC_list = []
LCC_rand = []
sample_frac=1
method='both'
for r in range(num_reps):
    focal_80 = focal_genes
    np.random.shuffle(focal_80)
    focal_80 = focal_80[:int(len(focal_80)*sample_frac)]
    seed_random = []
    for g in focal_80:
        degree_temp = nx.degree(Gint,g)
        genes_temp = bin_df.loc[actual_degree_to_bin_df_idx[degree_temp]]['genes_binned']
        np.random.shuffle(genes_temp)
        seed_random.append(genes_temp[0])
    #-----------------------
    if (method == 'numedges') or (method == 'both'):
        # number edges calc on focal set
        numedges_temp = len(nx.subgraph(Gint,focal_80).edges())
        numedges_list.append(numedges_temp)
        numedges_temp_rand = len(nx.subgraph(Gint,seed_random).edges())
        numedges_rand.append(numedges_temp_rand)
    if (method == 'LCC') or (method == 'both'):
        G_sub_temp = nx.Graph(nx.subgraph(Gint, focal_80))
        G_sub_temp = max(nx.connected_components(G_sub_temp), key=len)
        LCC_list.append(len(G_sub_temp))     
        G_sub_temp = nx.Graph(nx.subgraph(Gint, seed_random))
        G_sub_temp = max(nx.connected_components(G_sub_temp), key=len)
        LCC_rand.append(len(G_sub_temp))
#LCC
if method == 'numedges':
    analysis_list = numedges_list
    analysis_rand = numedges_rand
else:
    analysis_list = LCC_list
    analysis_rand = LCC_rand
    S=np.mean(LCC_list)#the size of S
    analysis_z = (np.mean(analysis_list) - np.mean(analysis_rand))/float(np.std(analysis_rand))#LCC_Zs
    print('lcc_pvalue',1 - ndtr(analysis_z))    
#edge
    analysis_list = numedges_list
    analysis_rand = numedges_rand
    number_edges = np.mean(numedges_list)
    analysis_z = (np.mean(analysis_list) - np.mean(analysis_rand))/float(np.std(analysis_rand))#edge_Zs
    edge_pvalue=1 - ndtr(analysis_z)