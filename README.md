# Network_Toxicology
This repository contains the source code for the article titled "TRAEC-guided disease module-based network toxicology and multi-omics validation reveal neurodevelopmental toxicity mechanisms of priority herbicides".

#1. Disease module network localization statistics:
The network localization statistics for disease modules use the code "Network_localization_statistics.py", which can provide "The Z-score of edges", "The P-value of edges", "The size of S", "The Z-score of S", and "The P-value of S".This code relies on the Python 3.9 environment.
Input files:
1) disease_genes_col.txt: Disease risk genes
2) Herbicide_target_col.txt: Herbicide targets
3) PPI.txt: Protein-protein interaction network

#2. Network Proximity Calculation
The network proximity analysis uses the code "proximity.py" to calculate the relative average shortest distance Zdc between the herbicide target network and the disease risk genes, to assess the neurodevelopmental toxicity risk of the herbicide. This code relies on the Python 2.7 environment.
Input files:
1) disease_genes_row.txt: Disease risk genes
2) Herbicide_target_row.txt: Herbicide targets
3) PPI.sif: Protein-protein interaction network
Output file:
proximity_results.csv: Network proximity calculation results

#3. Disease Module and Overlapping Network Construction
Using the code "disease_modules_and_overlapping_networks.py", calculate the connection relationships between nodes in the disease modules and overlapping networks, including the largest connected component, the nodes it contains, and the shortest distance ds between the nodes.This code relies on the Python 23.9 environment.
Input files:
1) Herbicide_target_col.txt: Herbicide targets. For constructing disease modules, input "disease_genes_col.txt". For overlapping network input, use "Herbicide_target_col.txt".
2) disease_genes_col.txt: Disease risk genes.
3) PPI.txt: Protein-protein interaction network.
Output files:
1) DM-PPI.txt: Connection relationships between nodes in the disease modules and overlapping networks.
2) DM-Smax1: Nodes included in the largest connected component.
3) DM_ds.tsv: Shortest distance ds between risk genes in the disease modules, or between disease genes and herbicide targets in the overlapping network.
