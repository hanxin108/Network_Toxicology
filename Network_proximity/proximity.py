

import sys
import os
import pandas as pd
from toolbox import wrappers
if sys.getdefaultencoding() != 'utf-8':
  sys.path.append(r"toolbox")
file_name = r"PPI.sif"
filename1 = open(r'Herbicide_target.txt')
filename2 = open(r'disease.txt')
list=[]
for j in filename2:
        L=j.strip("\r\n").split('\t')
disease_targets=L[1:]
for i in filename1:
            L=i.strip("\r\n").split('\t')
           drug_targets=L[1:]
            network = wrappers.get_network(file_name, only_lcc = True)
        print(drug_targets)
        print(disease_targets)
        nodes_from =drug_targets
        nodes_to =disease_targets
        d, z, (mean, sd) = wrappers.calculate_proximity(network, nodes_from, nodes_to, min_bin_size = 2)
        print (d, z, (mean, sd))
        list.append([L[0],d, z,mean,sd])
        df=pd.DataFrame(list,columns=["Herbicide","d","Zdc","mean","sd"])
       df.to_csv("proximity_results.csv",index=False)
       
