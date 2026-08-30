import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

dt=[]
ddt=[]
filename = open(r'C:\Users\11445\Desktop\Atrazine_target_col.txt')
for dt1 in filename:
    dt1 = str(dt1).replace('\n', '')
    dt.append(dt1)
    ddt.append(dt1)
#The number of herbicide targets (overlapping network) or the number of disease genes (constructing disease modules)
print("the number of dt :"+str(len(dt)))
dg=[]
filename2 = open(r'C:\Users\11445\Desktop\disease_genes_col.txt')
for dt2 in filename2:
    dt2 = str(dt2).replace('\n', '')
    dg.append(dt2)
    ddt.append(dt2)
print("the number of dg :"+str((len(dg))))#the number of disease genes

# input PPI
Gint = nx.read_edgelist(r"C:\Users\11445\Desktop\PPI_feng.txt", create_using=nx.DiGraph())
dt_inG = list(np.intersect1d(dt, Gint.nodes()))
G_sub_temp = nx.Graph(nx.subgraph(Gint, dt_inG))

#Disease risk gene networks (constructing disease modules) or overlapping networks of herbicide target networks and disease risk modules (constructing overlapping networks)
focal_genes = list(np.intersect1d(ddt, Gint.nodes())) # only use ddt which are in Gint
G_sub_temp = nx.Graph(nx.subgraph(Gint, focal_genes))
##The number of risk genes in the disease module (in the construction of the disease module) or the number of nodes in the overlapping network (in the construction of the overlapping network)
print ("the number of nodes in ddt :"+str(G_sub_temp.number_of_nodes()))
#The number of edges between risk genes in the disease module (disease module construction) or the number of edges between nodes in the overlapping network (overlapping network construction)
print ("the number of edge in ddt :"+str(G_sub_temp.number_of_edges()))
nx.draw(G_sub_temp, with_labels=True,node_color='r')
plt.show()
network1=[]
network1=nx.to_edgelist(G_sub_temp)
s = str(network1).replace('{','').replace('}','')
s = s.replace("'",'')
f=open(r'DM-PPI.txt','a')
f.write(str(s)+'\n')
f.close()

###########The largest connected component of the disease module or overlapping network
network=[]
network=nx.to_edgelist(G_sub_temp)
a=focal_genes
LCC_MAX =[]
G_sub_LCC =nx.connected_components(G_sub_temp)
LCC_MAX = max(nx.connected_components(G_sub_temp), key=len)
print("the number of lccmax :"+str(len(LCC_MAX)))#Check the size of S
Smax= str(LCC_MAX).replace('{','').replace('}','')
Smax = Smax.replace(",",'\n').replace('','')
f=open(r'DM-Smax1.txt','a')
f.write(str(Smax)+'\n')
f.close()

#Derive the node with the largest connected component S=1
LCC_list=[]
LCC_S1=[]
for i in G_sub_LCC:
    LCC_list.append(len(i))
    if len(i)==1:
        LCC_S1.append(i)
S1= str(LCC_S1).replace('{','').replace('}','')
S1 = S1.replace(",",'').replace("'",'')
f=open(r'DM-S1.txt','a')
f.write(str(S1)+'\n')
f.close()
print("the number of LCCS1 :"+str(len(LCC_S1)))#Check the size of LCC_S1

#Calculate the shortest distance between risk genes in the disease module (construct the disease module) or calculate the shortest distance from disease risk genes to the herbicide targets (overlapping network)
result=[]
for x in dg:
    count=[]
    node_dginds=[]
    
    if x in a:
        result.append(x)
        
        for y in dt:
            node_dginds.append(y)
            ds=[]
            ds1li=[]
            ds2li=[]
            ds3li=[]
            ds4li=[]
            if y in a:
                if x==y:
                    ds0=0
                    ds.append(ds0)
                else:
                    f1=nx.has_path(Gint, source=x, target=y)
                    f2=nx.has_path(Gint, source=y, target=x)
                    f3=nx.has_path(G_sub_temp, source=x, target=y)
                    if f3 == True:
                        ds4=nx.shortest_path_length(G_sub_temp,source=x, target=y)
                        p4=nx.shortest_path(G_sub_temp,source=x, target=y)
                        ds.append(ds4)
                        ds4li.append(ds4)
                    if f1 == True:
                        ds1=nx.shortest_path_length(Gint,source=x, target=y)
                        p1=nx.shortest_path(Gint,source=x, target=y)
                        ds.append(ds1)
                    if f2==True:
                        ds2=nx.shortest_path_length(Gint,source=y, target=x)
                        p2=nx.shortest_path(Gint,source=y, target=x)
                        ds.append(ds2)                       
                    else:                       
                        ds3=999   
                        ds.append(ds3)
                ds=min(ds) 
                count.append(ds)  
        result.append(count)
        result.append(min(count))
        result.append("/n")
        
result=str(result).replace("]",'').replace('[','')       
result=result.replace("'/n',",'\n').replace("'",'')
result=result.replace(' ','').replace(',','\t')
result=result.replace("/n",'\n')
#print(result)         
result=result.replace("'/n',",'\n').replace(',','\t')
result=result.replace("'",'').replace(' ','')
result=result.replace("/n",'\n')
f=open(r'DM_ds.tsv','a')
f.write(str(result)+'\n')
f.close()
print("It's ok")      

    