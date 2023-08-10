import pandas as pd
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
from collections import Counter
import networkx as nx



def return_degree_distribution(dataframe, return_hub = False, return_bipartite = False):
   

    bipartite_g = nx.Graph()
    bipartite_g.add_nodes_from(dataframe.TF,bipartite = 0)
    bipartite_g.add_nodes_from(dataframe.Genes,bipartite = 1)
    bipartite_g.add_weighted_edges_from([(row['TF'],row['Genes'],row['Score'])for idx, row in dataframe.iterrows()], weight = 'weight')
    top_nodes = {n for n, d in bipartite_g.nodes(data=True) if d["bipartite"] == 0}   # TF
    bottom_nodes = set(bipartite_g) - top_nodes #Genes

    if return_bipartite:
        return bipartite_g, top_nodes,bottom_nodes


    deg_genes, deg_tf = bipartite.degrees(bipartite_g,top_nodes)
    deg_genes_list = []


    for i in deg_genes:
        deg_genes_list.append(i[1])
    count = Counter(deg_genes_list)
    genes_degree_df = pd.DataFrame.from_dict(count, orient='index')

    deg_TF_list = []
    for i in deg_tf:
        deg_TF_list.append(i[1])
    count = Counter(deg_TF_list)
    tf_degree_df = pd.DataFrame.from_dict(count, orient='index')
    if return_hub:
        return deg_genes, deg_tf
    
    return genes_degree_df, tf_degree_df


def return_hub_genes_tf(dataframe):

    genes_degree, tf_degree = return_degree_distribution(dataframe, return_hub= True)
    genes_degree_df = pd.DataFrame(sorted(dict(genes_degree).items(), key=lambda x:x[1], reverse = True)[0:100], columns=['Genes','Degree'])
    genes_degree_df['Index'] = range(1,len(genes_degree_df)+1)
    tf_degree_df = pd.DataFrame(sorted(dict(tf_degree).items(), key=lambda x:x[1], reverse = True)[0:100], columns=['TF','Degree'])
    tf_degree_df['Index'] = range(1,len(tf_degree_df)+1)
    return genes_degree_df, tf_degree_df


def return_connected_components(dataframe):

    bipartite_g, top_nodes, bottom_nodes = return_degree_distribution(dataframe, return_hub = False, return_bipartite = True)

    dct = {} 
    for i,c in enumerate(sorted(nx.connected_components(bipartite_g), key=len, reverse=True)):
        dct[i] = str(c)[1:-1]
    
  
    return pd.DataFrame.from_dict(dct, orient = 'index')



def return_tf_projection(dataframe, projection = 'TF'):

    bipartite_g, top_nodes, bottom_nodes = return_degree_distribution(dataframe, return_hub = False, return_bipartite = True)
    
    if projection == 'TF':
        TF_projection = bipartite.projected_graph(bipartite_g, top_nodes)
        return TF_projection
    else:
        Gene_projection = bipartite.projected_graph(bipartite_g, bottom_nodes)
        return Gene_projection



def return_degree_distribution_graph_input(graph, return_hub= False):

    deg_list = []
    if return_hub:
        return graph.degree
    for i in (graph.degree):
        if i[1] == 0:
            continue
        deg_list.append(i[1])

    return pd.DataFrame.from_dict(Counter(deg_list), orient='index')


def return_hub_tf(graph):

    degree = return_degree_distribution_graph_input(graph, return_hub= True)
    
    degree_df = pd.DataFrame(sorted(dict(degree).items(), key=lambda x:x[1], reverse = True)[0:100], columns=['TF','Degree'])
    degree_df['Index'] = range(1,len(degree_df)+1)
    return degree_df

def return_hub_genes(graph):

    degree = return_degree_distribution_graph_input(graph, return_hub= True)
    
    degree_df = pd.DataFrame(sorted(dict(degree).items(), key=lambda x:x[1], reverse = True)[0:100], columns=['Genes','Degree'])
    degree_df['Index'] = range(1,len(degree_df)+1)
    return degree_df


def return_connected_components_projection(graph):

    dct = {} 
    for i,c in enumerate(sorted(nx.connected_components(graph), key=len, reverse=True)):
        dct[i] = str(c)[1:-1]
    
  
    return pd.DataFrame.from_dict(dct, orient = 'index')



