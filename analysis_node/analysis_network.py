import numpy as np


import networkx as nx
import math
#from networkx.algorithms.community.centrality import girvan_newman
#from networkx.algorithms.community.kclique import k_clique_communities
#from networkx.algorithms.community import greedy_modularity_communities
#from networkx.algorithms.community import label_propagation_communities

import analysis_node



#shannon entropy
#thermodynamic depth
#jensen-shannon divergence

#polytopal complexity

#harmonic geodesic distance
#centralization of stress centrality
#centralization of eigenvector centrality

def powerlaw_fit(matrix):
    import igraph as ig
    import powerlaw
    conn_indices = np.where(matrix)
    # get the weights corresponding to these indices
    weights = matrix[conn_indices]
    # a sequence of (i, j) tuples, each corresponding to an edge from i -> j
    edges = zip(*conn_indices)
    g = ig.Graph(edges)

    d = g.degree()
    result = powerlaw.Fit(d)
    return result



def detect_community(matrix, weight="weight",method="leidenalg"):
    #Find k-clique communities in graph using the percolation method.
    g = nx.from_numpy_matrix(matrix)
    dict_communities = {}
    #if method=="percolation":
        #communities = list(k_clique_communities(g, 3))
    #elif method=="Clauset-Newman-Moore":
        #communities = list(greedy_modularity_communities(g))
    #elif method == "infomap":
        #communities = findCommunities(g)
    #elif method == "label_propagation":
        #communities = label_propagation_communities(g)
    if method == "leidenalg": # run use python2
        import igraph as ig  # run use python2
        import leidenalg
        conn_indices = np.where(matrix)
        # get the weights corresponding to these indices
        weights = matrix[conn_indices]
        # a sequence of (i, j) tuples, each corresponding to an edge from i -> j
        edges = zip(*conn_indices)
        g = ig.Graph(edges)
        communities = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition)
        #print(communities.membership)
        #print(len(communities.membership))
        #print("network adj shape",matrix.shape)
        for i in range(0,len(communities.membership)):
            dict_communities[i] = communities.membership[i]
        return dict_communities

#neumann entropy
def get_HVN(matrix,weight="weight"):
    degree_dict = analysis_node.get_degree_dic(matrix)
    V=matrix.shape[0]
    degree_term_sum = 0
    for i in range(0,V):
        for j in range(0,V):
            if i > j:
                if matrix[i,j] > 0:
                    degree_i = degree_dict.get(i)
                    degree_j = degree_dict.get(j)
                    degree_term_sum = degree_term_sum+(1/(degree_i*degree_j))

    HVN = 1-(1/V)-(1/math.pow(V, 2))*degree_term_sum
    return HVN

#change of neumann entropy
def get_delta_HVN(matrix1,node_list1, matrix2, node_list2, weight="weight"):
    #get node union and only consider node degree not zero
    degree_dict1 = analysis_node.get_degree_dic(matrix1)
    V1 = matrix1.shape[0]
    degree_dict2 = analysis_node.get_degree_dic(matrix2)
    V2 = matrix2.shape[0]
    dHVN = None
    for i in range(0,V1):
        for j in range(0, V1):
            if i > j:
                if matrix1[i, j] > 0:
                    #if this edge in network2
                    pass


#average connectivity
def get_average_degree_connectivity(matrix,weight="weight"):
    #np.fill_diagonal(matrix, 0)
    g = nx.from_numpy_matrix(matrix)
    return nx.average_degree_connectivity(g, source='in+out', target='in+out', nodes=None, weight=None)  # nx


#average geodesic distance
#or average shortest path
def get_average_geodesic_distance(matrix, weight="weight"):
    # np.fill_diagonal(matrix, 0)
    g = nx.from_numpy_matrix(matrix)
    return nx.average_shortest_path_length(g)


#geodesic efficiency
#or the average global efficiency of the graph.
#The efficiency of a pair of nodes in a graph is the multiplicative inverse of the shortest path distance between the nodes.
def get_geodesic_efficiency(matrix, weight="weight"):
    # np.fill_diagonal(matrix, 0)
    g = nx.from_numpy_matrix(matrix)
    return nx.global_efficiency(g)


#centralization of degree
def get_centralization_of_degree(matrix, weight="weight"):
    # np.fill_diagonal(matrix, 0)
    g = nx.from_numpy_matrix(matrix)
    N=g.order()
    indegrees = g.in_degree().values()
    max_in = max(indegrees)
    centralization = float((N*max_in - sum(indegrees)))/(N-1)**2
    return centralization


#density
def get_density(matrix, weight="weight"):
    # np.fill_diagonal(matrix, 0)
    g = nx.from_numpy_matrix(matrix)
    return nx.density(g)


#avearge cluster coefficient
def get_average_cluster_coefficient(matrix, weight="weight"):
    # np.fill_diagonal(matrix, 0)
    g = nx.from_numpy_matrix(matrix)
    return nx.average_clustering(g, nodes=None, weight=None, count_zeros=True)

#transtivity
def get_transtivity(matrix,weight="weight"):
    # np.fill_diagonal(matrix, 0)
    g = nx.from_numpy_matrix(matrix)
    return nx.transitivity(g)

#connectedness
def get_graph_connectivity(matrix, weight="weight"):
    # np.fill_diagonal(matrix, 0)
    g = nx.from_numpy_matrix(matrix)
    return nx.node_connectivity(g, s=None, t=None, flow_func=None)

#following graph as igraph

def scipy_to_igraph(matrix, directed=True):
    sources, targets = matrix.nonzero()
    weights = matrix[sources, targets]
    return ig.Graph(zip(sources, targets), directed=directed, edge_attrs={'weight': weights})


def assortativity(graph, degrees=None):
    if degrees is None: degrees = graph.degree()
    degrees_sq = [deg ** 2 for deg in degrees]

    m = float(graph.ecount())
    num1, num2, den1 = 0, 0, 0
    for source, target in graph.get_edgelist():
        num1 += degrees[source] * degrees[target]
        num2 += degrees[source] + degrees[target]
        den1 += degrees_sq[source] + degrees_sq[target]

    num1 /= m
    den1 /= 2 * m
    num2 = (num2 / (2 * m)) ** 2

    return (num1 - num2) / (den1 - num2)

def betweenness_centralization(G):
    vnum = G.vcount()
    if vnum < 3:
        raise ValueError("graph must have at least three vertices")
    denom = (vnum - 1) * (vnum - 2)

    temparr = [2 * i / denom for i in G.betweenness()]
    max_temparr = max(temparr)
    return sum(max_temparr - i for i in temparr) / (vnum - 1)

def laplacian_centrality(graph, vs=None):
    if vs is None:
        vs = ig.xrange(graph.vcount())
    degrees = graph.degree(mode="all")
    result = []
    for v in vs:
        neis = graph.neighbors(v, mode="all")
        result.append(degrees[v]**2 + degrees[v] + 2 * sum(degrees[i] for i in neis))
    return result

#centralization of betweeness
def get_centralization_of_betweeness(matrix, weight="weight"):
    g = scipy_to_igraph(matrix, False)
    return betweenness_centralization((g))
