
import numpy as np
import networkx as nx


# connectivity
# stress centrality
# betweenness
# vulnerability

#harmonic centrality, return dict
#harmonic_centrality(G, nbunch=None, distance=None)
#closeness centrality (networkx



def get_degree_each(matrix, weight="weight"):
    #g = nx.from_numpy_matrix(matrix)
    #ans = g.degree(weight=weight)
    weighted_degree = matrix.sum(axis=1, dtype='float')
    #weighted_degree -= 1 #minus the 1 for node itself
    return weighted_degree


def get_degree_dic(matrix, weight="weight"):
    # remove the 1 for loop of node itself
    np.fill_diagonal(matrix, 0)
    g = nx.Graph(matrix)
    degree = nx.degree(g, weight=weight)
    return dict(degree)


def get_betweenness_centrality_dic(matrix, weight="weight"):
    # remove the 1 for loop of node itself
    np.fill_diagonal(matrix, 0)
    g = nx.from_numpy_matrix(matrix)
    betweenness_centrality = nx.betweenness_centrality(g, weight=weight)
    return dict(betweenness_centrality)


def get_clustering_coefficient_dic(matrix, weight="weight"):
    # remove the 1 for loop of node itself
    np.fill_diagonal(matrix, 0)
    g = nx.from_numpy_matrix(matrix)
    clustering = nx.clustering(g, weight=weight)
    return dict(clustering)


def get_eigenvector_centrality_dic(matrix, weight="weight"):
    # remove the 1 for loop of node itself
    np.fill_diagonal(matrix, 0)
    g = nx.from_numpy_matrix(matrix)
    eigenvector_centrality = nx.eigenvector_centrality(g, weight=weight)
    return dict(eigenvector_centrality)


def get_pagerank_dic(matrix, weight="weight"):
    # remove the 1 for loop of node itself
    np.fill_diagonal(matrix, 0)
    g = nx.from_numpy_matrix(matrix)
    try:
        pagerank = nx.pagerank(g, weight=weight)
        return dict(pagerank)
    except Exception as e:
        print(e)

