
from sklearn import preprocessing
from scipy import stats
from analysis_node import analysis_node
from analysis_node import analysis_network
import heapq


class Network:

    def __init__(self, adjacency_matrix, name, node_name):
        self.name = name
        self.adjacency_matrix = adjacency_matrix
        self.distance_matrix = None
        self.similarity_name = None
        self.distance_name = None
        if isinstance(node_name, list):
            self.node_name = node_name  # list
        else:
            self.node_name = node_name.tolist()
        self.degree_dic = None
        self.betweenness_centrality_dic = None
        self.clustering_coefficient_dic = None
        self.eigenvector_centrality_dic = None
        self.pagerank_dic = None
        self.community_color_dic = None

        self.normalized_degree_dic = None
        self.normalized_betweenness_centrality_dic = None
        self.normalized_clustering_coefficient_dic = None
        self.normalized_eigenvector_centrality_dic = None
        self.normalized_pagerank_dic = None

        self.density = None
        self.HVN = None
        self.average_degree_connectivity = None
        self.average_geodesic_distance = None
        self.geodesic_efficiency = None
        self.centralization_of_degree = None
        self.average_cluster_coefficient = None
        self.transtivity = None
        self.graph_connectivity = None

    def get_distance_matrix(self):
        return self.distance_matrix

    def put_distance_matrix(self, distance_matrix):
        self.distance_matrix = distance_matrix

    def get_adjacency_matrix(self):
        return self.adjacency_matrix

    def put_adjacency_matrix(self, adjacency_matrix):
        self.adjacency_matrix = adjacency_matrix

    def set_similarity_name(self,similarity_name):
        self.similarity_name = similarity_name

    def get_similarity_name(self):
        return self.similarity_name

    def set_distance_name(self, distance_name):
        self.distance_name = distance_name

    def get_distance_name(self):
        return self.distance_name

    def get_network_name(self):
        return self.name

    def set_node_name(self, node_name):
        self.node_name = node_name

    def get_node_name(self):
        return self.node_name

    def get_node_index_by_name(self, node_name):
        #index = np.where(self.node_name == node_name)
        index = -1
        if(node_name in self.node_name):
            index = self.node_name.index(node_name)
        return index

    def get_top_node_by_property(self, node_property="pagerank", top_percent=5):
        initial_dict = None
        top_dict = {}
        top_num = int(round(self.adjacency_matrix.shape[0]*(top_percent/100)))
        if node_property == "degree":
            self.degree_dic = analysis_node.get_degree_dic(self.adjacency_matrix)
            initial_dict = self.degree_dic
        elif node_property == "betweenness_centrality":
            self.betweenness_centrality_dic = analysis_node.get_betweenness_centrality_dic(self.adjacency_matrix)
            initial_dict = self.betweenness_centrality_dic
        elif node_property == "clustering_coefficient":
            self.clustering_coefficient_dic = analysis_node.get_clustering_coefficient_dic(self.adjacency_matrix)
            initial_dict = self.clustering_coefficient_dic
        elif node_property == "eigenvector_centrality":
            self.eigenvector_centrality_dic = analysis_node.get_eigenvector_centrality_dic(self.adjacency_matrix)
            initial_dict = self.eigenvector_centrality_dic
        elif node_property == "pagerank":
            self.pagerank_dic = analysis_node.get_pagerank_dic(self.adjacency_matrix)
            initial_dict = self.pagerank_dic
        if initial_dict is not None:
            nodes = heapq.nlargest(top_num, initial_dict, key=initial_dict.get)
            for val in nodes:
                top_dict.update({self.node_name[val]: initial_dict.get(val)})
            return top_dict
        else:
            return None

    def get_node_property(self,node_property="pagerank"):
        if node_property == "degree":
            self.degree_dic = analysis_node.get_degree_dic(self.adjacency_matrix)
            return self.degree_dic
        elif node_property == "betweenness_centrality":
            self.betweenness_centrality_dic = analysis_node.get_betweenness_centrality_dic(self.adjacency_matrix)
            return self.betweenness_centrality_dic
        elif node_property == "clustering_coefficient":
            self.clustering_coefficient_dic = analysis_node.get_clustering_coefficient_dic(self.adjacency_matrix)
            return self.clustering_coefficient_dic
        elif node_property == "eigenvector_centrality":
            self.eigenvector_centrality_dic = analysis_node.get_eigenvector_centrality_dic(self.adjacency_matrix)
            return self.eigenvector_centrality_dic
        elif node_property == "pagerank":
            self.pagerank_dic = analysis_node.get_pagerank_dic(self.adjacency_matrix)
            return self.pagerank_dic

    def calculate_node_topology_property(self, node_property="pagerank"):
        if node_property == "degree":
            self.degree_dic = analysis_node.get_degree_dic(self.adjacency_matrix)
        elif node_property == "betweenness_centrality":
            self.betweenness_centrality_dic = analysis_node.get_betweenness_centrality_dic(self.adjacency_matrix)
        elif node_property == "clustering_coefficient":
            self.clustering_coefficient_dic = analysis_node.get_clustering_coefficient_dic(self.adjacency_matrix)
        elif node_property == "eigenvector_centrality":
            self.eigenvector_centrality_dic = analysis_node.get_eigenvector_centrality_dic(self.adjacency_matrix)
        elif node_property == "pagerank":
            self.pagerank_dic = analysis_node.get_pagerank_dic(self.adjacency_matrix)

    def calculate_node_topology_properties_all(self):
        print("calculating")
        print("degree")
        self.degree_dic = analysis_node.get_degree_dic(self.adjacency_matrix)
        print("betweeness")
        self.betweenness_centrality_dic = analysis_node.get_betweenness_centrality_dic(self.adjacency_matrix)
        print("clustering coefficient")
        self.clustering_coefficient_dic = analysis_node.get_clustering_coefficient_dic(self.adjacency_matrix)
        print("eigenvector centrality")
        self.eigenvector_centrality_dic = analysis_node.get_eigenvector_centrality_dic(self.adjacency_matrix)
        print("pagerank")
        self.pagerank_dic = analysis_node.get_pagerank_dic(self.adjacency_matrix)
        print("end calculating")

    def normalize_node_topology_properties_all(self, normalize_method="min_max"):
        if normalize_method == "zscore":
            if self.degree_dic is not None:
                keys, vals = zip(*self.degree_dic.items())
                self.normalized_degree_dic = dict(zip(keys, stats.zscore(vals, ddof=1)))
            if self.betweenness_centrality_dic is not None:
                keys, vals = zip(*self.betweenness_centrality_dic.items())
                self.normalized_betweenness_centrality_dic = dict(zip(keys, stats.zscore(vals, ddof=1)))
            if self.clustering_coefficient_dic is not None:
                keys, vals = zip(*self.clustering_coefficient_dic.items())
                self.normalized_clustering_coefficient_dic = dict(zip(keys, stats.zscore(vals, ddof=1)))
            if self.eigenvector_centrality_dic is not None:
                keys, vals = zip(*self.eigenvector_centrality_dic.items())
                self.normalized_eigenvector_centrality_dic = dict(zip(keys, stats.zscore(vals, ddof=1)))
            if self.pagerank_dic is not None:
                keys, vals = zip(*self.pagerank_dic.items())
                self.normalized_pagerank_dic = dict(zip(keys, stats.zscore(vals, ddof=1)))
        elif normalize_method == "min_max":
            if self.degree_dic is not None:
                keys, vals = zip(*self.degree_dic.items())
                self.normalized_degree_dic = dict(zip(keys, preprocessing.minmax_scale(vals)))
            if self.betweenness_centrality_dic is not None:
                keys, vals = zip(*self.betweenness_centrality_dic.items())
                self.normalized_betweenness_centrality_dic = dict(zip(keys, preprocessing.minmax_scale(vals)))
            if self.clustering_coefficient_dic is not None:
                keys, vals = zip(*self.clustering_coefficient_dic.items())
                self.normalized_clustering_coefficient_dic = dict(zip(keys, preprocessing.minmax_scale(vals)))
            if self.eigenvector_centrality_dic is not None:
                keys, vals = zip(*self.eigenvector_centrality_dic.items())
                self.normalized_eigenvector_centrality_dic = dict(zip(keys, preprocessing.minmax_scale(vals)))
            if self.pagerank_dic is not None:
                keys, vals = zip(*self.pagerank_dic.items())
                self.normalized_pagerank_dic = dict(zip(keys, preprocessing.minmax_scale(vals)))

    def get_node_property_by_index(self, node_index, node_property):
        value = None
        if node_index != -1:
            if node_property == "degree":
                value = self.degree_dic.get(node_index)
            elif node_property == "betweenness_centrality":
                value = self.betweenness_centrality_dic.get(node_index)
            elif node_property == "clustering_coefficient":
                value = self.clustering_coefficient_dic.get(node_index)
            elif node_property == "eigenvector_centrality":
                value = self.eigenvector_centrality_dic.get(node_index)
            elif node_property == "pagerank":
                value = self.pagerank_dic.get(node_index)
            elif node_property == "community_color":
                value = self.community_color_dic.get(node_index)
        return value

    def get_node_property_by_name(self, node_name, node_property="pagerank"):
        index = self.get_node_index_by_name(node_name)
        if index < 0:
            return None
        else:
            value = self.get_node_property_by_index(index, node_property)
            return value

    def get_node_property_normalized_by_index(self, node_index, node_property):  # need modify
        value = None
        if node_index != -1:
            if node_property == "degree":
                value = self.normalized_degree_dic.get(node_index)
            elif node_property == "betweenness_centrality":
                value = self.normalized_betweenness_centrality_dic.get(node_index)
            elif node_property == "clustering_coefficient":
                value = self.normalized_clustering_coefficient_dic.get(node_index)
            elif node_property == "eigenvector_centrality":
                value = self.normalized_eigenvector_centrality_dic.get(node_index)
            elif node_property == "pagerank":
                value = self.normalized_pagerank_dic.get(node_index)
        return value

    def get_node_property_normalized_by_name(self, node_name, node_property="pagerank"):
        index = self.get_node_index_by_name(node_name)
        value = self.get_node_property_normalized_by_index(index, node_property)
        return value

    def calculate_network_statistics_all(self, property="density"):
        if property == "density":
            self.density = analysis_network.get_density(self.adjacency_matrix)
        elif property == "HVN":
            self.HVN = analysis_network.get_HVN(self.adjacency_matrix)
        elif property == "average_degree_connectivity":
            self.average_degree_connectivity = analysis_network.get_average_degree_connectivity(self.adjacency_matrix)
        elif property == "average_geodesic_distance":
            self.average_geodesic_distance = analysis_network.get_average_geodesic_distance(self.adjacency_matrix)
        elif property == "geodesic_efficiency":
            self.geodesic_efficiency = analysis_network.get_geodesic_efficiency(self.adjacency_matrix)
        elif property == "centralization_of_degree":
            self.centralization_of_degree = analysis_network.get_centralization_of_degree(self.adjacency_matrix)
        elif property == "average_cluster_coefficient":
            self.average_cluster_coefficient = analysis_network.get_average_cluster_coefficient(self.adjacency_matrix)
        elif property == "transtivity":
            self.transtivity = analysis_network.get_transtivity(self.adjacency_matrix)
        elif property == "graph_connectivity":
            self.graph_connectivity = analysis_network.get_graph_connectivity(self.adjacency_matrix)

    def get_network_statistics_by_property(self, property="density"):
        if self.density is None:
            self.calculate_network_statistics_all()
        if property == "density":
            return self.density
        elif property == "HVN":
            return self.HVN
        elif property == "average_degree_connectivity":
            return self.average_degree_connectivity
        elif property == "average_geodesic_distance":
            return self.average_geodesic_distance
        elif property == "geodesic_efficiency":
            return self.geodesic_efficiency
        elif property == "centralization_of_degree":
            return self.centralization_of_degree
        elif property == "average_cluster_coefficient":
            return self.average_cluster_coefficient
        elif property == "transtivity":
            return self.transtivity
        elif property == "graph_connectivity":
            return self.graph_connectivity

    def detect_community(self,k=None):
        self.community_color_dic = analysis_network.detect_community(self.adjacency_matrix)
        #print(self.community_color_dic)

    def calculate_network_statistics_for_each_community(self):
        #need finish
        pass
