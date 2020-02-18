import argparse

import sys
import os
import pandas as pd
import numpy as np

from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial import distance as ssd

from analysis_node import analysis_network
from analysis_node import analysis_node
from for_data import Network


class batch_network_analysis:
    def __init__(self, groupInfo_file, output_dir, network_col=0, network_id_col=1):
        self.file_list_groupInfo_file = groupInfo_file
        self.file_list_groupInfo = None
        self.output_path = output_dir
        self.network_statistics_table = None
        self.network_id_to_statistics = {}
        #for main_step3
        self.network_id_to_Network = {}
        self.species_pool = []
        self.unique_species_pool = []
        self.network_id_pool = []
        #for main_step4
        self.species_similarity = None
        self.species_distance = None
        self.aggregated_species_network = None
        #for main_step5
        self.community_color_table = None
        #input file setting for file_list_groupInfo_file
        self.network_col = network_col
        self.network_id_col = network_id_col

    def load_file_list_groupInfo(self,file_path):
        df = pd.read_csv(file_path, header='infer',dtype=str)
        if file_path.endswith(".txt"):
            df = pd.read_csv(file_path, header='infer', delimiter="\t", dtype=str)  #
        return df

    def get_node_name(self,fileName):
        if os.path.isfile(fileName):
            if os.stat(fileName).st_size == 0:
                return None
            else:
                result_table = pd.read_csv(fileName, header='infer')
                mat = result_table.values
                node_name = mat[:,0]
                return node_name #return node_name_an
        else:
            print("no file"+fileName)
            return None

    def load_sparse_network(self,filename):
        print(filename)
        if os.path.isfile(filename):
            if os.stat(filename).st_size == 0:
                print(filename+" has no content")
                return None
            else:
                m = pd.read_csv(filename, header=None)
                adj = m.values
                where_are_NaNs = np.isnan(adj)
                adj[where_are_NaNs] = 0
                return adj
        else:
            print("no file"+filename)
            return None

    def load_sparse_network_withHeader(self,filename):
        if os.path.isfile(filename):
            if os.stat(filename).st_size == 0:
                print(filename+" has no content")
                return None
            else:
                adj1 = pd.read_csv(filename, header='infer')
                nodes1 = np.asarray(adj1.iloc[:, 0].astype(str))  #
                adj1.fillna(value=0)
                m1 = adj1.values
                matrix1 = m1[:, 1:m1.shape[1]]
                matrix1 = np.nan_to_num(matrix1)

                matrix1 = pd.DataFrame(matrix1)
                node_list1 = pd.DataFrame(nodes1)

                matrix1.to_csv(filename+"_adj.csv", header=False, index=False)
                node_list1.to_csv(filename+"_nodename.csv", header=False, index=False)
                matrix1 = self.load_sparse_network(filename+"_adj.csv")
                return matrix1
        else:
            print("no file"+filename)
            return None

    def load_node_name(self,filename):
        if os.path.isfile(filename):
            if os.stat(filename).st_size == 0:
                print(filename+" has no content")
                return None
            else:
                adj1 = pd.read_csv(filename, header='infer')
                nodes1 = np.asarray(adj1.iloc[:, 0].astype(str))  #
                node_list1 = pd.DataFrame(nodes1)
                return node_list1
        else:
            print("no file"+filename)
            return None

    def calculate_network_statistics_for_one_network(self,network_adj_file,network_id):
        print(network_adj_file)
        matrix = self.load_sparse_network_withHeader(network_adj_file)
        statistics = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
        if matrix is not None:
            print(matrix.shape)
            if True:
            #if matrix.shape[0] >= 4000:
                HVN = analysis_network.get_HVN(matrix)
                average_cluster_coefficient = analysis_network.get_average_cluster_coefficient(matrix)
                #print(analysis_network.get_average_geodesic_distance(matrix))
                geodesic_efficiency = analysis_network.get_geodesic_efficiency(matrix)
                #print(analysis_network.get_centralization_of_degree(matrix))
                density = analysis_network.get_density(matrix)
                transtivity = analysis_network.get_transtivity(matrix)
                graph_connectivity = analysis_network.get_graph_connectivity(matrix)
                statistics = [matrix.shape[0], HVN, average_cluster_coefficient, geodesic_efficiency, density, transtivity, graph_connectivity]
        self.network_id_to_statistics[network_id] = statistics
        return statistics


    def calculate_network_statistics_for_list_network(self, out_file):
        self.file_list_groupInfo = self.load_file_list_groupInfo(self.file_list_groupInfo_file)
        f = open(out_file, 'a')
        f.write("network_id,size,HVN,average_cluster_coefficient,geodesic_efficiency,density,transtivity,graph_connectivity\n")
        f.close()
        for i in range(0, self.file_list_groupInfo.shape[0]):
            f = open(out_file, 'a')
            network_file = self.file_list_groupInfo.iloc[i][0]
            if self.network_id_col is None:
                print("require network id in the groupInfo table")
            network_id = self.file_list_groupInfo.iloc[i][self.network_id_col]
            statistics = self.calculate_network_statistics_for_one_network(network_file,network_id)
            print(statistics)
            f.write(network_id+","+str(statistics[0])+","+str(statistics[1])+","+str(statistics[2])+","+str(statistics[3])+","+str(statistics[4])+","+str(statistics[5])+","+str(statistics[6])+"\n")
            f.close()

    def output_network_node_statistics_for_network_list(self,out_file, node_property="degree"):
        result_mat = np.zeros((len(self.unique_species_pool),len(self.network_id_pool)))
        # one table for one property (all nodes * network_id, list in time order)
        for i in range(0, len(self.network_id_pool)):
            network_i = self.network_id_to_Network.get(self.network_id_pool[i])
            for j in range(0, len(self.unique_species_pool)):
                result_mat[j, i] = network_i.get_node_property_by_name(self.unique_species_pool[j],node_property)
        df = pd.DataFrame(result_mat,dtype=float,columns=self.network_id_pool,index=self.unique_species_pool)
        df.to_csv(out_file)
        return result_mat

    def calculate_network_node_statistics_for_list_network(self, out_file_header):
        self.file_list_groupInfo = self.load_file_list_groupInfo(self.file_list_groupInfo_file)
        for i in range(0, self.file_list_groupInfo.shape[0]):
            network_file = self.file_list_groupInfo.iloc[i][self.network_col]
            if self.network_id_col is not None:
                network_id = self.file_list_groupInfo.iloc[i][self.network_id_col]
            else:
                print("require network id")
            matrix = self.load_sparse_network_withHeader(network_file)
            node_name = self.get_node_name(network_file)
            if matrix is not None:
                network_i = Network.Network(matrix, network_id, node_name)
                network_i.calculate_node_topology_property("degree")
                self.network_id_to_Network[network_id] = network_i
                self.species_pool.extend(node_name)
                self.network_id_pool.append(network_id)
        self.unique_species_pool = list(set(self.species_pool))
        self.output_network_node_statistics_for_network_list(out_file_header + "_degree.csv", "degree")

        for i in range(0,len(self.network_id_pool)):
            network_i = self.network_id_to_Network.get(self.network_id_pool[i])
            network_i.calculate_node_topology_property("pagerank")
            self.network_id_to_Network[self.network_id_pool[i]] = network_i
        self.output_network_node_statistics_for_network_list(out_file_header + "_pagerank.csv", "pagerank")

        for i in range(0,len(self.network_id_pool)):
            network_i = self.network_id_to_Network.get(self.network_id_pool[i])
            network_i.calculate_node_topology_property("clustering_coefficient")
            self.network_id_to_Network[self.network_id_pool[i]] = network_i
        self.output_network_node_statistics_for_network_list(out_file_header + "_clustering_coefficient.csv", "clustering_coefficient")

        #for i in range(0,len(self.network_id_pool)):
        #    network_i = self.network_id_to_Network.get(self.network_id_pool[i])
        #    print(len(network_i.get_node_name()))
        #    network_i.calculate_node_topology_property("eigenvector_centrality")
        #    self.network_id_to_Network[network_id] = network_i
        #self.output_network_node_statistics_for_network_list(out_file_header + "_eigenvector_centrality.csv", "eigenvector_centrality")

        #for i in range(0,len(self.network_id_pool)):
        #    network_i = self.network_id_to_Network.get(self.network_id_pool[i])
        #    network_i.calculate_node_topology_property("betweenness_centrality")
        #    self.network_id_to_Network[network_id] = network_i
        #self.output_network_node_statistics_for_network_list(out_file_header + "_betweenness_centrality.csv", "betweenness_centrality")


    def load_species_network_mat_property_file(self,filename, withHeader=True):
        if os.path.isfile(filename):
            if os.stat(filename).st_size == 0:
                print(filename+" has no content")
                return None
            else:
                if withHeader:
                    mat = pd.read_csv(filename, header='infer',index_col=0)
                else:
                    mat = pd.read_csv(filename, header=None)
                print(mat)
                return mat.T
        else:
            print("no file"+filename)
            return None
        pass


    def get_species_similarity_from_multiple_property_mat(self, species_network_mat_property_file_list,HGT_network_mat=None):
        simi_mat = np.zeros([len(self.unique_species_pool),len(self.unique_species_pool)])
        for i in range(0,len(species_network_mat_property_file_list)):
            mat_dfi = self.load_species_network_mat_property_file(species_network_mat_property_file_list[i])
            mat_corr=mat_dfi.corr()
            print("mat_corr shape:", mat_corr.values.shape)
            simi_mat=simi_mat+np.abs(mat_corr.values)
            np.savetxt(self.output_path+"mat_corr.csv",mat_corr.values,delimiter=",")
            np.savetxt(self.output_path+"simi_mat.csv", simi_mat, delimiter=",")
            #print(mat_corr)
        number_of_property=len(species_network_mat_property_file_list)
        simi_mat=self.normalize_simi_mat_to_max(simi_mat)
        return simi_mat


    def merge_similarity_mat_by_weight(self,simi_mat1,weight1,simi_mat2,weight2):
        simi=simi_mat1*weight1+simi_mat2*weight2
        return simi

    def convert_similarity_to_distance(self,similarity_mat):
        distance_mat = 1-similarity_mat
        np.fill_diagonal(distance_mat,0)
        return distance_mat

    def recover_adj_mat_from_edges_file(self,file):
        adj=np.zeros([len(self.unique_species_pool),len(self.unique_species_pool)])
        species_to_index={}
        for i in range(0,len(self.unique_species_pool)):
            species_to_index[self.unique_species_pool[i]]=i
        with open(file) as fp:
            line = fp.readline()
            cnt = 1
            for line in fp:
                tmp=line.strip().split(',')
                cnt += 1
                species1=tmp[0]
                species2=tmp[1]
                index1=species_to_index.get(species1)
                index2=species_to_index.get(species2)
                adj[index1,index2]=tmp[2]
                adj[index2,index1]=tmp[2]
        return adj

    def normalize_simi_mat_to_max(self,mat):
        [n,m] = mat.shape
        simi_mat=mat
        max_value = mat.max()
        for i in range(0,n):
            simi_mat[i,i]=1
        simi_mat=mat/max_value
        return simi_mat


    def do_cluster_on_species(self,similarity_mat,out_file_header):
        distance_mat = self.convert_similarity_to_distance(similarity_mat)
        # Compare clustering alternatives
        Zd = linkage(ssd.squareform(distance_mat), method="complete")
        # Cluster based on distance matrix
        cld = fcluster(Zd, 10, criterion='distance')
        df = pd.DataFrame(cld,index=self.unique_species_pool)
        df.to_csv(self.output_path+out_file_header+"_cluster.csv")

    def append_network(self,big_network, new_network, normalize=False):
        if big_network is None:
            big_network = np.zeros([len(self.unique_species_pool),len(self.unique_species_pool)])
        if new_network is None:
            return big_network
        adj = new_network.get_adjacency_matrix()
        species_list = new_network.get_node_name()
        #print(species_list)
        if normalize:
            pass
        corresponding_index = []
        for i in range(0,len(species_list)):
            A = self.unique_species_pool.index(str(species_list[i]))
            corresponding_index.append(A)
        for i in range(0,len(corresponding_index)):
            A=corresponding_index[i]
            for j in range(0,len(corresponding_index)):
                B=corresponding_index[j]
                big_network[A,B] = big_network[A,B]+adj[i,j]
        return big_network

    def print_aggregated_network_to_file(self,network,node_names,file_header):
        if network is not None:
            A = np.where(network == 0)
            network[A] = np.nan
            df = pd.DataFrame(network,columns=node_names,index=node_names)
            df.to_csv(file_header+"_mat.csv")
            if False:
                f = open(file_header+"_edges.csv", 'a')
                f.write("source,target,weight\n")
                for i in range(0,network.shape[0]):
                    for j in range(0,network.shape[1]):
                        if network[i,j]!=0:
                            f.write(node_names[i]+","+node_names[j]+","+str(network[i,j])+"\n")
                f.close()

    def aggregate_network_all_for_visualize(self,out_file_header):
        self.init_unique_species_pool_network_pool()
        for i in range(0,len(self.network_id_pool)):
            print(i)
            network_i = self.network_id_to_Network.get(self.network_id_pool[i])
            self.aggregated_species_network=self.append_network(self.aggregated_species_network, network_i, False)
        self.print_aggregated_network_to_file(self.aggregated_species_network,self.unique_species_pool,self.output_path+out_file_header)

    def aggregate_network_for_rows(self,rows,out_file_header):
        family_network_id_pool=[]
        for i in rows:
            if self.network_id_col is not None:
                network_id = self.file_list_groupInfo.iloc[i][self.network_id_col]
            else:
                print("require network id col")
            family_network_id_pool.append(network_id)
        family_aggregated_species_network = None
        for i in range(0,len(family_network_id_pool)):
            network_i = self.network_id_to_Network.get(family_network_id_pool[i])
            if network_i is None:
                print(family_network_id_pool[i], ": network is None")
            else:
                print("aggregating network id:",network_i.get_network_name())
                family_aggregated_species_network=self.append_network(family_aggregated_species_network, network_i, False)
        if family_aggregated_species_network is not None:
            self.print_aggregated_network_to_file(family_aggregated_species_network,self.unique_species_pool,self.output_path+out_file_header)
            family_aggregated_species_network_network = Network.Network(family_aggregated_species_network,"",self.unique_species_pool)
        else:
            family_aggregated_species_network_network = None
        return family_aggregated_species_network_network

    def aggregate_network_for_each_group_and_all(self,out_file_header,group_col):
        self.init_unique_species_pool_network_pool()
        self.file_list_groupInfo = self.load_file_list_groupInfo(self.file_list_groupInfo_file)
        group_labels=list(self.file_list_groupInfo.iloc[:,group_col])
        unique_group_labels = list(set(self.file_list_groupInfo.iloc[:,group_col]))

        #print(group_labels)
        print("unique_group_labels",unique_group_labels)

        self.aggregated_species_network=None
        for i in range(0, len(unique_group_labels)):
            rows_in_group=[]
            for j in range(0,len(group_labels)):
                if group_labels[j]==unique_group_labels[i]:
                    rows_in_group.append(j)
            #print(unique_group_labels[i])
            #print(rows_in_group)
            family_aggregated_species_network = self.aggregate_network_for_rows(rows_in_group, out_file_header+str(unique_group_labels[i]))
            self.aggregated_species_network = self.append_network(self.aggregated_species_network, family_aggregated_species_network, False)
        self.print_aggregated_network_to_file(self.aggregated_species_network, self.unique_species_pool,
                                                  self.output_path + out_file_header)

    def init_unique_species_pool_from_node_property_file(self,fileName,exclude_first_col=False):
        if os.path.isfile(fileName):
            if exclude_first_col:
                result_table = pd.read_csv(fileName, header='infer',index_col=0)
            else:
                result_table = pd.read_csv(fileName, header='infer')
            mat = result_table.values
            node_name = mat[:,0]
            self.unique_species_pool=node_name
            return node_name #return node_name_an
        else:
            print("no file"+fileName)
            return None

    def init_unique_species_pool_from_adj_file(self,fileName):
        if os.path.isfile(fileName):
            self.unique_species_pool=self.load_node_name(fileName)
        else:
            print("no file"+fileName)
            return None

    def init_unique_species_pool_network_pool(self):
        self.aggregated_species_network = None
        self.file_list_groupInfo = self.load_file_list_groupInfo(self.file_list_groupInfo_file)
        #for i in range(0,10):#need finish, for test
        for i in range(0, self.file_list_groupInfo.shape[0]):
            network_file = self.file_list_groupInfo.iloc[i][self.network_col]
            if self.network_id_col is not None:
                network_id = self.file_list_groupInfo.iloc[i][self.network_id_col]
            else:
                print("require network id")
            matrix = self.load_sparse_network_withHeader(network_file)
            if matrix is not None:
                node_name = self.get_node_name(network_file)
                network_i = Network.Network(matrix, network_id, node_name)
                self.network_id_to_Network[network_id] = network_i
                self.species_pool.extend(node_name)
                self.network_id_pool.append(network_id)
        self.unique_species_pool = list(set(self.species_pool))


    def detect_community_for_all_network(self, out_file_header):
        self.file_list_groupInfo = self.load_file_list_groupInfo(self.file_list_groupInfo_file)
        self.init_unique_species_pool_network_pool()
        for i in range(0,len(self.network_id_pool)):
            network_i = self.network_id_to_Network.get(self.network_id_pool[i])
            network_i.detect_community()
        self.community_color_table = self.output_network_node_statistics_for_network_list(out_file_header + "_community_color.csv", "community_color")


    def scale_free_fit_for_one_network(self,network_adj_file,network_id):
        print(network_adj_file)
        matrix = self.load_sparse_network_withHeader(network_adj_file)
        statistics = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
        if matrix is not None:
            print(matrix.shape)
            if True:
            #if matrix.shape[0] >= 4000:
                result = analysis_network.powerlaw_fit(matrix)
                p1 = result.distribution_compare('power_law', 'lognormal_positive')
                p2 = result.distribution_compare('power_law', 'stretched_exponential', normalized_ratio=True)
                p3 = result.distribution_compare('power_law', 'exponential', normalized_ratio=True)
                statistics = [p1[0],p2[0],p3[0],result.xmin,result.power_law.alpha]
        self.network_id_to_statistics[network_id] = statistics
        return statistics


    def scale_free_fit_for_list_network(self, out_file):
        self.file_list_groupInfo = self.load_file_list_groupInfo(self.file_list_groupInfo_file)
        f = open(out_file, 'a')
        f.write("p_powerlaw_vs_lognormal_positive,p_powerlaw_vs_stretched_exponential,p_powerlaw_vs_exponential,xmin,alpha\n")
        f.close()
        for i in range(0, self.file_list_groupInfo.shape[0]):
            f = open(out_file, 'a')
            network_file = self.file_list_groupInfo.iloc[i][0]
            if self.network_id_col is None:
                print("require network id in the groupInfo table")
            network_id = self.file_list_groupInfo.iloc[i][self.network_id_col]
            statistics = self.scale_free_fit_for_one_network(network_file,network_id)
            print(statistics)
            f.write(network_id+","+str(statistics[0])+","+str(statistics[1])+","+str(statistics[2])+","+str(statistics[3])+","+str(statistics[4])+"\n")
            f.close()
