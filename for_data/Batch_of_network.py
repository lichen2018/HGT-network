from for_data import Network
from analysis_node import analysis_batch_of_network
import numpy as np
import pandas as pd
import os

class Batch_of_network: # need finish #call function in analysis_batch_of_network

    def __init__(self):
        self.id_to_network = {}
        self.id_to_file = {}
        self.file_to_id = {}
        self.id_to_label_table = None #a df
        self.label1_col = None
        self.label2_col = None
        self.label3_col = None
        self.time_col = None
        self.statistic_output_path = None
        self.network_comparison_output_path = None
        self.network_comparison_outfiles_list_file = None
        self.merge_network_comparison_output_path = None
        self.community_output_path = None

    def load_networks_from_file_info_table(self,filename,file_col,network_id_col,label1_col=None,label2_col=None,label3_col=None,time_col=None):
        pass

    def load_file_list_groupInfo(self,file_path):
        df = pd.read_csv(file_path, header='infer',dtype=str)
        if file_path.endswith(".txt"):
            df = pd.read_csv(file_path, header='infer', delimiter="\t", dtype=str)  #
        return df

    def get_node_name(self,fileName):
        if os.path.isfile(fileName):
            if os.stat(filename).st_size == 0:
                return None
            else:
                result_table = pd.read_csv(fileName, header='infer',index_col=0)
                mat = result_table.values
                node_name = mat[:,0]
                return node_name #return node_name_an
        else:
            print("no file"+fileName)
            return None

    def get_key_for_a_group(self, state):
        pass


    def key_to_path(self, key_row):
        pass

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

    def new_batch_of_network(self):
        pass

    def calculate_network_statistics_for_networks(self):
        #and output
        pass

    def calculate_node_properties_for_networks(self):
        #and output
        pass

    def network_comparison_vs_networks_by_group(self):
        #record
        pass

    def merge_network_comparison_out_from_listfile(self):
        # by network_comparison_outfiles_list_file
        pass

    def community_detection_for_networks(self):
        pass

    def community_evolution_for_networks(self):
        pass
