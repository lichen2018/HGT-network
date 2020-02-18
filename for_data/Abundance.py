import pandas as pd
import numpy as np
import os.path


class Abundance:

    def __init__(self, file_path=None, name=None):
        self.name = None
        self.file_path = None
        self.sample_id = None
        self.node_name = None
        self.abundance_matrix = None
        if file_path is not None:
            if os.path.isfile(file_path):
                self.name = name
                self.file_path = file_path
                df = pd.read_csv(file_path, header='infer')
                if self.file_path.endswith(".txt"):
                    df = pd.read_csv(file_path, header='infer', delimiter="\t")
                #m1 = df.as_matrix()
                m1 = df.values #test
                sample_id = df.columns.values
                self.sample_id = sample_id[1:len(sample_id)]
                self.node_name = m1[:, 0]
                self.abundance_matrix = np.delete(m1, 0, axis=1)
            else:
                print("Do not exist file " + file_path + " ,please check")


    def new_abundance(self, name, sample_id, node_name, abundance_matrix):
        self.name = name
        self.sample_id = sample_id
        self.node_name = node_name
        self.abundance_matrix = abundance_matrix

    def select_according_to_node_order(self, node_name):
        nodelist1 = [a.lower() for a in node_name]
        nodelist2 = [a.lower() for a in self.node_name]
        node_name = np.asarray(nodelist1)
        node_name_in_abundance = np.array(nodelist2)
        select = []
        for i in range(0, len(node_name)):
            index = np.where(node_name_in_abundance == node_name[i])
            select.append(index)
        abundance_matrix = self.abundance_matrix[select,:]
        return abundance_matrix

    def reorder_according_to_node_order(self, node_name):
        nodelist1 = [a.lower() for a in node_name]
        nodelist2 = [a.lower() for a in self.node_name]
        node_name = np.asarray(nodelist1)
        node_name_in_abundance = np.array(nodelist2)
        reorderlist = []
        for i in range(0, len(node_name)):
            index = np.where(np.asarray(node_name_in_abundance) == node_name[i])
            reorderlist.extend(index[0])
        self.node_name = [self.node_name[j] for j in reorderlist]
        self.abundance_matrix = self.abundance_matrix[reorderlist, :]

    def reorder_according_to_sample_order(self, sample_id):
        sample_id = np.asarray(sample_id)
        sample_id_in_abundance = np.array(self.sample_id)
        if sample_id.dtype is not str:
            sample_id_in_abundance = np.asarray(sample_id_in_abundance, dtype=sample_id.dtype)
        if len(self.sample_id) != len(sample_id):
            print("sample ID in abundance table and sampleInfo may not be one to one map")
        reorderlist = []
        for i in range(0, len(sample_id)):
            index = np.where(np.asarray(sample_id_in_abundance) == sample_id[i])
            reorderlist.extend(index[0])
        self.sample_id = [self.sample_id[j] for j in reorderlist]
        self.abundance_matrix = self.abundance_matrix[:, reorderlist]

    def select_samples_and_return_abundance_matrix(self, sample_id):
        sample_id = np.asarray(sample_id)
        sample_id_in_abundance = np.array(self.sample_id)
        select = []
        if sample_id.dtype is not str:
            sample_id_in_abundance = np.asarray(sample_id_in_abundance, dtype=sample_id.dtype)
        for i in range(0, len(sample_id)):
            index = np.where(sample_id_in_abundance == sample_id[i])
            select.append(index)
        abundance_matrix = self.abundance_matrix[:, select]
        return abundance_matrix

    def __del__(self):
        print("delete Abundance instance "+str(self.name))

    def get_network_name(self):
        return self.name

    def get_abundance_matrix(self):
        return self.abundance_matrix

    def get_node_name(self):
        return self.node_name

    def get_sample_id(self):
        return self.sample_id

    def set_network_name(self, name):
        self.name = name

    def set_abundance_matrix(self, matrix):
        self.abundance_matrix = matrix

    def set_node_name(self, node_name):
        self.node_name = node_name

    def set_sample_id(self, sample_id):
        self.sample_id = sample_id

    def output_abundance(self,file):
        self.abundance_matrix.astype(float)
        df = pd.DataFrame(self.abundance_matrix, index=self.node_name, columns=self.sample_id)
        df.to_csv(file, sep=',', header=True, index=True)



def test():
    test_abundance = "/Users/jiaxchen2/Desktop/v1Conserve/v1Conserve/testInput/rpkm_L.csv"
    abundance_1 = Abundance(test_abundance, "L")
    print(abundance_1.get_sample_id())
    print(abundance_1.get_abundance_matrix())


if __name__ == '__main__':
    test()
