import argparse

import sys
import os
import pandas as pd
import numpy as np

from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial import distance as ssd

import batch_network_analysis
from analysis_node import analysis_network
from analysis_node import analysis_node
from for_data import Network


####################################
##input:    file 1. groupInfo_table 2. a batch of networks (format as adj with..., note that we will transfer na value to 0)
#           parameters: indicate network_col (for file path), network_id_col (for network id)
#                       groupInfo_table_file_path, output_path
#
####################################

parser = argparse.ArgumentParser(description="example")

#input files
parser.add_argument('-step', required=True, help='')
parser.add_argument('-groupInfo_file', default="None", help='')
parser.add_argument('-output_dir',default=".", help='')


parser.add_argument('-network_col',type=int,default=0,help='')
parser.add_argument('-network_id_col',type=int,default=1, help='')
parser.add_argument('-group_col',type=int,default=1, help='')

#for construct network from abundance matrix
parser.add_argument('-abundance_file',  help='')
parser.add_argument('-keep_whole_adj', default="T", help='')#T or F

##for node cluster by property
parser.add_argument('-species_network_mat_property_file_list',help='')
parser.add_argument('-aggregated_network_file',help='')

##for community evolution
parser.add_argument('-community_color_file', default="None", help='')
parser.add_argument('-network_individual_col',default="NA",help='')
parser.add_argument('-network_time_col',default="NA", help='')

##for community cluster
parser.add_argument('-comm_jaccard_simi_mat_file', default="None", help='')
parser.add_argument('-beta',type=int,default=1,help='')
parser.add_argument('-cut_height',type=float,default=0.6, help='')
parser.add_argument('-minClusterSize',type=int,default=10, help='')

##for hgt event cluster
parser.add_argument('-hgt_event_file', default="None", help='')
parser.add_argument('-hgt_event_filter',type=int,default=5, help='')


args = parser.parse_args()


sys.path.append(os.path.realpath('..'))


def main_network_reconstruction(abundanceFile, output_dir,similarityType="correlation",distanceType="root_absolute",keep_whole_adj="T", cutoff=0.3):
    # run function_for_distance.R, with parameter expr file and similarity type, distance type
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if not output_dir.endswith('/'):
        output_dir=str(output_dir)+'/'

    Rscript_other_similarity = "for_reconstruction/run_calculate_distance_by_R.R"
    inputFile = abundanceFile
    output_similarity_file=output_dir+"adj.csv"
    output_distance_file=output_dir+"dist.csv"

    command = "Rscript " + Rscript_other_similarity + " " + inputFile + " " + output_similarity_file + " "\
              + output_distance_file + " " + similarityType + " " + distanceType + " " + keep_whole_adj + " " + str(cutoff)
    print("running script:" + command)
    os.system(command)
    pass


def main_calculate_node_statistics(groupInfo_file,output_dir,network_col=0,network_id_col=1):
    bn = batch_network_analysis.batch_network_analysis(groupInfo_file,output_dir,int(network_col), int(network_id_col))
    # run calculate network node statistics
    #bn.file_list_groupInfo_file = "/Users/jiaxchen2/Desktop/HGT_network/Data_by_Lichen/file_groupInfo_whole.csv"
    #bn.output_path = "/Users/jiaxchen2/Desktop/HGT_network/experiment/output_network_statistics/"

    if not os.path.isdir(bn.output_path):
        os.mkdir(bn.output_path)
    if not bn.output_path.endswith('/'):
        bn.output_path=str(bn.output_path)+'/'
    # calculate_statistics
    bn.calculate_network_node_statistics_for_list_network(bn.output_path + "network_node_statistics")

def main_calculate_network_statistics(groupInfo_file,output_dir,network_col=0,network_id_col=1):
    bn = batch_network_analysis.batch_network_analysis(groupInfo_file,output_dir,int(network_col), int(network_id_col))
    #bn.file_list_groupInfo_file = "/Users/jiaxchen2/Desktop/HGT_network/Data_by_Lichen/file_groupInfo_whole.csv"
    #bn.output_path = "/Users/jiaxchen2/Desktop/HGT_network/experiment/output_network_statistics/"
    bn.network_statistics_table = None

    if not os.path.isdir(bn.output_path):
        os.mkdir(bn.output_path)
    if not bn.output_path.endswith('/'):
        bn.output_path=str(bn.output_path)+'/'
    # calculate_statistics
    bn.calculate_network_statistics_for_list_network(bn.output_path + "network_statistics_table_large.csv")

def main_aggregate_network_by_group(groupInfo_file,output_dir,network_col=0,network_id_col=1,group_col=5):
    bn = batch_network_analysis.batch_network_analysis(groupInfo_file,output_dir,int(network_col), int(network_id_col))
    #bn.file_list_groupInfo_file = "/Users/jiaxchen2/Desktop/HGT_network/Data_by_Lichen_IBD/file_groupInfo_IBD.csv"
    #bn.output_path = "/Users/jiaxchen2/Desktop/HGT_network/experiment_IBD/output_network_species_cluster_propertySim/"
    if not os.path.isdir(bn.output_path):
        os.mkdir(bn.output_path)
    if not bn.output_path.endswith('/'):
        bn.output_path=str(bn.output_path)+'/'
    bn.aggregate_network_for_each_group_and_all("aggregated_network_", int(group_col))

def main_cluster_node_by_property(groupInfo_file,output_dir,species_network_mat_property_files, aggregated_network_file,
                                  is_adj_matrix=True,network_col=0,network_id_col=1):
    bn = batch_network_analysis.batch_network_analysis(groupInfo_file,output_dir,int(network_col), int(network_id_col))
    # run calculate network similarity and cluster
    #bn.file_list_groupInfo_file = "/Users/jiaxchen2/Desktop/HGT_network/Data_by_Lichen/file_groupInfo_whole.csv"
    #bn.output_path = "/Users/jiaxchen2/Desktop/HGT_network/experiment/output_network_species_cluster_propertySim/"

    if not os.path.isdir(bn.output_path):
        os.mkdir(bn.output_path)
    if not bn.output_path.endswith('/'):
        bn.output_path=str(bn.output_path)+'/'
    # load from file if any

    if species_network_mat_property_files is not None:
        species_network_mat_property_file_list=[]
        tmp=species_network_mat_property_files.split(',')
        for i in range(len(tmp)):
            species_network_mat_property_file_list.append(tmp[i])

    #species_network_mat_property_file_list = ['/Users/jiaxchen2/Desktop/HGT_network/experiment/output_network_statistics/network_node_statistics_degree.csv','/Users/jiaxchen2/Desktop/HGT_network/experiment/output_network_statistics/network_node_statistics_pagerank.csv','/Users/jiaxchen2/Desktop/HGT_network/experiment/output_network_statistics/network_node_statistics_clustering_coefficient.csv']
    #aggregated_network_file = "/Users/jiaxchen2/Desktop/HGT_network/experiment/output_network_species_cluster_propertySim/aggregated_network__edges.csv"

    if aggregated_network_file is not None:
        # run main_step4_network_aggregate_forvisualize
        if is_adj_matrix:
            bn.init_unique_species_pool_from_adj_file(aggregated_network_file)
            bn.aggregated_species_network = bn.load_sparse_network_withHeader(aggregated_network_file)
        else:
            bn.init_unique_species_pool_from_node_property_file(species_network_mat_property_file_list[0])
            bn.aggregated_species_network = bn.recover_adj_mat_from_edges_file(aggregated_network_file)
        bn.aggregated_species_network_simi = bn.normalize_simi_mat_to_max(bn.aggregated_species_network)
        bn.do_cluster_on_species(bn.aggregated_species_network_simi, "aggregated_HGT_network_")
    if True:
        if species_network_mat_property_file_list is not None:
            # run main_step3_node_statistics first
            # load from file and merge similarity
            bn.species_property_similarity = bn.get_species_similarity_from_multiple_property_mat(
                species_network_mat_property_file_list)
            df1 = pd.DataFrame(bn.species_property_similarity)
            df1.to_csv(bn.output_path+"species_property_similarity_mat.csv")

            bn.do_cluster_on_species(bn.species_property_similarity, "species_network_property")
            # do cluster
            #merged_similarity = bn.merge_similarity_mat_by_weight(bn.aggregated_species_network_simi, 0.5,
             #                                                       bn.species_property_similarity, 0.5)
            #bn.do_cluster_on_species(merged_similarity, "merged_")

def main_community_detection(groupInfo_file,output_dir,network_col=0,network_id_col=1):  #run use python2
    bn = batch_network_analysis.batch_network_analysis(groupInfo_file,output_dir,int(network_col), int(network_id_col))
    # community detection for each network
    # find roubust community in family and changed community
    #bn.file_list_groupInfo_file = "/Users/jiaxchen2/Desktop/HGT_network/Data_by_Lichen/file_groupInfo_whole.csv"
    #bn.output_path = "/Users/jiaxchen2/Desktop/HGT_network/experiment/output_network_community_detection/"

    if not os.path.isdir(bn.output_path):
        os.mkdir(bn.output_path)
    if not bn.output_path.endswith('/'):
        bn.output_path=str(bn.output_path)+'/'
    # find roubust community and changed community
    bn.detect_community_for_all_network(bn.output_path + "community_leidenal")

####################################
def main_community_evolution(file_list_groupInfo_file,community_color_file,output_dir,network_id_col, network_file_col,network_individual_col,network_time_col):#need modify to input according to args
    #require step: main_community_detection
    #see community_evolution_jaccard_infant.R
    #do: community info collect and calculate comm similarity and output community_link_net
    # input:"/Users/jiaxchen2/Desktop/HGT_network/Data_by_Lichen/file_groupInfo_whole.csv"
    # input: "/Users/jiaxchen2/Desktop/HGT_network/experiment/output_network_community_detection/community_leidenal_community_color.csv"
    #output: "comm_jaccard_simi_mat.csv","comm_length.csv","comm_info.csv"
    #need modify, to define the column for time, family, ...

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if not output_dir.endswith('/'):
        output_dir=str(output_dir)+'/'


    Rscript_main_do_community_detection = "Rmain_community_evolution.R"
    #community_color_file = "/Users/jiaxchen2/Desktop/HGT_network/experiment/output_network_community_detection/community_leidenal_community_color.csv"
    #file_list_groupInfo_file = "/Users/jiaxchen2/Desktop/HGT_network/Data_by_Lichen/file_groupInfo_whole.csv"
    #output_dir = "./"

    command = "Rscript " + Rscript_main_do_community_detection + " " + community_color_file + " " + \
              file_list_groupInfo_file + " " + output_dir + " " + str(network_id_col) + " " + str(network_file_col) + " "+ \
              str(network_individual_col) + " " + str(network_time_col)
    print("running script:" + command)
    os.system(command)

def main_community_cluster(comm_jaccard_simi_mat_file, output_dir, power=1, cut_height=0.6, minClusterSize =10):#need modify to input according to args
    #require step: main_community_evolution
    #see run8_infant_comm_cluster.R
    #input:"comm_jaccard_simi_mat.csv"
    #output:community cluster and partition
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if not output_dir.endswith('/'):
        output_dir = str(output_dir)+'/'

    Rscript_main_do_community_cluster = "Rmain_community_cluster.R"
    #comm_jaccard_simi_mat_file = "comm_jaccard_simi_mat.csv"
    command = "Rscript " + Rscript_main_do_community_cluster + " " + comm_jaccard_simi_mat_file+ " "+ output_dir + " "+str(power) +" "+ str(cut_height)+ " "+ str(minClusterSize)
    print("running script:" + command)
    os.system(command)

def main_hgt_event_cluster(file_list_groupInfo_file, hgt_file, output_dir, hgt_event_filter = 5, cut_height = 0.6, minClusterSize = 10):#need modify to input according to args
    #require step, lichen's hgt annotation?
    #input: "infant_cluster.hgt.bkp.txt", "file_groupInfo_whole_withid.csv"
    #set: cutoff
    #output: "hgt_hits_clubkp_filtered_5.csv",
    #output: cluster result, cluster partition...
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if not output_dir.endswith('/'):
        output_dir = str(output_dir)+'/'

    Rscript_main_do_hgt_event_cluster = "Rmain_cluster_hgt_event.R"
    #hgt_file = "infant_cluster.hgt.bkp.txt"
    #file_list_groupInfo_file = "file_groupInfo_whole_withid.csv"

    command = "Rscript " + Rscript_main_do_hgt_event_cluster + " " + \
              hgt_file + " "+ file_list_groupInfo_file + " " + output_dir +" "+ \
              str(hgt_event_filter) + " " + str(cut_height) + " " + str(minClusterSize)
    print("running script:" + command)
    os.system(command)

def main_scale_free_fit(groupInfo_file,output_dir,network_col,network_id_col):
    bn = batch_network_analysis.batch_network_analysis(groupInfo_file, output_dir, int(network_col),
                                                       int(network_id_col))
    bn.network_statistics_table = None

    if not os.path.isdir(bn.output_path):
        os.mkdir(bn.output_path)
    if not bn.output_path.endswith('/'):
        bn.output_path = str(bn.output_path) + '/'

    bn.scale_free_fit_for_list_network(bn.output_path + "network_scale_free_fit.csv")


#####################################
def main_visualization_network_similarity():#need modify to input according to args
    #require:node statistic calculate
    #note: divide group by???
    #see:1-infant_network_similarity,2_key_driver_and_network_similarity_IBD.R
    #input:network_node_statistics_pagerank,degree...
    #output:similarity plot
    pass

def main_visualization_key_nodes():
    #require:
    #input:
    #output:
    pass

def main_visualization_complexity_curve():
    #see: 3-infant_network_complexity_plot
    #note: divide group and assign time order
    #require: calculate
    #input:network_statistics_table_complete
    #output: complexity curve
    pass

def main_community_annotation():
    #require: community_evolution, community_cluster, hgt_event_cluster
    #see:
    #input:community_info, cluster_partition, acc2tex.txt, cluster.hgt.bkp.txt(hgt event and samples it appear)
    #output:
    ## corresponding visualization?
    pass


def main():
    if args.step == "network_reconstruction":
        print("running, constructing network...")
        main_network_reconstruction(abundanceFile=args.abundance_file, output_dir=args.output_dir,
                                    similarityType="correlation",distanceType="root_absolute",keep_whole_adj=args.keep_whole_adj, cutoff=0.3)

    if args.step == "node_statistics":
        print("running, calculate node statistics...")
        main_calculate_node_statistics(groupInfo_file=args.groupInfo_file,output_dir=args.output_dir,
                                       network_col=args.network_col,network_id_col=args.network_id_col)
    if args.step == "network_statistics":
        print("running, calculate network statistics...")
        main_calculate_network_statistics(groupInfo_file=args.groupInfo_file, output_dir=args.output_dir,
                                          network_col=args.network_col, network_id_col=args.network_id_col)

    if args.step == "aggregate_network":
        print("running, aggregate networks...")
        main_aggregate_network_by_group(groupInfo_file=args.groupInfo_file, output_dir=args.output_dir,
                                        network_col=args.network_col, network_id_col=args.network_id_col, group_col=args.group_col)

    if args.step == "species_cluster":
        print("running, cluster species...")
        main_cluster_node_by_property(groupInfo_file=args.groupInfo_file, output_dir=args.output_dir,
                                      species_network_mat_property_files=args.species_network_mat_property_file_list,
                                      aggregated_network_file=args.aggregated_network_file,is_adj_matrix=True,
                                      network_col=args.network_col, network_id_col=args.network_id_col)

    if args.step == "community_detection":
        print("running, community detection...")
        main_community_detection(groupInfo_file=args.groupInfo_file, output_dir=args.output_dir,
                                 network_col=args.network_col, network_id_col=args.network_id_col)

    if args.step == "community_evolution":
        print("running, community evolution...")
        main_community_evolution(file_list_groupInfo_file=args.groupInfo_file,
                                 community_color_file=args.community_color_file,output_dir=args.output_dir,
                                 network_id_col=args.network_id_col, network_file_col=args.network_col,
                                 network_individual_col=args.network_individual_col, network_time_col=args.network_time_col)

    if args.step == "community_cluster":
        print("running, community cluster...")
        main_community_cluster(comm_jaccard_simi_mat_file=args.comm_jaccard_simi_mat_file, output_dir=args.output_dir,
                               power=args.beta, cut_height=args.cut_height, minClusterSize =args.minClusterSize)

    if args.step == "hgt_event_cluster":
        print("running, HGT event cluster...")
        main_hgt_event_cluster(file_list_groupInfo_file=args.groupInfo_file, hgt_file=args.hgt_event_file, output_dir=args.output_dir,
                               hgt_event_filter = args.hgt_event_filter, cut_height=args.cut_height, minClusterSize=args.minClusterSize)

    if args.step == "scale_free_fit":
        print("running, scale free fit for all networks...")
        main_scale_free_fit(groupInfo_file=args.groupInfo_file, output_dir=args.output_dir,
                                          network_col=args.network_col, network_id_col=args.network_id_col)

    if args.step == "pipeline":
        #require groupInfo_file, output_dir, network_col, network_id_col
        #optional: network_individual_col, network_time_col
        print("running, calculate node statistics...")
        main_calculate_node_statistics(groupInfo_file=args.groupInfo_file, output_dir=args.output_dir+"out_node_statistics",
                                       network_col=args.network_col, network_id_col=args.network_id_col)
        print("running, calculate network statistics...")
        main_calculate_network_statistics(groupInfo_file=args.groupInfo_file, output_dir=args.output_dir+"out_network_statistics",
                                          network_col=args.network_col, network_id_col=args.network_id_col)
        print("running, community detection...")
        main_community_detection(groupInfo_file=args.groupInfo_file, output_dir=args.output_dir+"out_community_detection",
                                 network_col=args.network_col, network_id_col=args.network_id_col)
        print("running, community evolution...")
        main_community_evolution(file_list_groupInfo_file=args.groupInfo_file,
                                 community_color_file=args.output_dir+"out_community_detection"+"/community_leidenal_community_color.csv",
                                 output_dir=args.output_dir+"out_community_evolution",
                                 network_id_col=args.network_id_col, network_file_col=args.network_col,
                                 network_individual_col=args.network_individual_col,
                                 network_time_col=args.network_time_col)


if __name__ == '__main__':
    main()