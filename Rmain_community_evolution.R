# Title     : TODO
# Objective : TODO
# Created by: jiaxchen2
# Created on: 2020/1/4



source("Rfunctions_for_community_evolution.R")



work_community_evolution=function(community_color_file, file_list_groupInfo_file, output_dir, network_id_col, network_file_col, network_individual_col=NA, network_time_col=NA, jaccard_cut=0.1)
{

    community_color_data=read.csv(community_color_file,header = T, row.names = 1) #gloabal variable
    #change the file community_color_data
    #node_property_mat_file ="/Users/jiaxchen2/Desktop/HGT_network/experiment_IBD/output_network_statistics/network_node_statistics_degree.csv"
    #aggregated_network_file = "/Users/jiaxchen2/Desktop/HGT_network/experiment_IBD/output_network_species_cluster_propertySim/aggregated_network__edges.csv"
    #mat1 = read.csv(node_property_mat_file, header = T,row.names = 1)
    node_name_list = rownames(community_color_data)

    file_info_group=read.csv(file_list_groupInfo_file,header = T)

    ###run jaccard######
    #1.jaccard, 2.similarity of community statistics.
    #visualize, find common domain, for one individual
    #cluster similar, across individual and label


    network_n=dim(community_color_data)[2]
    species_n=dim(community_color_data)[1]
    comm_list=initial_comm_list(network_n, community_color_data)
    comm_id_list=comm_list$comm_id_list
    simi_mat=matrix(NA,length(comm_id_list),length(comm_id_list))
    for(i in 1:network_n)
    {
      for(j in 1:network_n)
      {
        if(j > i)
        {
          print(paste("i",i,"j",j))
          simi_mat=simi_for_communities_in_two_networks(i,j,simi_mat,comm_id_list,community_color_data, network_file_col)
        }
      }
    }
    colnames(simi_mat)=comm_id_list
    rownames(simi_mat)=comm_id_list
    write.csv(simi_mat,file=paste(output_dir,"comm_jaccard_simi_mat.csv",sep=""),quote=F)
    write.table(comm_list,file=paste(output_dir,"comm_length.csv", sep=""),quote=F,col.names = F,row.names = F,sep=",")
    #######
    collect_comm_info(comm_list,file_info_group,comm_info_out_file=paste(output_dir,"comm_info.csv", sep=""),community_color_data,node_name_list, network_file_col,network_id_col)

    if(!is.na(network_individual_col))
    {
        if(!is.na(network_time_col))
        {
            run_for_each_individual(simi_mat,comm_id_list,file_info_group,output_dir,network_id_col, network_individual_col, network_time_col,jaccard_cut=jaccard_cut)
        }
    }

}

########################end define functions#########################

args = commandArgs(TRUE)
community_color_file = args[1]
file_list_groupInfo_file = args[2]
output_dir = args[3]

network_id_col= as.integer(args[4])+1
network_file_col = as.integer(args[5])+1
network_individual_col = as.character(args[6])
network_time_col = as.character(args[7])

if(network_individual_col == 'NA')
{
    network_individual_col=NA
}else{
    network_individual_col=as.integer(network_individual_col)+1
}

if(network_time_col == 'NA')
{
    network_time_col=NA
}else{
    network_time_col=as.integer(network_time_col)+1
}


work_community_evolution(community_color_file,file_list_groupInfo_file,output_dir, network_id_col,network_file_col, network_individual_col, network_time_col)


