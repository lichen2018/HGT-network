# Title     : TODO
# Objective : TODO
# Created by: jiaxchen2
# Created on: 2020/1/22
source("Rfunction_for_cluster.R")

###########work#################

args = commandArgs(TRUE)
simi_mat_file = args[1]
output_dir = args[2]
power = as.integer(args[3])
cutHeight= as.double(args[4])
minClusterSize=as.integer(args[5])


simi_mat=read.csv(simi_mat_file,header = T,row.names = 1)

A=which(is.na(simi_mat),arr.ind=T)
simi_mat[A]=0
dist_mat=1-simi_mat
#add size distance by weight??
#d=as.dist(dist_mat)
dist_mat=dist_mat^power


do_cluster_from_distance(dist_mat,rownames(dist_mat),paste(output_dir,"comm_cluster_height_",cutHeight,"_minClusterSize_",minClusterSize,sep=""),"dynamic_tree_cut","complete",cutHeight=cutHeight, minClusterSize=minClusterSize)



