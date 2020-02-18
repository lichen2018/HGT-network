# Title     : TODO
# Objective : TODO
# Created by: jiaxchen2
# Created on: 2020/1/22

##################
source("Rfunction_for_cluster.R")


args = commandArgs(TRUE)
hgt_file = args[1]
file_list_groupInfo_file = args[2]
output_dir = args[3]
hgt_event_filter=as.integer(args[4])
cutHeight= as.double(args[5])
minClusterSize=as.integer(args[6])

network_info=read.csv(file_list_groupInfo_file,header = T)
hgt_data=read.table(hgt_file, header=F,sep="\n")
network_list=network_info[,2]
hgt_hits=c()
hgt_id_list=c()

for(i in 1:dim(hgt_data)[1])
{
  tmp=unlist(strsplit(as.character(hgt_data[i,1]),split=","))

  if(length(tmp)>(4+hgt_event_filter)) #need to redefine hgt event file format
  {
    print(i)
    hgt_id=paste(tmp[1],tmp[2],tmp[3],tmp[4],sep=";")
    hits=tmp[5:length(tmp)]
    one_hgt_hits=rep(0,length(network_list))
    for(j in 1:length(hits))
    {
      A=which(network_list==hits[j],arr.ind = T)
      one_hgt_hits[A]=1
    }
    hgt_hits=rbind(hgt_hits,one_hgt_hits)
    hgt_id_list=c(hgt_id_list,hgt_id)
  }

}
#copy and work from here
rownames(hgt_hits)=hgt_id_list
colnames(hgt_hits)=network_list
B=which(hgt_hits==0,arr.ind=T)
store_hgt_hits=hgt_hits
store_hgt_hits[B]=NA
write.csv(store_hgt_hits,file=paste(output_dir,"hgt_hits_clubkp_filtered.csv",sep=""),quote=F)

d=dist(hgt_hits,method="binary")

do_cluster_from_distance(as.matrix(d),hgt_id_list,paste(output_dir,"hgt_event_cluster_height_",cutHeight,"_minClusterSize_",minClusterSize,sep=""),"dynamic_tree_cut","complete",cutHeight=cutHeight,minClusterSize=minClusterSize)
