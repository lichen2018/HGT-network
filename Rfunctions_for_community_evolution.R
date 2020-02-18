# Title     : TODO
# Objective : TODO
# Created by: jiaxchen2
# Created on: 2020/1/4



#######define functions#####################################
get_jaccard_similarity=function(binary_vector1,binary_vector2)
{
  if(length(binary_vector1)!=length(binary_vector2))
  {
    print("the length of two vector should be the same")
  }
  A=which(binary_vector1!=0,arr.ind = T)
  B=which(binary_vector2!=0,arr.ind = T)
  C=intersect(A,B)
  D=union(A,B)
  simi=length(C)/length(D)
  list(simi=simi,num_inter=length(C),num_A=length(A),num_B=length(B),num_union=length(D))
}



get_related_network=function(network_id,ci,file_info_group,community_color_data,node_name_list, network_file_col,network_id_col)
{
  network_id_list=file_info_group[,network_id_col]
  i=which(network_id_list==network_id, arr.ind=T)

  network_id_network_list_in_community_color_data=colnames(community_color_data)
  i2=which(network_id_network_list_in_community_color_data==network_id,arr.ind=T)
  v1=community_color_data[,i2]
  A_=which(v1==ci,arr.ind = T)

  network1=read.csv(as.character(file_info_group[i,network_file_col]),header = T,row.names = 1)
  node_network1=as.array(as.matrix(rownames(network1)))
  node_in_comm=node_name_list[A_]
  indexi=c()
  for(ni in node_in_comm)
  {
    index_nodei=which(node_network1==ni,arr.ind = T)
    indexi=c(indexi,index_nodei[1])
  }

  result=network1[indexi,indexi]
  colnames(result)=node_in_comm
  rownames(result)=node_in_comm
  result
}




simi_for_communities_in_two_networks=function(i,j,simi_mat,comm_id_list,community_color_data, network_file_col)
{

  #return map with max similarity, with score
  #for i
  #for j

  v1=community_color_data[,i]
  v2=community_color_data[,j]
  network_list_in_community_color_data=colnames(community_color_data)
  unique_c_in_v1=unique(v1)
  unique_c_in_v2=unique(v2)
  #print(unique_c_in_v1)

  for(ci in unique_c_in_v1)
  {
    if(!is.na(ci))
    {

      A_=which(v1==ci,arr.ind = T)
      v_for_ci=rep(0,length(v1))
      v_for_ci[A_]=1
      id_ci=paste(network_list_in_community_color_data[i],ci,sep="_comm_")
      if(length(A_)>1)
      {

        for(cj in unique_c_in_v2)
        {
          if(!is.na(cj))
          {
            id_cj=paste(network_list_in_community_color_data[j],cj,sep="_comm_")
            v_for_cj=rep(0,length(v2))

            B_=which(v2==cj,arr.ind = T)

            v_for_cj[B_]=1
            if(length(B_>1))
            {

              simi=get_jaccard_similarity(v_for_ci,v_for_cj)
              if(simi$simi>0)
              {
                #print(paste("i",i,"j",j,"ci",ci,"cj",cj,"simi",simi$simi,"num",simi$num_inter,"num_A",simi$num_A,"num_B",simi$num_B,"num_union",simi$num_union))
                index_ci=get_index_of_comm_id(id_ci,comm_id_list)
                index_cj=get_index_of_comm_id(id_cj,comm_id_list)
                simi_mat[index_ci,index_cj]=simi$simi
                simi_mat[index_cj,index_ci]=simi$simi
              }

            }

          }
        }
      }

    }

  }

  simi_mat
}



get_index_of_comm_id=function(id,comm_id_list)
{
  A=which(comm_id_list==id,arr.ind = T)
  A[1]
}



initial_comm_list=function(network_n, community_color_data)
{
  comm_id_list=c()
  comm_length_list=c()
  accum_community_num=0
  for(i in 1:network_n)
  {
    networki=colnames(community_color_data)[i]
    #print(file_info_group[i,2])
    #print(colnames(community_color_data)[i])
    v1=community_color_data[,i]
    unique_c_in_v1=unique(v1)
    accum_community_num=accum_community_num+(length(unique(v1)))
    #print(accum_community_num)
    for(ci in unique_c_in_v1)
    {
      if(!is.na(ci))
      {
        v_for_ci=v1

        A_=which(v1==ci,arr.ind = T)
        id_ci=paste(networki,ci,sep="_comm_")
        comm_id_list=c(comm_id_list,id_ci)
        comm_length_list=c(comm_length_list,length(A_))
      }
    }

  }
  list(comm_id_list=comm_id_list,comm_length_list)
}


##############functions for decide map between two network####
get_comm_index_list_given_networki=function(comm_id_list,networki)
{
  result=c()

  for(i in 1:length(comm_id_list))
  {
    tmp=unlist(strsplit(comm_id_list[i],split="_comm_"))
    neti=tmp[1]
    ci=tmp[2]
    if(networki==neti)
    {
      result=c(result,i)
    }
  }
  result
}


collect_comm_info=function(comm_list,file_info_group,comm_info_out_file="comm_info.csv",community_color_data,node_name_list, network_file_col,network_id_col)
{

  result=c()
  comm_id_list=comm_list$comm_id_list
  #comm_id to network id to network info to network file to nodes in comm
  for(i in 1:length(comm_id_list))
  {
    comm_id=comm_id_list[i]
    tmp=unlist(strsplit(as.character(comm_id),split = "_comm_"))
    network_id=tmp[1]
    network_index_in_file_info_group=which(file_info_group[,network_id_col]==network_id,arr.ind=T)
    comm_index=tmp[2]
    adj=get_related_network(network_id,comm_index,file_info_group,community_color_data,node_name_list, network_file_col,network_id_col)
    nodes_in_comm=''
    for(j in 1:dim(adj)[1])
    {
      nodes_in_comm=paste(nodes_in_comm,rownames(adj)[j],";",sep="")
    }
    onerow=c(comm_id,network_id,comm_index)

    for(k in 1:dim(file_info_group)[2])
    {
        onerow=c(onerow, file_info_group[network_index_in_file_info_group,k])
    }
    onerow=c(onerow,nodes_in_comm)
    result=rbind(result,onerow)
  }
  result_colnames=c("comm_id","network_id","comm_label_in_network")
    for(k in 1:dim(file_info_group)[2])
    {
        result_colnames=c(result_colnames, colnames(file_info_group)[k])
    }
    result_colnames=c(result_colnames, "node_in_comm")
    colnames(result)=result_colnames
  write.csv(result,file=comm_info_out_file,quote=F,row.names = F)

}


#######################################
###############
work_on_select_network_list=function(simi_mat,comm_id_list,select_network_list,out_edges_file,jaccard_cut=0)
{
  result=c()

  for(i in 2:length(select_network_list))
  {

    comm_index_list1=get_comm_index_list_given_networki(comm_id_list,select_network_list[i-1])
    comm_index_list2=get_comm_index_list_given_networki(comm_id_list,select_network_list[i])
    select_simi_mat=simi_mat[comm_index_list1,comm_index_list2]

    rownames(select_simi_mat)=comm_id_list[comm_index_list1]
    colnames(select_simi_mat)=comm_id_list[comm_index_list2]
    #get_comm_map(select_simi_mat)
    #cutoff of jaccard, of comm size
    #for each two next time network, print edges for communities
    A=which(select_simi_mat>jaccard_cut,arr.ind = T)
    if(length(A)>0)
    {
      for(j in c(1:dim(A)[1]))
      {
        print(A[j,])
        node1=rownames(select_simi_mat)[A[j,1]]
        node2=colnames(select_simi_mat)[A[j,2]]
        weight=select_simi_mat[A[j,1],A[j,2]]
        tmp=c(node1,node2,weight)
        result=rbind(result,tmp)

      }
    }

  }

  colnames(result)=c("source","target","weight")
  write.table(result,file=out_edges_file,quote=F,row.names = F,sep=",")
}



#another function, run by each individual
run_for_each_individual=function(simi_mat,comm_id_list,file_info_group,output_dir,network_id_col,network_individual_col, network_time_col,jaccard_cut=0.1)
{
  ###modify from family id to individual id , consider infant and m
  ####need modify!!!!
  unique_individuals=unique(file_info_group[,network_individual_col])
  for(i in 1:length(unique_individuals))
  {
    individual_i=unique_individuals[i]
    A=which(file_info_group[,network_individual_col]==individual_i,arr.ind = T)
    selected_index=A
    print(paste("select_index",selected_index))
    time_point_for_select=file_info_group[selected_index, network_time_col]

    the_order=order(time_point_for_select,decreasing = F)
    ordered_selected_index=selected_index[the_order]
    if(length(ordered_selected_index)>1)
    {
      ordered_selected_network_id=file_info_group[ordered_selected_index,network_id_col]
      work_on_select_network_list(simi_mat,comm_id_list,ordered_selected_network_id,paste(output_dir,individual_i,"_links.csv",sep=""),jaccard_cut = jaccard_cut)

    }

  }
}