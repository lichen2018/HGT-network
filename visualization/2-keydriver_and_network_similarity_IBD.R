library("viridis")   
library(ggplot2)


calculate_pagerank_corr=function(mat,fileHeader)
{
  network_id_list=colnames(mat)
  species_id_list = rownames(mat)
  network_value_corr_similarity=matrix(NA, length(network_id_list),length(network_id_list))
  network_jaccard_similarity_mat=matrix(NA, length(network_id_list),length(network_id_list))
  network_similarity=matrix(NA, length(network_id_list),length(network_id_list))
  for(i in 1:length(network_id_list))
  {
    pagerank_list_i=mat[,i]
    A=which(!is.na(pagerank_list_i),arr.ind = T)
    for(j in 1:length(network_id_list))
    {
      if(i>j)
      {
        pagerank_list_j=mat[,j]
        B=which(!is.na(pagerank_list_j),arr.ind = T)
        co_occurance_node=intersect(A,B)
        jaccard_similarity=length(co_occurance_node)/length(union(A,B))
        value_similarity=cor(pagerank_list_i[co_occurance_node],pagerank_list_j[co_occurance_node],method="spearman")
        network_similarity[i,j]=jaccard_similarity*value_similarity
        network_similarity[j,i]=jaccard_similarity*value_similarity
        network_jaccard_similarity_mat[i,j]=jaccard_similarity
        network_jaccard_similarity_mat[j,i]=jaccard_similarity
        network_value_corr_similarity[i,j]=value_similarity
        network_value_corr_similarity[j,i]=value_similarity
        
      }
      if(i==j)
      {
        network_similarity[i,j]=1
        network_jaccard_similarity_mat[j,i]=1
        network_value_corr_similarity[i,j]=1
      }
    }
  }
  colnames(network_similarity)=network_id_list
  rownames(network_similarity)=network_id_list
  colnames(network_jaccard_similarity_mat)=network_id_list
  rownames(network_jaccard_similarity_mat)=network_id_list
  colnames(network_value_corr_similarity)=network_id_list
  rownames(network_value_corr_similarity)=network_id_list
  write.csv(network_similarity,file=paste(fileHeader,"_network_simi.csv",sep=""),quote=F)
  write.csv(network_jaccard_similarity_mat,file=paste(fileHeader,"_network_jaccard_simi.csv",sep=""),quote=F)
  write.csv(network_value_corr_similarity,file=paste(fileHeader,"_network_property_value_simi.csv",sep=""),quote=F)
  list(network_similarity=network_similarity,network_jaccard_similarity_mat=network_jaccard_similarity_mat,network_value_corr_similarity=network_value_corr_similarity)
}


print_simi_by_comparison_cate=function(simi_mat,comparison_cate_mat,fileHeader)
{
 
  
 
 
  
 
  A=which(!is.na(simi_mat),arr.ind=T)
  print(A)
  A_left=A[,1]
  A_right=A[,2]
  B=which((A_left-A_right)>0,arr.ind = )
  select_terms=A[B,]
  
  
  values=simi_mat[select_terms]
  cate=comparison_cate_mat[select_terms]
  
  D=which(!is.na(cate),arr.ind = T)
  values=values[D]
  cate=cate[D]
  
  C=which(values==1,arr.ind = T)
  print(cate[C])
  print(select_terms[C,])
  for(i in 1:length(C))
  {
    
    network_id_list=colnames(simi_mat)
    print(network_id_list[select_terms[C[i],][1]])
    print(network_id_list[select_terms[C[i],][2]])
  }
  
  
  df=data.frame(values,cate)
  pdf(paste(fileHeader,'_by_comparison_cate_',".pdf",sep=""),width = 10,height = 12)
  p1=ggplot(df)+geom_boxplot(aes(cate,values,fill=cate))+theme(axis.text.x = element_text(angle = 45,hjust=1))+theme_classic()
  p1=p1+scale_fill_manual(values=c("mediumpurple1","lightsalmon","skyblue"))+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),legend.title  = element_blank())
  
  print(p1)
  
  dev.off()
}

filter_mat_by_size=function(degree_mat,pagerank_mat,clustering_coefficient_mat)
{
  #pagerank_mat = read.csv("network_node_statistics_pagerank.csv",header=T, row.names = 1)
  #degree_mat = read.csv("network_node_statistics_degree.csv",header = T, row.names = 1)
  #clustering_coefficient_mat = read.csv("network_node_statistics_clustering_coefficient.csv",header = T, row.names = 1)
  
  sizes=c()
  for(i in 1:dim(degree_mat)[2])
  {
    degree_neti = as.matrix(degree_mat[,i])
    A=which(!is.na(degree_neti),arr.ind = T)
    sizei=length(A)
    sizes=c(sizes,sizei)
  }
  B=which(sizes>=100,arr.ind = T)
  pagerank_mat=pagerank_mat[,B]
  degree_mat=degree_mat[,B]
  clustering_coefficient_mat=clustering_coefficient_mat[,B]
  result=list(degree_mat=degree_mat,pagerank_mat=pagerank_mat,clustering_coefficient_mat=clustering_coefficient_mat)
}


############################### end function definition###############################

###############################function definition for work 3 #########################
get_network_comparison_cate_IBD=function(simi_mat)
{
  
  #cate=c("c_nextTime_sameFamily","c_sameFamily","M_to_C_sameFamily","M_nextTime_sameFamily","M_sameFamily","other_sameFamily","different_family")
  network_id_list=colnames(simi_mat)
  #node name by individual_time_label
  states=c("non","CD","UC")
  comparision_cate=matrix(NA,length(network_id_list),length(network_id_list))
  for(i in 1:length(network_id_list))
  {
    network_id_i=network_id_list[i]
    tmpi=unlist(strsplit(as.character(network_id_i),split="_"))
    for(j in 1:length(network_id_list))
    {
      network_id_j=network_id_list[j]
      tmpj=unlist(strsplit(as.character(network_id_j),split="_"))
      statei=which(states==tmpi[3],arr.ind = T)
      statej=which(states==tmpj[3],arr.ind = T)
      if(i>j)
      {
        if(tmpi[1]==tmpj[1])
        {
          if(tmpi[3]== tmpj[3])
          {
            if(tmpi[3]=="non")
            {
              comparision_cate[i,j]=paste("sameIndividual_sameState_1non",sep="")
              comparision_cate[j,i]=paste("sameIndividual_sameState_1non",sep="")
              
            }
            if(tmpi[3]=="CD")
            {
              comparision_cate[i,j]=paste("sameIndividual_sameState_2CD",sep="")
              comparision_cate[j,i]=paste("sameIndividual_sameState_2CD",sep="")
              
            }
            if(tmpi[3]=="UC")
            {
              comparision_cate[i,j]=paste("sameIndividual_sameState_3UC",sep="")
              comparision_cate[j,i]=paste("sameIndividual_sameState_3UC",sep="")
            }
            
            
          }else{
            #comparision_cate[i,j]=paste("sameIndividual_diffState",states[min(statei,statej)],states[max(statei,statej)],sep="_")
            #comparision_cate[j,i]=paste("sameIndividual_diffState",states[min(statei,statej)],states[max(statei,statej)],sep="_")
          }
          
        }else{
          if(tmpi[3]== tmpj[3])
          {
            #comparision_cate[i,j]=paste("diffIndividual_sameState_",tmpi[3],sep="")
            #comparision_cate[j,i]=paste("diffIndividual_sameState_",tmpi[3],sep="")
            
          }else{
            #comparision_cate[i,j]=paste("diffIndividual_diffState",states[min(statei,statej)],states[max(statei,statej)],sep="_")
            #comparision_cate[j,i]=paste("diffIndividual_diffState",states[min(statei,statej)],states[max(statei,statej)],sep="_")
          }
          
        }
      }
    }
  }
  comparision_cate
}

############################### work 3 ################################################
pagerank_mat = read.csv("/Users/jiaxchen2/Desktop/HGT_network/experiment_IBD/output_network_statistics/network_node_statistics_pagerank.csv",header=T, row.names = 1)
degree_mat = read.csv("/Users/jiaxchen2/Desktop/HGT_network/experiment_IBD/output_network_statistics/network_node_statistics_degree.csv",header = T, row.names = 1)
clustering_coefficient_mat = read.csv("/Users/jiaxchen2/Desktop/HGT_network/experiment_IBD/output_network_statistics/network_node_statistics_clustering_coefficient.csv",header = T, row.names = 1)
#do filter...
#filtered=filter_mat_by_size(degree_mat,pagerank_mat,clustering_coefficient_mat)
#degree_mat=filtered$degree_mat
#pagerank_mat=filtered$pagerank_mat
#clustering_coefficient_mat=filtered$clustering_coefficient_mat

result_pagerank=calculate_pagerank_corr(pagerank_mat,"IBD_pagerank_corr")
result_degree=calculate_pagerank_corr(degree_mat,"IBD_degree_corr")
result_clustering_coefficient=calculate_pagerank_corr(clustering_coefficient_mat,"IBD_clustering_coefficient_corr")

comparison_cate = get_network_comparison_cate_IBD(result_pagerank$network_similarity)

print_simi_by_comparison_cate(result_degree$network_similarity,comparison_cate,"IBD_degree_overall")

overall_pagerank = result_pagerank$network_similarity*result_degree$network_jaccard_similarity_mat
#print_simi_by_comparison_cate(overall_pagerank,comparison_cate,"IBD_pagerank_overall")
print_simi_by_comparison_cate(result_pagerank$network_value_corr_similarity,comparison_cate,"IBD_pagerank_corr")
print_simi_by_comparison_cate(result_degree$network_value_corr_similarity,comparison_cate,"IBD_degree_corr")
print_simi_by_comparison_cate(result_degree$network_jaccard_similarity_mat,comparison_cate,"IBD_jaccard")
#print_simi_by_comparison_cate(result_degree$network_similarity,comparison_cate,"IBD_degree_overall")
print_simi_by_comparison_cate(result_clustering_coefficient$network_value_corr_similarity,comparison_cate,"IBD_clustering_coefficient_corr")
#print_simi_by_comparison_cate(result_clustering_coefficient$network_value_corr_similarity*result_degree$network_jaccard_similarity_mat,comparison_cate,"IBD_clustering_coefficient_overall")



