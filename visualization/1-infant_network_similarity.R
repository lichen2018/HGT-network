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

get_time_diff=function(time1,time2)
{
  time_diff_cate=NA
  
  if(time1 == "Gest")
  {
    if(time2 == "Birth")
    {
      time_diff_cate="1_Gest_Birth"
    }
  }
  if(time1 == "Birth")
  {
    if(time2 == "Gest")
    {
      time_diff_cate="1_Gest_Birth"
    }
    
    if(time2 == "14")
    {
      time_diff_cate="2_Birth_14days"
    }
    
    if(time2 == "3")
    {
      time_diff_cate="2_Birth_3month"
    }
    
    
  }
  if(time1 == "14")
  {
    if(time2 == "Birth")
    {
      time_diff_cate="2_Birth_14days"
    }
    if(time2 == "1")
    {
      time_diff_cate="3_14days_1month"
    }
  }
  if(time1 == "1")
  {
   
    if(time2 == "14")
    {
      time_diff_cate="3_14days_1month"
    }
    if(time2 == "2")
    {
      time_diff_cate="4_1month_2month"
    }
  }
  if(time1 == "2")
  {
   
    if(time2 == "1")
    {
      time_diff_cate="4_1month_2month"
    }
    if(time2 == "3")
    {
      time_diff_cate="5_2month_3month"
    }
  }
  if(time1 == "3")
  {
    if(time2 == "Birth")
    {
      time_diff_cate="2_Birth_3month"
    }
    if(time2 == "2")
    {
      time_diff_cate="5_2month_3month"
    }
  }
  

  time_diff_cate
}


get_time_diff_MC=function(time1,time2)
{
  time_diff_cate=NA
  
  if(time1 == "Gest")
  {
    #if(time2 == "Birth")
    #{
    #  time_diff_cate="1_Gest_Birth"
    #}
    if(time2 == "3")
    {
      time_diff_cate="5_Gest_3month"
    }
    
  }
  if(time1 == "Birth")
  {
    #if(time2 == "Gest")
    #{
    #  time_diff_cate="1_Gest_Birth"
    #}
    if(time2 == "Birth")
    {
      time_diff_cate="0_Birth_Birth"
    }
    
    if(time2 == "14")
    {
      time_diff_cate="2_Birth_14days"
    }
    if(time2 == "1")
    {
      time_diff_cate="2_Birth_1month"
    }
    if(time2 == "2")
    {
      time_diff_cate="2_Birth_2month"
    }
    if(time2 == "3")
    {
      time_diff_cate="2_Birth_3month"
    }
    
    
  }
  if(time1 == "14")
  {
    if(time2 == "Birth")
    {
      time_diff_cate="2_Birth_14days"
    }
    if(time2 == "1")
    {
      time_diff_cate="3_14days_1month"
    }
  }
  if(time1 == "1")
  {
    if(time2 == "Birth")
    {
      time_diff_cate="2_Birth_1month"
    }
    if(time2 == "14")
    {
      time_diff_cate="3_14days_1month"
    }
    if(time2 == "2")
    {
      time_diff_cate="4_1month_2month"
    }
  }
  if(time1 == "2")
  {
    if(time2 == "Birth")
    {
      time_diff_cate="2_Birth_2month"
    }
    if(time2 == "1")
    {
      time_diff_cate="4_1month_2month"
    }
    #if(time2 == "3")
    #{
    #  time_diff_cate="5_2month_3month"
    #}
  }
  if(time1 == "3")
  {
    if(time2 == "Birth")
    {
      time_diff_cate="2_Birth_3month"
    }
    #if(time2 == "2")
    #{
    #  time_diff_cate="5_2month_3month"
    #}
    if(time2 == "3")
    {
      time_diff_cate="5_3month_3month"
    }
    if(time2 == "Gest")
    {
      time_diff_cate="5_Gest_3month"
    }
    
  }
  
  time_diff_cate
}



get_network_comparison_cate=function(simi_mat)
{
  
  #cate=c("c_nextTime_sameFamily","c_sameFamily","M_to_C_sameFamily","M_nextTime_sameFamily","M_sameFamily","other_sameFamily","different_family")
  network_id_list=colnames(simi_mat)
  comparision_cate=matrix(NA,length(network_id_list),length(network_id_list))
  
  comparision_cate_CC=matrix(NA,length(network_id_list),length(network_id_list))
  comparision_cate_MM=matrix(NA,length(network_id_list),length(network_id_list))
  comparision_cate_MC=matrix(NA,length(network_id_list),length(network_id_list))
  comparision_color=matrix(NA,length(network_id_list),length(network_id_list))
  for(i in 1:length(network_id_list))
  {
    network_id_i=network_id_list[i]
    tmpi=unlist(strsplit(as.character(network_id_i),split="_"))
    for(j in 1:length(network_id_list))
    {
      network_id_j=network_id_list[j]
      tmpj=unlist(strsplit(as.character(network_id_j),split="_"))
      if(i>j)
      {
        #tmp[1] family tmp[2] time tmp[3] C or M
        if(tmpi[1]==tmpj[1])
        {
          if(tmpi[3]=='C' && tmpj[3]=='C') #same family, CC
          {
            comparision_cate[i,j]="9 sameFamily_C_different_time"
            comparision_cate[j,i]="9 sameFamily_C_different_time"
            
            
            time_diff_cate=get_time_diff(tmpi[2],tmpj[2])
            if(!is.na(time_diff_cate))
            {
              if(time_diff_cate=="2_Birth_3month")
              {
                comparision_cate_CC[i,j]=NA
                comparision_cate_CC[j,i]=NA
              }else{
                comparision_cate_CC[i,j]=time_diff_cate
                comparision_cate_CC[j,i]=time_diff_cate
              }
            }
            
            comparision_color[i,j]="same family"
            comparision_color[j,i]="same family"
            
            
          }else{
            if(tmpi[3]=='M' && tmpj[3]=='M') #same family MM
            {
              comparision_cate[i,j]="8 sameFamily_M_different_time"
              comparision_cate[j,i]="8 sameFamily_M_different_time"
              
              time_diff_cate=get_time_diff(tmpi[2],tmpj[2])
              comparision_cate_MM[i,j]=time_diff_cate
              comparision_cate_MM[j,i]=time_diff_cate
              
              comparision_color[i,j]="same family"
              comparision_color[j,i]="same family"
              
            }else{ #same family C M 
              if(tmpi[2]==tmpj[2]) 
              {
                comparision_cate[i,j]="7 sameFamily_M_to_C_same_time"
                comparision_cate[j,i]="7 sameFamily_M_to_C_same_time"
                
                time_diff_cate=get_time_diff_MC(tmpi[2],tmpj[2])
                comparision_cate_MC[i,j]=time_diff_cate
                comparision_cate_MC[j,i]=time_diff_cate
                
                comparision_color[i,j]="same family"
                comparision_color[j,i]="same family"
              }else{
                comparision_cate[i,j]="6 sameFamily_M_to_C_different_time"
                comparision_cate[j,i]="6 sameFamily_M_to_C_different_time"
                
                time_diff_cate=get_time_diff_MC(tmpi[2],tmpj[2])
                if(!is.na(time_diff_cate))
                {
                  
                  if(time_diff_cate=="2_Birth_3month")
                  {
                    print(time_diff_cate)
                    print(paste(tmpi[2],tmpi[3],tmpj[2],tmpj[3]))
                    if(tmpi[3]=='M')
                    {
                      if(tmpi[2]=='Birth')
                      {
                        #time_diff_cate="2_MBirth_C3month"
                      }else{
                        time_diff_cate=NA
                      }
                    }else{ 
                      if(tmpi[2]=='Birth')
                      {
                        time_diff_cate=NA
                      }else{
                        #time_diff_cate="2_MBirth_C3month"
                      }
                    }
                  }
                }
                comparision_cate_MC[i,j]=time_diff_cate
                comparision_cate_MC[j,i]=time_diff_cate
                
                comparision_color[i,j]="same family"
                comparision_color[j,i]="same family"
              }
              
            }
          }
          
        }else{
          if(tmpi[2]==tmpj[2]) 
          {
            if(tmpi[3]=='C' && tmpj[3]=='C') #diff family CC
            {
              comparision_cate[i,j]="1 differentFamily_C_same_time"
              comparision_cate[j,i]="1 differentFamily_C_same_time"
              
              comparision_color[i,j]="diff family"
              comparision_color[j,i]="diff family"
              
            }else{
              if(tmpi[3]=='M' && tmpj[3]=='M') #diff family MM
              {
                comparision_cate[i,j]="3 differentFamily_M_same_time"
                comparision_cate[j,i]="3 differentFamily_M_same_time"
                
                comparision_color[i,j]="diff family"
                comparision_color[j,i]="diff family"
              }else{ #diff family MC 
                comparision_cate[i,j]="5 differentFamily_M_to_C_same_time"
                comparision_cate[j,i]="5 differentFamily_M_to_C_same_time"
                
                time_diff_cate=get_time_diff_MC(tmpi[2],tmpj[2])
                if(!is.na(time_diff_cate))
                {
                  
                  if(time_diff_cate=="2_Birth_3month")
                  {
                    print(time_diff_cate)
                    print(paste(tmpi[2],tmpi[3],tmpj[2],tmpj[3]))
                    if(tmpi[3]=='M')
                    {
                      if(tmpi[2]=='Birth')
                      {
                        #time_diff_cate="2_MBirth_C3month"
                      }else{
                        time_diff_cate=NA#"2_M3month_CBirth"
                      }
                    }else{ 
                      if(tmpi[2]=='Birth')
                      {
                        time_diff_cate=NA
                      }else{
                        #time_diff_cate="2_MBirth_C3month"
                      }
                    }
                  }
                }
                comparision_cate_MC[i,j]=time_diff_cate
                comparision_cate_MC[j,i]=time_diff_cate
                
                comparision_color[i,j]="diff family"
                comparision_color[j,i]="diff family"
              }
              
            }
            
          }else{
            if(tmpi[3]=='C' && tmpj[3]=='C')
            {
              comparision_cate[i,j]="0 differentFamily_C_different_time"
              comparision_cate[j,i]="0 differentFamily_C_different_time"
              
              time_diff_cate=get_time_diff(tmpi[2],tmpj[2])
              if(!is.na(time_diff_cate))
              {
                if(time_diff_cate=="2_Birth_3month")
                {
                  comparision_cate_CC[i,j]=NA
                  comparision_cate_CC[j,i]=NA
                }else{
                  comparision_cate_CC[i,j]=time_diff_cate
                  comparision_cate_CC[j,i]=time_diff_cate
                }
              }
              
              comparision_color[i,j]="diff family"
              comparision_color[j,i]="diff family"
              
            }else{
              if(tmpi[3]=='M' && tmpj[3]=='M')
              {
                comparision_cate[i,j]="2 differentFamily_M_different_time"
                comparision_cate[j,i]="2 differentFamily_M_different_time"
                
                time_diff_cate=get_time_diff(tmpi[2],tmpj[2])
                comparision_cate_MM[i,j]=time_diff_cate
                comparision_cate_MM[j,i]=time_diff_cate
                
                comparision_color[i,j]="diff family"
                comparision_color[j,i]="diff family"
              }else{
                comparision_cate[i,j]="4 differentFamily_M_to_C_different_time"
                comparision_cate[j,i]="4 differentFamily_M_to_C_different_time"
                
                time_diff_cate=get_time_diff_MC(tmpi[2],tmpj[2])
                if(!is.na(time_diff_cate))
                {
                  
                  if(time_diff_cate=="2_Birth_3month")
                  {
                    print(time_diff_cate)
                    print(paste(tmpi[2],tmpi[3],tmpj[2],tmpj[3]))
                    if(tmpi[3]=='M')
                    {
                      if(tmpi[2]=='Birth')
                      {
                        #time_diff_cate="2_MBirth_C3month"
                      }else{
                        time_diff_cate=NA
                      }
                    }else{ 
                      if(tmpi[2]=='Birth')
                      {
                        time_diff_cate=NA
                      }else{
                        #time_diff_cate="2_MBirth_C3month"
                      }
                    }
                  }
                }
                comparision_cate_MC[i,j]=time_diff_cate
                comparision_cate_MC[j,i]=time_diff_cate
                
                comparision_color[i,j]="diff family"
                comparision_color[j,i]="diff family"
              }
              
            }
            
          }
          
        }
      }
    }
  }
  result=list(comparison_cate=comparision_cate,comparison_cate_CC=comparision_cate_CC,
       comparison_cate_MM=comparision_cate_MM,comparison_cate_MC=comparision_cate_MC,
       comparison_color=comparision_color)
}

print_simi_by_comparison_cate=function(simi_mat,comparison_cate_mat,comparison_color_mat,fileHeader)
{
 
  A=which(!is.na(simi_mat),arr.ind=T)
  A_left=A[,1]
  A_right=A[,2]
  B=which((A_left-A_right)>0,arr.ind = )
  select_terms=A[B,]
  
  values=simi_mat[select_terms]
  cate=comparison_cate_mat[select_terms]
  color_v=comparison_color_mat[select_terms]
  
  D=which(!is.na(cate),arr.ind = T)
  values=values[D]
  cate=cate[D]
  color_v=color_v[D]
  
  C=which(values==1,arr.ind = T)
  for(i in 1:length(C))
  {
    
    network_id_list=colnames(simi_mat)
    print(network_id_list[select_terms[C[i],][1]])
    print(network_id_list[select_terms[C[i],][2]])
  }
  
  
  df=data.frame(values,cate,color_v)
  
  tmp_pvalues=c()
  unique_cate=unique(cate)
  for(catei in unique_cate)
  {
    index_in_catei=which(cate==catei,arr.ind = T)
    values_in_catei=values[index_in_catei]
    color_in_catei=color_v[index_in_catei]
    unique_color=unique(color_in_catei)
    if(length(unique_color)==2)
    {
      A1=which(color_in_catei==unique_color[1],arr.ind = T)
      A2=which(color_in_catei==unique_color[2],arr.ind = T)
      
      x=values_in_catei[A1]
      y=values_in_catei[A2]
      print(catei)
      tmp_pvalues=rbind(tmp_pvalues, catei)
      r1=wilcox.test(x,y)
      print(r1$p.value)
      tmp_pvalues=rbind(tmp_pvalues, r1$p.value)
      r2=t.test(x,y)
      print(r2$p.value)
      tmp_pvalues=rbind(tmp_pvalues, r2$p.value)
    }
  }
  
  write.csv(tmp_pvalues,file=paste(fileHeader,'_pvalues.csv',sep=""))
  require(RColorBrewer)
  pdf(paste(fileHeader,'_by_comparison_cate_',".pdf",sep=""),width = 10,height = 12)
  p1=ggplot(df)+geom_boxplot(aes(cate,values,fill=color_v))+theme(axis.text.x = element_text(angle = 45,hjust=1))
  p1=p1+coord_cartesian(ylim=c(-0.1,0.4))+theme_classic()
  #p1=p1+coord_cartesian()+theme_classic()
  p1=p1+scale_fill_manual(values=c("lightsalmon","skyblue"))+theme(axis.title.x=element_blank(),legend.title  = element_blank())
  
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
  #B=which(sizes>=100,arr.ind = T)
  B=which(sizes>=0,arr.ind = T)
  pagerank_mat=pagerank_mat[,B]
  degree_mat=degree_mat[,B]
  clustering_coefficient_mat=clustering_coefficient_mat[,B]
  result=list(degree_mat=degree_mat,pagerank_mat=pagerank_mat,clustering_coefficient_mat=clustering_coefficient_mat)
}


############################### end function definition###############################
pagerank_mat = read.csv("network_node_statistics_pagerank.csv",header=T, row.names = 1)
degree_mat = read.csv("network_node_statistics_degree.csv",header = T, row.names = 1)
clustering_coefficient_mat = read.csv("network_node_statistics_clustering_coefficient.csv",header = T, row.names = 1)


result_pagerank=calculate_pagerank_corr(pagerank_mat,"pagerank_corr_filter")
result_degree=calculate_pagerank_corr(degree_mat,"degree_corr_filter")
result_clustering_coefficient=calculate_pagerank_corr(clustering_coefficient_mat,"clustering_coefficient_corr_filter")

overall_pagerank_degree = result_pagerank$network_similarity*result_degree$network_jaccard_similarity_mat

comparison_cate = get_network_comparison_cate(result_pagerank$network_similarity)


#overall_pagerank_degree
#result_pagerank$network_value_corr_similarity
#result_clustering_coefficient$network_value_corr_similarity
#result_degree$network_value_corr_similarity
#result_degree$network_similarity
#result_degree$network_jaccard_similarity_mat

similarity_mat=result_degree$network_similarity
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_CC,comparison_cate$comparison_color,"CC_degree")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MM,comparison_cate$comparison_color,"MM_degree")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MC,comparison_cate$comparison_color,"MC_degree")

tmp=function()
{
similarity_mat=result_degree$network_value_corr_similarity
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_CC,comparison_cate$comparison_color,"CC_degree_value")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MM,comparison_cate$comparison_color,"MM_degree_value")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MC,comparison_cate$comparison_color,"MC_degree_value")



similarity_mat=result_degree$network_jaccard_similarity_mat
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_CC,comparison_cate$comparison_color,"CC_jaccard")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MM,comparison_cate$comparison_color,"MM_jaccard")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MC,comparison_cate$comparison_color,"MC_jaccard")

similarity_mat=overall_pagerank_degree
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_CC,comparison_cate$comparison_color,"CC_pagerank_degree")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MM,comparison_cate$comparison_color,"MM_pagerank_degree")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MC,comparison_cate$comparison_color,"MC_pagerank_degree")


similarity_mat=result_pagerank$network_value_corr_similarity
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_CC,comparison_cate$comparison_color,"CC_pagerank_value")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MM,comparison_cate$comparison_color,"MM_pagerank_value")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MC,comparison_cate$comparison_color,"MC_pagerank_value")


similarity_mat=result_clustering_coefficient$network_value_corr_similarity
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_CC,comparison_cate$comparison_color,"CC_coef_value")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MM,comparison_cate$comparison_color,"MM_coef_value")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MC,comparison_cate$comparison_color,"MC_coef_value")
##########################################
#do filter...
filtered=filter_mat_by_size(degree_mat,pagerank_mat,clustering_coefficient_mat)
degree_mat=filtered$degree_mat
pagerank_mat=filtered$pagerank_mat
clustering_coefficient_mat=filtered$clustering_coefficient_mat

result_pagerank=calculate_pagerank_corr(pagerank_mat,"pagerank_corr_filter")
result_degree=calculate_pagerank_corr(degree_mat,"degree_corr_filter")
result_clustering_coefficient=calculate_pagerank_corr(clustering_coefficient_mat,"clustering_coefficient_corr_filter")

overall_pagerank_degree = result_pagerank$network_similarity*result_degree$network_jaccard_similarity_mat

comparison_cate = get_network_comparison_cate(result_pagerank$network_similarity)
##########################
similarity_mat=result_degree$network_value_corr_similarity
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_CC,comparison_cate$comparison_color,"CC_degree_value_filter")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MM,comparison_cate$comparison_color,"MM_degree_value_filter")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MC,comparison_cate$comparison_color,"MC_degree_value_filter")

similarity_mat=result_degree$network_similarity
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_CC,comparison_cate$comparison_color,"CC_degree_filter")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MM,comparison_cate$comparison_color,"MM_degree_filter")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MC,comparison_cate$comparison_color,"MC_degree_filter")

similarity_mat=result_degree$network_jaccard_similarity_mat
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_CC,comparison_cate$comparison_color,"CC_jaccard_filter")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MM,comparison_cate$comparison_color,"MM_jaccard_filter")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MC,comparison_cate$comparison_color,"MC_jaccard_filter")

similarity_mat=overall_pagerank_degree
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_CC,comparison_cate$comparison_color,"CC_pagerank_degree_filter")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MM,comparison_cate$comparison_color,"MM_pagerank_degree_filter")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MC,comparison_cate$comparison_color,"MC_pagerank_degree_filter")


similarity_mat=result_pagerank$network_value_corr_similarity
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_CC,comparison_cate$comparison_color,"CC_pagerank_value_filter")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MM,comparison_cate$comparison_color,"MM_pagerank_value_filter")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MC,comparison_cate$comparison_color,"MC_pagerank_value_filter")


similarity_mat=result_clustering_coefficient$network_value_corr_similarity
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_CC,comparison_cate$comparison_color,"CC_coef_value_filter")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MM,comparison_cate$comparison_color,"MM_coef_value_filter")
print_simi_by_comparison_cate(similarity_mat,comparison_cate$comparison_cate_MC,comparison_cate$comparison_color,"MC_coef_value_filter")
}


#overall_pagerank = result_pagerank$network_similarity*result_degree$network_jaccard_similarity_mat
#print_simi_by_comparison_cate(overall_pagerank,comparison_cate,"pagerank_overall_filter")
#print_simi_by_comparison_cate(result_pagerank$network_value_corr_similarity,comparison_cate,"pagerank_corr_filter")
#print_simi_by_comparison_cate(result_degree$network_value_corr_similarity,comparison_cate,"degree_corr_filter")
#print_simi_by_comparison_cate(result_degree$network_jaccard_similarity_mat,comparison_cate,"jaccard_filter")
#print_simi_by_comparison_cate(result_degree$network_similarity,comparison_cate,"degree_overall_filter")
#print_simi_by_comparison_cate(result_clustering_coefficient$network_value_corr_similarity,comparison_cate,"clustering_coefficient_corr_filter")


