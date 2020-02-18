
get_one_family_by_time=function(family_id,select_m_or_c,family_ids,m_or_cs,time_ids,time_order, value_list)
{
  
  A=which(family_ids==family_id,arr.ind = T)
  B=which(m_or_cs==select_m_or_c,arr.ind = T)
  D=intersect(A,B)
  select_time=time_ids[D]
  candidate_value=value_list[D]
  select_value=c()

  for(timei in time_order)
  {
    
    E=which(select_time==timei,arr.ind = T)

    if(length(E)==0)
    {
      select_value=c(select_value, NA)
    }else{
      select_value=c(select_value, candidate_value[E])
      
    }
    
  }
  #print(select_value)
  select_value
}


plot_one_property=function(result_m, result_c,fileHeader,plot_c=T, plot_m=T)
{
  require(ggplot2)
  x_all=c()
  y_all=c()
  family_id_all=c()
  m_or_c=c()
  size=c()
  color=c()
  
  time_order=colnames(result_c)
  family_ids=rownames(result_c)
  x=c(1:length(time_order))
  
  for_eigen_c=c()
  for_eigen_m=c()
  
  if(plot_c)
  {
    for(i in c(1:dim(result_c)[1]))
    {
      tmp=which(!is.na(result_c[i,]), arr.ind = T)
      print(result_c[i,])
      print(tmp)
      print(length(tmp))
      if(length(tmp)==5)
      {
        for_eigen_c=rbind(for_eigen_c,result_c[i,])
        for(j in c(1:dim(result_c)[2]))
        {
          
          x=j
          y=result_c[i,j]
          family_id=family_ids[i]
          
          if(!is.na(y))
          {
            x_all=c(x_all,x)
            y_all=c(y_all,y)
            m_or_c=c(m_or_c,'c')
            family_id_all=c(family_id_all,family_id)
            size=c(size,0.1)
            color=c(color,'grey')
          }
          
        }
      }
    }
    v_c=get_eigen_vector(for_eigen_c,'c')
    v_c=apply(for_eigen_c,2,median)
  }
  
  
  
  if(plot_m)
  {
    for(i in c(1:dim(result_m)[1]))
    {
      tmp=which(!is.na(result_m[i,]),arr.ind=T)
      if(length(tmp)==3)
      {
        for_eigen_m=rbind(for_eigen_m,result_m[i,])
        for(j in c(1:dim(result_m)[2]))
        {
          x=j
          y=result_m[i,j]
          family_id=family_ids[i]
          
          if(!is.na(y))
          {
            x_all=c(x_all,x)
            y_all=c(y_all,y)
            m_or_c=c(m_or_c,'m')
            family_id_all=c(family_id_all,family_id)
            size=c(size,0.1)
            color=c(color,'grey')
          }
          
        }
      }
      
    }
    v_m=get_eigen_vector(for_eigen_m,'m')
    v_m=apply(for_eigen_m,2,median)
  }
  
  line_plot_x_all=c()
  line_plot_y_all=c()
  line_plot_m_or_c=c()
  line_plot_color=c()
  if(plot_c)
  {
    for(i in c(1:length(v_c)))
    {
      if(!is.na(v_c[i]))
      {
        line_plot_x_all=c(line_plot_x_all,i)
        line_plot_y_all=c(line_plot_y_all,v_c[i])
        line_plot_m_or_c=c(line_plot_m_or_c,"c_summary")
        line_plot_color=c(line_plot_color,'red')
      }
    }
  }
  
  
  if(plot_m)
  {
    for(i in c(1:length(v_m)))
    {
      if(!is.na(v_m[i]))
      {
        line_plot_x_all=c(line_plot_x_all,i)
        line_plot_y_all=c(line_plot_y_all,v_m[i])
        line_plot_m_or_c=c(line_plot_m_or_c,"m_summary")
        line_plot_color=c(line_plot_color,'red')
      }
    }
  }
  
  #x_all=(x_all)
  
  df=data.frame(x_all,y_all,family_id_all,m_or_c,size,color)
  print(df)
  
  df_line_plot=data.frame(line_plot_x_all,line_plot_y_all,line_plot_m_or_c,line_plot_color)
  
  print(df_line_plot)
  
  ###for lines version####
  #p1=ggplot(df, aes(x_all,y_all,line_type=family_id_all,color=m_or_c,size=size))+
   #           geom_line()+theme_light()+
  #  scale_colour_manual(values=c(c="grey",m="grey",c_summary="red",m_summary="blue"))
  ###for boxplot + line version####
  p1=ggplot()+geom_boxplot(data=df, aes(x=as.character(x_all), y=y_all,color=m_or_c), position = "dodge2")
  p1=p1+geom_line(data=df_line_plot,aes(x=as.numeric(as.character(line_plot_x_all)),y=line_plot_y_all,color=line_plot_m_or_c),size=2)
  p1=p1+theme_light()
  p1=p1+ylim(0,4000)
  if(plot_c)
  {
    if(plot_m)
    {
      pdf(paste(fileHeader,".pdf",sep=""))
      print(p1)
      dev.off()
    }else{
      pdf(paste(fileHeader,"_c.pdf",sep=""))
      print(p1)
      dev.off()
    }
  }else{
    if(plot_m)
    {
      pdf(paste(fileHeader,"_m.pdf",sep=""))
      print(p1)
      dev.off()
    }
  }
  
}


get_eigen_vector=function(for_eigen,m_or_c='m')
{
  if(m_or_c=='m')
  {
    select_col=c(1,2,6)
  }else{
    select_col=c(2,3,4,5,6)
  }
  
  mat=t(for_eigen[,select_col])
  pca_result <- prcomp(mat,
                       center = T,
                       scale. = T,retx=T) 
  
  print(pca_result$x)
  
  result=rep(NA,dim(for_eigen)[2])
  result[select_col]=pca_result$x[,1] 
  result
}


all_family=function(network_statistic, time_order_C, time_order_M,one_property)
{
  network_ids=network_statistic$network_id
  
  family_ids=c()
  time_ids=c()
  m_or_cs=c()
  for(i in c(1:length(network_ids)))
  {
    tmp=unlist(strsplit(as.character(network_ids[i]), split="_"))
    family_ids=c(family_ids,tmp[1])
    time_ids=c(time_ids,tmp[2])
    m_or_cs=c(m_or_cs,tmp[3])
  }

  value_index=A=which(colnames(network_statistic)==one_property, arr.ind = T)
  value_list=network_statistic[,value_index]
  print(value_list)
  unique_family=unique(family_ids)
  
  result_m=c()
  result_c=c()
  for(family_id in unique_family)
  {
    print(family_id)
    
    select_value_list_c=get_one_family_by_time(family_id,"C",family_ids,m_or_cs,time_ids,time_order_C, value_list)
    select_value_list_m=get_one_family_by_time(family_id,"M",family_ids,m_or_cs,time_ids,time_order_M, value_list)
    
    result_m=rbind(result_m, select_value_list_m)
    result_c=rbind(result_c, select_value_list_c)
  }
  rownames(result_m)=unique_family
  rownames(result_c)=unique_family
  colnames(result_c)=time_order_C
  colnames(result_m)=time_order_M
  
  plot_one_property(result_m,result_c,one_property,T,F)
  plot_one_property(result_m,result_c,one_property,F,T)
  plot_one_property(result_m,result_c,one_property,T,T)
  
}




#################end function definition###########################

network_statistic=read.csv("network_statistics_table_complete.csv") 

time_order_C=c("Gest","Birth","14","1","2","3")
time_order_M=c("Gest","Birth","14","1","2","3")


all_family(network_statistic, time_order_C, time_order_M,"size")

all_family(network_statistic, time_order_C, time_order_M,"HVN")
