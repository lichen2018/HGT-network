# Title     : TODO
# Objective : TODO
# Created by: jiaxchen2
# Created on: 2020/1/22


plot_cluster_heatmap=function(distanceMatrices,order,labels,out_plot_file)
{
  dist=distanceMatrices[order,order]
  pdf(out_plot_file)
  ggplot_heatmap_continues(dist,rep(' ',length(labels)),rep(' ',length(labels)),"distance","firebrick4","white",T,"blue")
  dev.off()
}


ggplot_heatmap_continues=function(mat,xlabel,ylabel,fill_legend,color_low="white",color_high="steelblue",is_coord_equal=F,mid_color=NULL)
{
  require("ggplot2")
  require("reshape2")
  require("gplots")
  require(ggrepel)
  require(scales)

  colnames(mat)=c(1:dim(mat)[2])
  rownames(mat)=c(1:dim(mat)[1])

  mat <- data.frame(mat)
  mat$id<-rownames(mat)
  melt_mat <- melt(mat, id.var="id")

  if(is.null(mid_color))
  {
    p1=ggplot(melt_mat,aes(as.integer(variable), as.integer(id)))+
      geom_tile(aes(fill=melt_mat$value),colour="white")+
      scale_fill_gradient(fill_legend,low=color_low,high=color_high)
  }else{
    p1=ggplot(melt_mat,aes(as.integer(variable), as.integer(id)))+
      geom_tile(aes(fill=melt_mat$value),colour="white")+
      scale_fill_gradient2(fill_legend,low=color_low,mid=mid_color,high=color_high)
  }


  p1=p1 + scale_x_continuous(breaks=seq(1:length(xlabel)),labels=as.array(xlabel))
  p1=p1 + scale_y_continuous(breaks=seq(1:length(ylabel)),labels=as.array(ylabel))
  p1=p1+theme(axis.text.x = element_text(angle=75, hjust=1)
              ,axis.title.x=element_blank(),axis.title.y = element_blank())
  if(is_coord_equal)
  {
    p1=p1+coord_equal()
  }
  print(p1)
}

hcluster_and_dynamicTreeCut=function(distanceMatrices, out_partition_file, out_dendrogram_pdf, hcluster_mode="ave",deepSplit=F, cutHeight=0.6, minClusterSize=20)
{
  library(WGCNA)
  hc=hclust(as.dist(distanceMatrices),hcluster_mode)
  ################################################################
    write.csv(hc$order,paste(out_dendrogram_pdf,"hc_order.csv",sep=""),quote=F)
  write.csv(hc$labels,paste(out_dendrogram_pdf,"hc_labels.csv",sep=""),quote=F)
  write.csv(hc$merge,paste(out_dendrogram_pdf,"hc_merge.csv",sep=""),quote=F)

  colorh1=cutreeDynamic(hc, deepSplit=F, cutHeight=cutHeight, minClusterSize=minClusterSize)
  pdf(out_dendrogram_pdf)
  par(mfrow=c(2,1))
  plot(hc, main="cluster tree", labels=F, xlab="", sub="");
  abline(h=cutHeight,col="deepskyblue")
  plotColorUnderTree(hc,colors=data.frame(colorh1),rowText=colorh1)
  #plot_cluster_heatmap(distanceMatrices,hc$order,hc$labels,paste(out_dendrogram_pdf,"_dist_heatmap.pdf"))
  title("Colored by UNMERGED dynamic modules")
  dev.off()
  write.csv(cbind(hc$labels,colorh1),out_partition_file,quote = F,row.names = F)
  colorh1
}


do_cluster_from_distance=function(distanceMatrices1,node_name,output_cluster_result_fileHeader,cluster_method,cluster_mode,cluster_num=20,cutHeight=0.6,minClusterSize=20)
{
  dendrogram_file=paste(output_cluster_result_fileHeader,"_dendrogram_",cluster_method,".pdf",sep="")
  partition_file= paste(output_cluster_result_fileHeader,"_clusterPartition_",cluster_method,".csv",sep="")
  if(cluster_method == "dynamic_tree_cut")
  {
    hcluster_and_dynamicTreeCut(distanceMatrices1, partition_file,dendrogram_file, hcluster_mode=cluster_mode,deepSplit=F, cutHeight=cutHeight, minClusterSize=minClusterSize)
  }else{
    memb1=doCluster(distanceMatrices1,cluster_method,cluster_num,dendrogram_file)
    write.table(cbind(as.character(node_name),memb1),file =partition_file,quote=F,sep=",",row.names = F,col.names = F)

  }
}
