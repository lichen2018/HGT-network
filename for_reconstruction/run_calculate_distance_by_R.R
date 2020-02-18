# Title     : TODO
# Objective : TODO
# Created by: jiaxchen2
# Created on: 19/3/6

require("lsa")


produceDistanceMatrices=function(expr,method)
{
  texpr=t(expr)
  distanceMatrices=NULL
  if(method=="absolute")#absolute Pearson correlation distance
  {
    corrMatrix=cor(texpr)
    distanceMatrices=1-abs(corrMatrix)
  }
  if(method=="root_absolute")#square root absolute correlation distance
  {
    corrMatrix=cor(texpr)
    distanceMatrices=1-abs(corrMatrix)
    distanceMatrices=sqrt(distanceMatrices)
  }
  if(method=="root_square")#square root absolute correlation distance
  {
    corrMatrix=cor(texpr)
    distanceMatrices=1-corrMatrix^2
    distanceMatrices=sqrt(distanceMatrices)
  }

  if(method=="spearman")#absolute Pearson correlation distance
  {
    corrMatrix=cor(texpr,method="spearman")
    distanceMatrices=1-abs(corrMatrix)
  }
  if(method=="root_spearman")#square root absolute correlation distance
  {
    corrMatrix=cor(texpr,method="spearman")
    distanceMatrices=1-abs(corrMatrix)
    distanceMatrices=sqrt(distanceMatrices)
  }

  if(method=="kendall")#absolute Pearson correlation distance
  {
    corrMatrix=cor(texpr,method="kendall")
    distanceMatrices=1-abs(corrMatrix)
  }

  if(method=="uncentered")#absolute Pearson correlation distance
  {
    corrMatrix=cosine(texpr)
    distanceMatrices=1-abs(corrMatrix)
  }
  if(method=="root_uncentered")#square root absolute correlation distance
  {
    corrMatrix=cosine(texpr)
    distanceMatrices=1-abs(corrMatrix)
    distanceMatrices=sqrt(distanceMatrices)
  }


  #return distance matrices
  distanceMatrices
}


produce_similarity=function(expr,method)
{
  #corr, spearman, cosine,
  #absolute corr
  #CLR debug, MI, MMI?, DI,
  #power adj, WGCNA
  texpr=t(expr)
  similarityMatrices=NULL
  if(method=="correlation")#absolute Pearson correlation distance
  {
    similarityMatrices=cor(texpr)
  }
  similarityMatrices
}

adj_to_edgeFile=function(mat0,geneName,filename,isSymm=T,cutoff=0)
{
  result=c()
  mat0=as.matrix(mat0)
  mat=abs(mat0)
  A=which(mat>cutoff,arr.ind = T)
  for(i in c(1:dim(A)[1]))
  {
    gene1=A[i,1]
    gene2=A[i,2]
    if(isSymm)
    {
      if(gene1>gene2)
      {
        geneName1=geneName[gene1]
        geneName2=geneName[gene2]

        weight=mat[gene1,gene2]
        the_sigh=0
        if(mat0[gene1,gene2]>0)
        {
          the_sigh=1
        }else{
          if(mat0[gene1,gene2]<0)
          {
            the_sigh=-1
          }
        }
        tmp=c(geneName1,geneName2,weight,the_sigh)
        result=rbind(result,tmp)
      }
    }else{
      geneName1=geneName[gene1]
      geneName2=geneName[gene2]
      weight=mat[gene1,gene2]
      tmp=c(geneName1,geneName2,weight,the_sigh)
      result=rbind(result,tmp)
    }


  }
  colnames(result)=c("source","target","weight","sigh")
  write.table(result,file=filename,quote=F,row.names = F,sep=",")
}

run_calculate_similarity=function(inputFile,outputFile,similarityType, keep_whole_adj,cutoff)
{
    data = as.matrix(read.csv(inputFile,header=T,row.names=1))
    similarityMatrices=produce_similarity(data,similarityType)
    if(keep_whole_adj=="T")
    {
      write.table(similarityMatrices,file=outputFile,quote=F,row.names=T,col.names=T,sep=",")
    }else{
        adj_to_edgeFile(similarityMatrices,colnames(data),outputFile,isSymm=T,cutoff=cutoff)
    }

}

run_calculate_distance=function(inputFile,outputFile,distanceType, keep_whole_adj)
{
    data = read.csv(inputFile,header=T,row.names=1)
    distanceMatrices=produceDistanceMatrices(data,distanceType)
    if(keep_whole_adj=="T")
    {
      write.table(distanceMatrices,file=outputFile,quote=F,row.names=T,col.names=T,sep=",")
    }

}

run_calculate_similarity_and_distance=function(inputFile,output_similarity_file,output_distance_file,similarityType,distanceType, keep_whole_adj,cutoff)
{
    print("running")
    print("outfile")
    print(output_similarity_file)
    print("inputFile")
    print(inputFile)
    print("keep_whole_adj")
    print(keep_whole_adj)

    run_calculate_similarity(inputFile,output_similarity_file,similarityType,keep_whole_adj,cutoff)
    run_calculate_distance(inputFile,output_distance_file,distanceType, keep_whole_adj)

}


############################end define functions##############################

args = commandArgs(TRUE)
inputFile = args[1]
output_similarity_file = args[2]
output_distance_file = args[3]
similarityType = args[4]
distanceType = args[5]
keep_whole_adj = args[6]
cutoff = args[7]

options(stringsAsFactors=F)

run_calculate_similarity_and_distance(inputFile,output_similarity_file,output_distance_file,similarityType,distanceType, keep_whole_adj,cutoff)
