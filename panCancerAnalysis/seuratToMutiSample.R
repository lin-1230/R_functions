seuratToMitiSample <- function(seuratObject){
  
  ## seuratObject,seurat对象，请保证有patient列
  ## 根据patient列分样本
  
  ## 取出所有病人编号
  patients <- unique(seuratObject$patient)
  
  ## 创建结果List
  matrixList <- list()
  for(i in 1:length(patients)){
    ## 取出每个病人的细胞
    sampleSeurat <- subset(seuratObject,subset = (patient == patients[i]))
    ## 取出矩阵
    matrixList[[i]] <- as.matrix(sampleSeurat@assays$RNA@counts)
  }
  
  print('每个样本的细胞数量：')
  print(table(seuratObject$patient))
  
  ## 返回要输出的List
  return(matrixList)
}