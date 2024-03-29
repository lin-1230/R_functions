seuratUmap <-  function(seuratOb,
                       pdfOutDir='',
                       #RDSDir='',
                       umapResolution=0.5,
                       autoSelectPcaDim=T,
                       selectPcaDim=0,
                       nfeatures=2000,
                       dimShold=0.1){

  ## 这是一个专门用来跑Seurat对象标准化到
  ## UMAP流程的函数
  
  # umapResolution 数值型  Seurat umap降纬使用的分辨率
  # autoSelectPcaDim 布尔值，是否自动选择Seurat umap降纬聚类使用的纬度
  # selectPcaDim 数值型，手动选择纬度，只有在autoSelectPcaDim=F的时候起作用
  # nfeatures 数值型，Seurat选择可变异基因的数量
  # dimShold 数值型，自动选择纬度时，前后变化的阈值
  
  
  library(stringr)
  library(Seurat)
  
  seuratOb <- NormalizeData(seuratOb)
  all.genes <- rownames(seuratOb)
  seuratOb <- ScaleData(seuratOb, features = all.genes)
  
  ## 识别高变异的基因(默认2000)
  seuratOb <- FindVariableFeatures(seuratOb, selection.method = "vst", nfeatures = nfeatures)
  ## PCA
  seuratOb <- RunPCA(seuratOb, features = VariableFeatures(object = seuratOb))
  
  ## 当开启自动选择纬度的时候，选择前后两次变化小于阈值（默认0.1）的纬度
  if (autoSelectPcaDim){
    n <- length(seuratOb@reductions[["pca"]]@stdev)
    a <- seuratOb@reductions[["pca"]]@stdev[1:n-1]-seuratOb@reductions[["pca"]]@stdev[2:n]
    selectPcaDim <- which(a < dimShold)[1]
  }
  print(paste0('选择的PCA纬度是：',selectPcaDim))
  
  
  
  ## UMAP降纬
  seuratOb <- FindNeighbors(seuratOb, dims = 1:selectPcaDim)
  seuratOb <- FindClusters(seuratOb, resolution = umapResolution)
  seuratOb <- RunUMAP(seuratOb, dims = 1:selectPcaDim)
  
  ## @date 2024/03/29 由于保存后的RDS太大了
  ## 效率不如重新跑一遍，所以算了
  #print('开始保存umap降纬后的单细胞数据')
  #saveRDS(seuratOb,RDSDir)
  #print('保存结束')
  
  pdf(pdfOutDir,width = 15)
  plot1 <- DimPlot(seuratOb, reduction = "umap")
  print(plot1)
  dev.off()
  
  
  return(seuratOb)
}
