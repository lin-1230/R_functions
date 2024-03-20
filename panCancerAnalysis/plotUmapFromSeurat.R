plotUmapFromSeurat <- function(seuratOb,
                               pdfOutDir='',
                               umapResolution=0.5,
                               autoSelectPcaDim=T,
                               selectPcaDim=0,
                               nfeatures=2000,
                               dimShold=0.1){
  
  
  # umap 图片输出位置，包括文件名
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
  
  sampleName <- colnames(seuratOb)
  sampleName <- str_split(sampleName,'_',simplify = T)[,1]
  ## 提取样本编号
  seuratOb$sampleName <- sampleName
  ## 提取病人
  seuratOb$patient <- str_extract_all(sampleName,'\\d+',simplify = T)[,1]
  
  
  pdf(pdfOutDir,width = 15)
  plot1 <- DimPlot(seuratOb, reduction = "umap")
  plot2 <- DimPlot(seuratOb, reduction = "umap",group.by = 'group')
  plot3 <- DimPlot(seuratOb, reduction = "umap",group.by = 'patient')
  print(plot1+plot3+plot2) 
  dev.off()
  
  return(seuratOb)
}
