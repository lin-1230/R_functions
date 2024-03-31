
zhuShi <- function(seuratOB,
                   outPutDir,
                   markList,
                   ifFindMark=T,
                   ifSingleR=T,
                   minPct=0.25,
                   logfcThreshold=0.25,
                   returnThresh=0.05,
                   topN=20){
  
  ## seuratOB,seurat对象
  ## minPct，Seurat包中findAllMarker函数中的
  ## min.pct参数，默认0.25
  ## 用来制定至少在多少比例中的细胞中表达的基因才会
  ## 被选为marker基因
  
  ## markList,List,包含每个细胞的marker 的List
  ## 每一列是一种细胞,命名规则统一为：Tcell、Bcell、NKcell
  ## markList创建示例：
  ## markList <- list(Tcell=c("CD3D",'CD3E','CD8A','CD8B','CD7','IL7R','NKG7','GNLY','GZMA','FGFBP2','PRF1'),
  ##                         Bcell=c("MS4A1","CD79A",'CD79B','CD19','IGLL1','IL4R','MZB1','IGHG1','SDC1','JCHAIN','IGHG3'),
  ##                         NKcell=c('KLRC1',"GNLY","NKG7",'FGFBP2','PRF1'),
  ##                         Monocyte=c('LYZ',"CD14","CD68",'FCN1','APOBEC3A','THBS1'),
  ##                         DCs=c('FCER1A','CLEC10A','LILRA4'),
  ##                         Neutrophils=c('MPO','AZU1','ELANE','CST3','FCGR3B','CSF3R','LTF'),
  ##                         Macrophage=c('FCGR3A',"CD163","CD68"),
  ##                         Erythroblast=c('HBA','HBB'),
  ##                         HSCs=c('SPINK2','CRHBP','LAPTM4B','CLEC3B','ATP8B4','CRYGD','ELN','EXD2'))
  
  
  ## outPutDir,character,所有结果输出的位置
  ## ifFindMark、ifSingleR，是否开启findMark、singleR流程
  ## 默认开启
  
  ## logfcThreshold、returnThresh,numeric,Seurat包中findAllMarker函数中的参数
  ## 分别用来限制logfc和p的阈值
  
  ## topN,每个Cluster提取的top Marker基因数量
  
  library(Seurat)
  library(dplyr)
  library(ggplot2)

  
  


##########寻找差异表达基因
  if (ifFindMark){
    markers <- FindAllMarkers(seuratOB, only.pos = TRUE, min.pct = minPct, logfc.threshold = logfcThreshold,return.thresh = returnThresh)
    #每个cluster的top基因
    top <- markers %>%
      group_by(cluster) %>%
      top_n(n = topN, wt = avg_log2FC)
    
    
    pdf(file=paste0(outPutDir,'heatMap.pdf'),width=15,height=10)
    print(DoHeatmap(seuratOB,                             #seurat对象(数据)
                    features = top$gene,            #绘制的基因
                    label = F,                        #图中是否添加label信息
                    slot = "scale.data",              #绘图使用的表达矩阵
                    group.by = 'seurat_clusters',     #分组名称
                    group.bar = T))                 #是否展示bar
    dev.off()
    write.table(top,paste0(outPutDir,'heatMap.xls'),sep='\t',quote = T,row.names=F)
    
  }
  


#####################singleR注释
  if (ifSingleR){
    library(SingleR)
    library(celldex)
    library(Seurat)
    library(pheatmap)
    
    #seuratOB的meta文件，包含了seurat的聚类结果
    meta<-seuratOB@meta.data 
    
    #载入参考数据集
    ref_use <- HumanPrimaryCellAtlasData()
    #SingleR注释
    pred <- SingleR(test=as.matrix(seuratOB@assays$RNA@count),       #输入表达矩阵
                    ref=ref_use,                                #参考数据
                    labels=ref_use$label.fine,                  #标签列
                    clusters = seuratOB$seurat_clusters,         #细胞聚类信息（只有method = "cluster"才使用）
                    method = "cluster")                         #注释类型，按照cluster还是按照cell进行注释
    
    #保存SingleR注释结果
    write.table(pred,paste0(outPutDir,'singleR.xls'),sep='\t',quote = F)
  }


#####################根据marker注释
  
  ## for循环使用marker注释每一种类型的细胞
  for (i in 1:length(markList)){
    ## 取出目前要注释的细胞类型的名称和marker
    cellAnnName <- names(markList)[i]
    markers <- markList[[i]]
    
    ## 画点图和feature图
    DotPlot(seuratOB, features = markers)+coord_flip()+
      theme_bw()+
      theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
      labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
      scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
    
    ggsave(paste0(outPutDir,'dot_markers_',cellAnnName,'.pdf'),width = 8,height = 7)
    FeaturePlot(seuratOB,features = markers,raster=FALSE)
    ggsave(paste0(outPutDir,'feature_markers_',cellAnnName,'.pdf'),width = 15,height = 10)
    
  }
} 