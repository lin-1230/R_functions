## 筛选细胞
## 将中间态细胞重新分配给肿瘤细胞
## 重新分类细胞身份（tumor or 中间态）
filiterCellsFromSeurat <- function(seuratObject,
                                   method='radio',
                                   threshold=0.7){
  
  ## seuratObject，Seurat对象
  
  ## method，character，radio 是指根据每个簇中肿瘤细胞的比例来决定
  ## 簇中的中间态细胞是否归归入到肿瘤细胞
  ### sample，是指根据每个簇中的细胞是否为同一个样本，如果是的话（认为这些
  ### 细胞根据病人聚类，异质性强，所以认为是肿瘤细胞），这个簇的中间态细胞归为肿瘤细胞
  
  ## @date：2024年03月19号，原radio参数更改为thresholds参数，用来同时表示
  ## 两种方法的阈值
  ## thresholds，numeric，当method = radio时，代表的是每个聚类中，癌细胞占的比例阈值；
  ## 当当method = sample时，代表的是每个聚类中，样本占的比例阈值
  
  
  # @date：2024年03月19号，因为需要增加一种新的筛选细胞的方式
  ## 所以决定新添加一个method参数，用于选择筛选细胞的方式，
  ## 根据选择的不同方法，运行不同的筛选方式
  library(data.table)
  if(method=='radio'){
    metaDF <- as.data.table(seuratObject@meta.data) 
    metaDF2 <- metaDF[,.(s=sum(group=='tumor')/length(group)),by = seurat_clusters]
    
    ## 将中间态细胞分配一下
    metaDF3 <- metaDF2[s > threshold]
    Cluster <- metaDF3$seurat_clusters
    print('符合条件的簇：')
    print(Cluster)
    
    metaDF <- metaDF[seurat_clusters %in% Cluster,group:='tumor']
    seuratObject$group <- metaDF$group
    return(seuratObject)
  }
  else if (method == 'sample'){
    metaDF <- as.data.table(seuratObject@meta.data)
    ## 计算每个簇的细胞数量，并赋值给c列
    ## 计算每个簇中细胞数量最多的病人，并将该病人的细胞数目
    ## 赋值给m
    metaDF2 <- metaDF[,.(c=.N,m=max(table(patient))),by=seurat_clusters]
    ## s列是每个簇中细胞数目最多的样本占所有细胞的比例
    metaDF2$s <- metaDF2$m/metaDF2$c
    
    
    ## 将中间态细胞分配一下
    metaDF3 <- metaDF2[s > threshold]
    Cluster <- metaDF3$seurat_clusters
    print('符合条件的簇：')
    print(Cluster)
    
    ## 将满足上述条件的簇重新分配分组,
    ## 只保留满足条件的簇，并且保留的簇全部归为肿瘤细胞
    metaDF[,group:='']
    metaDF <- metaDF[seurat_clusters %in% Cluster,group:='tumor']
    seuratObject$group <- metaDF$group
    return(seuratObject)
  }
  else{
    print('method参数请选择radio或者sample')
    return(seuratObject)
  }
  
}
