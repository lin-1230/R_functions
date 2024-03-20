## 批量matrix转话为Seurat对象
matrixToSeurat <- function(matrixList){
  
  ## matrixList，带有matrix 的list，该函数将list中的matirx
  ## 都转换为Seurat对象
  
  for (i in 1:length(matrixList)){
    ## 如果matrixList第i个元素是空的话，就跳过
    ## 如果不是空的，就创建Seurat对象
    if (!is.null(matrixList[[i]])){
      matrixList[[i]] <- CreateSeuratObject(counts = matrixList[[i]])
    }
  }
  return(matrixList)
}

## 合并Seurat列表中的每一个Seurat对象为一整个Seurat对象
heBing <- function(SeuratList,cellType=''){
  
  ## SeuratList，包含Seurat对象的List
  ## cellType，字符类型，表示细胞的类型，用于Seurat merge函数的add.cell.ids参数，
  ## 用于区分不同样本
  
  ## @date 2024/03/20, 增加了对输入的list只有一个样本的情况
  ## 直接输出List，不需要合并
  if(length(SeuratList) == 1){
    ## @date 2024/03/20, 增加了索引[[1]],这样子返回的就不是list，
    ## 而是Seurat对象
    return(SeuratList[[1]])
  }
  
  
  a <- c()
  for (i in 1:length(SeuratList)) {
    if (!is.null(SeuratList[[i]])){
      a <- c(a,SeuratList[[i]])
    }
  }
  
  x <- a[[1]]
  y <- a[2:length(a)]
  if (cellType == ''){
    ## @date 2024年03月19日,debug, 这里原本的merge(x = SeuratList[[1]]，
    ## 在当第一个样本为NULL的时候会报错，所以进行修改
    ## 修改为：
    SeuratList <- merge(x = x,y = y)
  }
  else{
    cellIds <- paste0(cellType,1:(length(a)))
    ## @date 2024年03月19日,debug, 这里原本的merge(x = SeuratList[[1]]，
    ## 在当第一个样本为NULL的时候会报错，所以进行修改
    ## 修改为：
    SeuratList <- merge(x = x,y = y,add.cell.ids=cellIds)
  }
  
  return(SeuratList)
}



## 将癌症细胞矩阵和中间态细胞矩阵合并，
## 并转化为Seurat对象
heBingMatrixToSeurat <- function(tumorMatrix,
                                 transitionMatrix){
  
  
  ## tumorMatrix，肿瘤细胞矩阵list
  ## transitionMatrix，中间态细胞矩阵list
  
  library(Seurat)
  # 将matrix转化为Seurat对象
  tumorSeuratList <- matrixToSeurat(tumorMatrix)
  transitionSeuratList <- matrixToSeurat(transitionMatrix)
  
  
  ## 使用Seurat包merge函数合并tumor和中间态所有样本
  tumorSeurat <- heBing(tumorSeuratList,'tumor')
  transitionSeurat <- heBing(transitionSeuratList,'transition')
  
  ##  最后合并两个中间态和肿瘤的样本
  seuratObject <- heBing(c(tumorSeurat,transitionSeurat))
  
  seuratObject$group <- c(rep('tumor',ncol(tumorSeurat)),rep('transition',ncol(transitionSeurat)))
  #colnames(seuratObject) <- paste0(seuratObject$group,'_',colnames(seuratObject))
  return(seuratObject)
}

