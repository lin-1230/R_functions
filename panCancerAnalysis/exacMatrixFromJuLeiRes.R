
exacMatrixFromJuLeiRes <- function(juLeiResult,
                                   seuratOB){
  
  ##这是一个从聚类的结果（juLei那个函数）中
  ##分别提取出肿瘤细胞、正常细胞、中间态细胞矩阵的函数
  
 
  #juLeiResult，List，juLei函数返回的结果
  #@ date 2024/03/26 添加
  #seuratOB,seurat对象，需要提取不同细胞的单细胞数据的seurat对象
  
  ## 创建存放三个结果的List
  normalMatrix <- list()
  transitionMatrix <- list()
  tumorMatrix <- list()
  
  ## 遍历每个样本的结果
  for (i in 1:length(juLeiResult)){
    ## @date 2024/03/28 修改bug
    ## 因为if判断，当juLeiResult[[i]]有三个元素，但是判断标准是
    ## !='没有差异基因'，这个时候会报错，所以进行如下修改
    a <- juLeiResult[[i]]
    ## 如果该样本没有结果就跳过,当a的长度大于1的时候
    ## 说明有结果
    if (length(a) > 1){
      print(i)
      
      
      
      ## 如果list有三个元素，说明这个样本是经历过死循环的样本（详见juLei函数）
      ## 有死循环的样本，因为去掉了中间态细胞重新跑了一次，
      ## 所以提取方式有所不同，所以做一个判断
      if (is.null(a[[3]])){
        
        ## 根据Kmeans结果中的centers来判断哪一个簇是肿瘤细胞
        ## 在疾病相对正常上调细胞评分最高的cluster序号就是cancer
        ## 在疾病相对正常上调细胞评分最低的cluster序号就是normal
        ## 剩下的那个就是中间态细胞
        cancerCluster <- which.max(as.data.frame(a[[1]]$centers)$upGeneScore)
        normalCluster <- which.min(as.data.frame(a[[1]]$centers)$upGeneScore)
        transitionCluster <- setdiff(c(1,2,3),c(cancerCluster,normalCluster))
        
        ## 根据序号取出每种类别的表达谱数据
        ## @date 2024/03/26 将代码中使用的data_GSE182434（DLBCL的数据集，应该是忘记更换了）
        ## 替换成了seuratOB
        countData <- seuratOB[[i]]@assays$RNA@counts
        tumorMatrix[[i]] <- as.matrix(countData[,unname(which(a[[1]]$cluster==cancerCluster))])
        normalMatrix[[i]] <- as.matrix(countData[,unname(which(a[[1]]$cluster==normalCluster))])
        transitionMatrix[[i]] <- as.matrix(countData[,unname(which(a[[1]]$cluster==transitionCluster))])
        
      }
      else{
        ## 判断每种类别分别对应了哪个序号
        cancerCluster <- which.max(as.data.frame(a[[1]]$centers)$upGeneScore)
        normalCluster <- which.min(as.data.frame(a[[1]]$centers)$upGeneScore)
        
        ## 根据序号取出每种类别的表达谱数据
        countData <- seuratOB[[i]]@assays$RNA@counts
        
        ## 由于样本中取出了一些中间态细胞，所以不能直接用顺序取细胞
        ## 而需要通过匹配
        tumorMatrix[[i]] <- as.matrix(countData[,match(names(which(a[[1]]$cluster==cancerCluster)),colnames(countData))])
        normalMatrix[[i]] <- as.matrix(countData[,match(names(which(a[[1]]$cluster==normalCluster)),colnames(countData))])
        transitionMatrix[[i]] <- as.matrix(countData[,match(names(a[[3]]),colnames(countData))])
        
        
      }
      
    }
  }
  
  return(list(normalMatrix,transitionMatrix,tumorMatrix))
}