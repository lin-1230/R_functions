matrixToCsv <- function(matrixData,
                        outPutDir,
                        fileNamePre){
  
  ## matrixData，List，matrix数据的list
  ## outPutDir，character，存放输出csv的目录
  ## fileNamePre，character，输出csv文件的前缀，
  ## 最后输出的csv名称会是前缀加上数字1、2、3....
  
  ## 遍历list
  for (i in 1:length(matrixData)){
    ## 如果第i个元素是null就跳过
    if (!is.null(matrixData[[i]])){
      write.csv(matrixData[[i]],paste0(outPutDir,fileNamePre,'_',i,'.csv'))
    }
  }
}

RDSToCSV <- function(inputDir,
                     outPutDir){
  
  ## intPutDir，character，存放输入RDS的目录
  ## outPutDir，character，存放输出csv的目录
  
  
  ## 获得inputDir目录下所有文件的地址和名称
  ## (RDS下数据应该为GSE编号_tumoreMatrix.RDS样式)
  RDSFileDir <- list.files(inputDir,full.names = T)
  RDSFileName <- list.files(inputDir)
  
  
  # 设置csv文件前缀
  library(stringr)
  fileNamePre <- str_split(RDSFileName,'_',simplify = T)[,1]
  
  ## 遍历每个文件的地址
  ## 读取并转化为csv，并保存
  for (i in 1:length(RDSFileDir)){
    matrixData <- readRDS(RDSFileDir[i])
    matrixToCsv(matrixData = matrixData,
                outPutDir = outPutDir,fileNamePre = fileNamePre[i])
  }
  
  
  
}
