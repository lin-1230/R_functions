## 设置分组信息
limmaDE <- function(data,
                    dataLabel,
                    group1,
                    group2){
  
  ## data,表达谱
  ## dataLabel，character，data的标签向量，按照样本顺序
  ## group1，group2，character，dataLabel中的两种组别名称
  ## 记住，这里做的差异是group1比group2，即如果
  ## group1是DLBCL，group2是Normal
  ## 则结果中的上调就是在DLBCL中高表达的
  
  library(edgeR)
  
  ## 制作design矩阵
  group <- as.factor(dataLabel)
  design <- model.matrix(~0+group)
  row.names(design) <- colnames(data)
  
  
  ## limma中makeContrasts的contrasts参数要求这种格式
  group1 <- paste0('group',group1)
  group2 <- paste0('group',group2)
  contrasts <- paste0(group1,'-',group2)
  ## 开始做差异
  DGElist <- DGEList(counts = data, group = group )
  keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2 
  DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
  DGElist <- calcNormFactors( DGElist )
  v <- voom(DGElist, design, plot = F, normalize = "quantile")
  fit <- lmFit(v, design)
  cont.matrix <- makeContrasts(contrasts = c(contrasts), levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  nrDEG_limma_voom = topTable(fit2, coef = contrasts, n = Inf)
  nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
  nrDEG_limma_voom <- nrDEG_limma_voom[nrDEG_limma_voom[,5] < 0.05,]
  nrDEG_limma_voom <- nrDEG_limma_voom[order(nrDEG_limma_voom[,1]),]
  return(nrDEG_limma_voom)
}