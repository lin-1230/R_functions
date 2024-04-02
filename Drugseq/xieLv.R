## 计算斜率
computeSlope <- function(seuratOb,mouseType,
                         geneDF,
                         timeSortedVetors,
                         #col,
                         k=10,
                         xd=0.001){
  
  library(stringr)
  library(mgcv)
  
  
  slope <- matrix(nrow=ncol(geneDF),ncol = length(timeSortedVetors))
  colnames(slope) <- timeSortedVetors
  row.names(slope) <- colnames(geneDF)
  for (i in 1:ncol(geneDF)){
    geneList <- geneDF[,i]
    geneList <- geneList[geneList != '']
    geneList <- list(geneList)
    score <- AddModuleScore(object = seuratOb,features = geneList)
    
    metaDF <- as.data.frame(score@meta.data)
    metaDF$hours <- as.numeric(str_replace(metaDF$orig.ident,pattern = 'H',''))
    
    
    mod <- gam(Cluster1 ~ s(hours, bs = "cs", fx = TRUE, k = 5), data = metaDF)
    for ( j in 1:length(timeSortedVetors)) {
      x1 <- timeSortedVetors[j]
      #x1 <- as.numeric(x1)[1]
      y1 <- predict.gam(object = mod,newdata=data.frame(hours=x1))
      y2 <- predict.gam(object = mod,newdata=data.frame(hours=x1+xd))
      
      slope[i,j] <- (y2-y1)/xd
    }
  }
  return (slope)
}


## 绘制斜率函数
slopeZheXianPlot <- function(slope,mouseType,fileName,pdfWidth){
  library(data.table)
  library(ggplot2)
  time <- rep(colnames(slope),each=nrow(slope))
  group <- rep(row.names(slope),ncol(slope))
  value <- as.numeric(slope)
  
  zhexianDf <- data.table(value=value,time=time,group=group)
  zhexianDf$time <- as.numeric(zhexianDf$time)
  
  pdf(paste0('Rresult/',mouseType,'/干性打分/slope/',fileName,'.pdf'),width = pdfWidth)
  p1 <- ggplot(data = zhexianDf,aes(x=time,y=value,group=group,color=group))+
    geom_line(size=1,linetype=1)+
    scale_fill_brewer(palette = "Set3")+
    scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 2))+
    theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
          plot.title = element_text(hjust = 0.5),
          axis.line=element_line(colour="black",linewidth=0.25),
          axis.title=element_text(size=13,face="plain",color="black"),
          axis.text = element_text(size=13,face="plain",color="black"),
          legend.text=element_text(face="plain", colour="black",size=13),
          legend.title=element_text(face="plain", colour="black",size=13),
          panel.grid.major = element_blank(),   
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(hjust = 0.5))
  
  print(p1)
  dev.off()
}