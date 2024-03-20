## 定义一个类
#setClass('plotHeatMap',representation(name='character'))

# 定义一个通用函数
#setGeneric("rePlotFunctionHeatMap", function(object) {
#  standardGeneric("rePlotFunctionHeatMap")
#})



#setMethod("rePlotFunctionHeatMap", "plotHeatMap", function() {
  
#})



rePlotFunctionHeatMap <- function(inputFiles,
                                  outputFiles,
                                  #selectFunction,
                                  height=3,
                                  width=10){
  library(data.table)
  library(ggplot2)
  
  data <- fread(inputFiles)
  plotData <- data
  #plotData <- data[Description %in% selectFunction]
  plotData$p <- -1*log10(plotData$qvalue+0.00000000000001)
  setorder(plotData,p)
  plotData$Description <- factor(plotData$Description,levels = unique(plotData$Description))
  plotData <- plotData[!duplicated(plotData$Description),] 
  
  pdf(outputFiles,height = height,width = width)
  p1 <- ggplot(plotData, aes(x = Description, y = p)) +
    geom_bar(stat = "identity",fill='#f09e45') +
    coord_flip()+
    labs(y= '-log10(p)')+
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
