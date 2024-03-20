seficialGeneDiffScoreByTime <- function(seuratOb,
                             mouseType,
                             timeVetors,
                             selectedGenes,
                             #inputFiles,
                             outputFiles){
  
  
  
  ## 处理数据成数据框
  library(stringr)
  library(ggplot2)
  library(ggsignif)
  
  
  DF <- seuratOb@assays$RNA@scale.data
  
  
  pdf(outputFiles)
  print(
    ggplot(scoreDF, aes(x = group, y = score, fill = group))+
      geom_boxplot(fill=c('#dd7051','#597cb2'))+
      labs(y='score')+
      geom_signif(comparisons = list(c(groupOneName,groupTwoName)),  ##添加显著性
                  map_signif_level = TRUE, 
                  textsize = 4, 
                  vjust = -0.5)+
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
  )
    
}
