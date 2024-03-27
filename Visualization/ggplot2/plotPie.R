plotPie <- function(DF,
                    outputDir,
                    top=10){
  ## 这是一个用来简单画饼图的函数
  
  ## DF，data.frame或者data.table对象，
  ## 仅有两列，第一列是名称，第二列是每一种种类的数量
  
  ## outputDir，character，图片输出地址
  ## top,integer,用来制定绘制的饼图最多包含几个种类
  ## 其他的会被归为Others类里，默认是10
  ## 如果是All，则默认绘制全部
  
  if (top!='All' && top >= nrow(DF)){
    print('top不能大于全部类型数量')
    return()
  }
  
  
  
  
  ## 转换类型
  library(data.table)
  DF <- as.data.table(DF)
  
  ## 给数据框重新命名，方便后续处理
  colnames(DF) <- c('class','num')
  
  ## 按照种类数量重新排序
  setorder(DF,-num)
  
  if(top != 'All'){
    a <- DF[(top+1):.N,sum(num)]
    pieData <- DF[1:top,.(class,num)]
    pieData <- as.data.frame(pieData)
    pieData[top+1,] <- c('others',a)
    pieData <- as.data.table(pieData)
  }else{
    ## 当需要绘制全部内容时，无需作出任何变化
    pieData <- DF
  }
  
  ## 将数量列强制为int型
  pieData[,num:=as.integer(num)]
  
  ## 计算占比
  pieData[,zhanBi:=scales::percent( num / sum(num)) ]
  pieData[,labs:=paste0(class,'(',zhanBi,')')]
  pieData[,labs:=factor(labs,levels = unique(labs))]
  
  pdf(outputDir)
  p1 <- ggplot(pieData,aes(x = '', y = num, fill = labs)) + 
    geom_bar(stat = 'identity', width = 1) + 
    
    theme_bw() + 
    labs(x = '', y = '',title = '') #清除x-y轴的标题，设置主标签。
  
  p2 <- p1 + coord_polar(theta = 'y', start = 0, direction = -1) #direction = 1,顺时针，direction = -1， 逆时针方向。
  print(p2)
  
  dev.off()
  
  
}