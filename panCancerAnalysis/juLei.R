juLei <- function(SCseuratOb,
                  upGeneBulk,
                  downGeneBulk,
                  bulkSeuratOb,
                  bulkLable,
                  disease,
                  logDir,
                  unClassedCellLoc=NULL,
                  k=3,
                  stable_threshold=0.8,
                  radio_threshold=0.05,
                  nbin = 10){
  
  ## SCseuratOb,单细胞seurat 对象
  
  ## upGeneBulk，character,bulk数据中疾病相对正常上调的基因名称向量，长度为50（没有硬性要求，但是我和师姐是这么讨论的）
  ## downGeneBulk，character,bulk数据中疾病相对正常下调的基因名称向量，长度为50（没有硬性要求，但是我和师姐是这么讨论的）
  
  ## bulkSeuratOb，bulk的seurat对象
  ## bulkLable，character，根据bulk数据样本顺序的样本标签，一般带有疾病缩写和Normal标签
  ## 疾病缩写标签需要与disease参数一致
  
  ## disease，character,疾病缩写，如AML
  ## k,numeric,最后一次聚类的时候需要聚成几类，是kemeans参数，如3
  
  ## 下面两个是循环停止的阈值
  ##stable_threshold,numeric，细胞稳定性的阈值，即前一次聚类和这一次聚类，被分到
  ##同一个簇里的细胞比例，默认0.8
  ##radio_threshold=0.05,细胞比例参数，即这一次循环和上一次循环，分的两个簇比例的变化情况
  ##例如上一次循环比例为0.8个0.2，这一次循环分的是0.79，0.21
  ## 0.8和0.79相差小于0.05，即达到停止条件
  
  ## nbin, numeric,addModuleScore函数的参数，详见该参数的帮助。如：10
  
  ## logDir,character,log文件存放的位置，带文件名称
  
  ## 引用Seurat
  library(Seurat)
  library(edgeR)
  
  ## @date 2024/03/28 将输出重定向到文件中，方便查错
  sink(logDir)
  
  
  upGeneScore <- AddModuleScore(object = SCseuratOb,features = list(upGeneBulk),nbin = nbin)
  downGeneScore <- AddModuleScore(object = SCseuratOb,features = list(downGeneBulk),nbin = nbin)
  
  
  ## 用来存放细胞比例的空向量，后面用来判断循环是否陷入死循环
  rCellVector <- vector()
  rCellVector <- as.numeric(rCellVector)
  
  
  ## 开始循环，直到两群细胞稳定以后结束循环
  i <- 0
  while(T){
    i <- i +1
    print(paste0('-------------第',i,'次循环-----------'))
    
    ## 做差异
    ## 设置分组信息
    ## 算评分
    
    
    ## 使用limma包做差异
    ## 在i=1，也就是第一次循环的时候，使用的是bulk数据的
    ## 差异基因，所以第一次循环不需要做差异
    if (i != 1){
      
      ## limma包差异分析流程
      ## 制作design矩阵
      group <- as.factor(unname(result_p$cluster))
      design <- model.matrix(~0+group)
      row.names(design) <- colnames(SCseuratOb)
      
      DGElist <- DGEList( counts = SCseuratOb@assays$RNA@counts, group = group )
      keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2 
      DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
      DGElist <- calcNormFactors( DGElist )
      v <- voom(DGElist, design, plot = F, normalize = "quantile")
      fit <- lmFit(v, design)
      
      
      cont.matrix <- makeContrasts(contrasts = c('group1-group2'), levels = design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
      nrDEG_limma_voom = topTable(fit2, coef = 'group1-group2', n = Inf)
      nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
      ## 取出矫正后p值小于0.05的结果
      nrDEG_limma_voom <- nrDEG_limma_voom[nrDEG_limma_voom[,5] < 0.05,]
      ## 排序
      nrDEG_limma_voom <- nrDEG_limma_voom[order(nrDEG_limma_voom[,1]),]
      
      
      ## 取出上下调
      nrDEG_limma_voom_up <- nrDEG_limma_voom[nrDEG_limma_voom[,1]>0,]
      nrDEG_limma_voom_down <- nrDEG_limma_voom[nrDEG_limma_voom[,1] < 0,]
      
      
      ## 如果差异基因没有50个，就取所有的显著差异基因
      ## 如果大于50个，就取FC最大的50个
      if (nrow(nrDEG_limma_voom_up) > 50 ){
        upGene <- row.names(nrDEG_limma_voom)[(nrow(nrDEG_limma_voom)-49):nrow(nrDEG_limma_voom)]
      }else{
        upGene <- row.names(nrDEG_limma_voom_up)
      }
      
      if (nrow(nrDEG_limma_voom_down) > 50 ){
        downGene <- row.names(nrDEG_limma_voom)[1:50]
      }else{
        downGene <- row.names(nrDEG_limma_voom_down)
      }
      
      ## 如果没有差异基因，就去掉这个样本
      if (length(downGene) == 0 || length(upGene) == 0){
        return('没有差异基因')
      }
      
      print(upGene)
      print(downGene)
      
      
      
      ## 算评分
      upGeneScore <- AddModuleScore(object = SCseuratOb,features = list(upGene),nbin = nbin)
      downGeneScore <- AddModuleScore(object = SCseuratOb,features = list(downGene),nbin = nbin)
      
      
    }
    
    
    ## Kmeans聚类，根据bulk（第一次循环）或者单细胞（其他循环）差异得到的
    ## 差异基因获得的评分对细胞进行kmeans聚类
    
    ## 将上调打分和下调打分组合成一个数据框
    scoreDF <- data.frame(upGeneScore=upGeneScore$Cluster1,downGeneScore=downGeneScore$Cluster1)
    ## 将细胞名称赋值到新数据框里
    row.names(scoreDF) <- names(upGeneScore$Cluster1)
    ## 用kmeans聚类，将细胞聚成两类
    result <- kmeans(scoreDF, centers = 2)
    ## 取出聚类后，被分到第一类和第二类的细胞的索引
    
    resultC1 <- unname(which(result$cluster == 1))
    resultC2 <- unname(which(result$cluster == 2))
    
    ## 得到两类细胞比例，并获得比例更大的那群细胞的比例
    ## 用于后面判断程序是否陷入死循环的
    ## 死循环情况就是：细胞比例会在几个数字之间来回跳动
    ## 例如：c(0.6,0.2,0.6,0.2,0.6,0.2),
    ## 这种情况说明细胞分群根本没得到进展，只是其中一群细胞
    ## 非常不稳定的在两个聚类中来回跳动
    rCell <- length(resultC1)/(length(resultC1)+length(resultC2))
    rCell <- max(rCell,(1-rCell))
    
    ## 将这次循环的细胞比例加入细胞比例向量
    ## rCellVector 是一个有i长度的向量，也就是i是多少
    ## rCellVector就应该有多长
    rCellVector <- append(rCellVector,rCell)
    
    
    ## 判断两群细胞是否稳定下来了
    ## 第一次循环不需要判断是否结束迭代
    if (i != 1){
      ## 分类,并判断迭代条件
      
      ## 这里是一个之前变量命名不规范的地方，但是现在不好修改了
      ## 说明一下：
      ## resultC1 和resultC2是当前循环的聚类结果
      ## resultPC1 和resultPC2是前一次循环的聚类结果
      ## 取出上一次循环聚类后，被分到第一类和第二类的细胞的索引
      resultPC1 <- unname(which(result_p$cluster == 1))
      resultPC2 <- unname(which(result_p$cluster == 2))
      
      
      ## 上一次循环和这一次循环均有两类
      ## 所以有四种情况，如果将上一次和这一次的聚类结果用c(x,y)
      ## 的形式表示，即x代表上一次的聚类结果，y代表这一次的聚类结果
      ## 则有c(1,1),c(2,1),c(1,2),c(2,2)四种情况
      ## r1-r4是这四种聚类情况的细胞的比例
      r1 <- length(intersect(resultC1,resultPC1))/length(union(resultC1,resultPC1))
      r2 <- length(intersect(resultC1,resultPC2))/length(union(resultC1,resultPC2))
      r3 <- length(intersect(resultC2,resultPC1))/length(union(resultC2,resultPC1))
      r4 <- length(intersect(resultC2,resultPC2))/length(union(resultC2,resultPC2))
      
      
      ## 打印信息，获得反馈
      ##打印这次聚类两类细胞的比例
      print(paste0('细胞比率1:',rCell))
      print(paste0('细胞比率2:',1-rCell))
      ##打印上一次聚类其中一类细胞的比例
      print(paste0('细胞比率_p:',rCell_p))
      ## 打印r1-r4
      print(paste0('r1:',r1))
      print(paste0('r2:',r2))
      print(paste0('r3:',r3))
      print(paste0('r4:',r4))
      
      ## 迭代停止条件：
      ## 一：两次迭代的时候细胞比率变化小
      ## 二：细胞身份的稳定比率达到stable_threshold阈值
      if(abs(rCell-rCell_p) < radio_threshold){ 
        if ((r1 > stable_threshold && r4 > stable_threshold) || (r2 > stable_threshold && r3 > stable_threshold )){
          
          ## 达到迭代停止条件后，最后再将细胞分为3类（默认情况下k=3）
          ## 恶性细胞、正常细胞、中间态细胞（不确定的细胞）
          result <- kmeans(scoreDF, centers = k)
          
          ## 由于这边的up和down并不一定是疾病相对正常
          ## 所以使用对bulk数据进行评分
          ## 来判断哪些基因是疾病上调，哪些是疾病下调
          ## 认为使用最后一次循环的差异基因对bulk数据进行评分
          ## 如果评分结果，上调基因的评分中，疾病细胞的评分高于正常细胞
          ## 则认为这些上调基因就是疾病相对正常的上调
          ## 反之则认为是疾病相对正常的下调
          bulkUpGeneScore <- AddModuleScore(object = bulkSeuratOb,features = list(upGene),nbin = 10)
          bulkDownGeneScore <- AddModuleScore(object = bulkSeuratOb,features = list(downGene),nbin = 10)
          
          
          ## 对疾病和正常细胞的得分分别取平均
          bulkUpGeneScoreMean_disease <- mean(bulkUpGeneScore$Cluster1[bulkLable==disease])
          bulkUpGeneScoreMean_normal <- mean(bulkUpGeneScore$Cluster1[bulkLable=='Normal'])
          print(paste0('bulkUpGeneScoreMean_disease:',bulkUpGeneScoreMean_disease,';bulkUpGeneScoreMean_normal:',bulkUpGeneScoreMean_normal))
          
          ## 如果最后一次循环的上调基因评分中，疾病的评分小于正常细胞的评分
          ## 说明该上调基因实际上应该是疾病相对正常的下调基因
          ## 则修改得分数据库的列名，将上调基因更改为下调
          ## 将下调更改为上调
          if (bulkUpGeneScoreMean_disease < bulkUpGeneScoreMean_normal){
            colnames(scoreDF) <- c('downGeneScore','upGeneScore')
            # @ date 2024/03/21，增加了对kmeans结果（result变量）
            # 的修改，这样子才不会在exacMatrixFromJuLeiRes函数中
            # 导致提取出的正常和疾病矩阵相反
            colnames(result$centers) <- c('downGeneScore','upGeneScore')
          }
          
          ## 返回结果，返回的结果分3个，第一个元素是最后一次循环kMeans的结果
          ## 第二个元素是最后一次循环上下调基因的打分数据框，用于可视化点图
          ## 最后一个元素是当出现死循环时，未能被分类的细胞（中间态）
          ## 不出死循环时，默认为NULL
          return (list(result,scoreDF,unClassedCellLoc))
        }
      }
    }
    
    ## 判断循环是否陷入了来回跳，无法找到最优解的情况
    ## 即有一群细胞每次循环都跑到另外一群细胞中的情况
    ## 如果有，则去掉这群细胞，再跑一次，看看是否可以循环
    ## 当第8次循环以后才开始判断
    if (i > 7){
      print('细胞比率向量:')
      print(rCellVector)
      ## 基本思想
      ## 代码能运行到这里，说明rCellVector中的细胞比例
      ## 两两之间是超过0.05的（默认情况）
      ## 接下来就是判断间隔位置是否相差很小，例如1、3、5位置和2、4、6位置
      ## 每次判断范围为：[i-5:i]
      if (abs(rCellVector[i-1]-rCellVector[i-3]) < radio_threshold && abs(rCellVector[i-3]-rCellVector[i-5]) < radio_threshold){
        if (abs(rCellVector[i]-rCellVector[i-2]) < radio_threshold && abs(rCellVector[i-2]-rCellVector[i-4]) < radio_threshold){
          ## 即说明循环陷入了死循环
          ## 取出反复横跳的那群细胞
          print('出现死循环了')
          
          
          ## 去掉导致死循环的细胞，也就是分类分不出来的细胞
          ## 由于来回跳的情况，就是一群细胞在两个簇中来回跳
          ## 那么，在两次循环中，被分到的簇是不一样的就是这群来回
          ## 跳的细胞
          classedCellLoc <- unname(which(result_p$cluster == result$cluster))
          unClassedCellLoc <- which(!result_p$cluster == result$cluster)
          
          print('去除中间群细胞，再次跑代码')
          
          ## 递归调用自身，再次分配一次细胞身份，这里k为2
          ## 这群不能被分类的细胞，我们认为是中间态细胞
          ## 所以最后循环只要分两次就行
          return (juLei(SCseuratOb[,classedCellLoc],upGeneBulk,downGeneBulk,bulkSeuratOb,bulkLable,disease,unClassedCellLoc,k=2))
          
          
        }
      }
      
      
      
    }
    
    ## 超过20次循环的话，结束迭代
    if ( i > 20){
      return('超出迭代范围')
    }
    
    
    ## result是当前这次循环的结果
    ## result_p是上一次循环的结果
    ## 这次的结果到下一个循环就是前一次循环的结果了
    ## 所以在这里进行变量的更新
    result_p <- result
    rCell_p <- rCell
  }
  
  sink()
  return (list(result_p,scoreDF,unClassedCellLoc))
}
