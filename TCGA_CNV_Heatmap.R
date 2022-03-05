## Ref: https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  ##### 
  #install.packages("tidyverse")
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  
##### Load data ##### 
  TCGA_CNV.df <- read.delim("TCGA_Gistic2_CopyNumber_Gistic2_all_data_by_genes",
                      header = F,sep = "\t")
  
  colnames(TCGA_CNV.df) <- TCGA_CNV.df[1,]
  row.names(TCGA_CNV.df) <- TCGA_CNV.df[,1]
  TCGA_CNV.df <- TCGA_CNV.df[-1,-1]
  Colname <- colnames(TCGA_CNV.df)
  Rowname <- row.names(TCGA_CNV.df)
  
  TCGA_CNV.df <- lapply(TCGA_CNV.df,as.numeric) %>% as.data.frame()
  colnames(TCGA_CNV.df) <- Colname
  row.names(TCGA_CNV.df) <- Rowname
    
  # ## Check ##  
  # SumTestC1 <- sum(abs(TCGA_CNV.df[,1]))
  # SumTestR1 <- sum(abs(TCGA_CNV.df[1,]))
  # SumTestC.All <-apply(TCGA_CNV.df, 2, function(a)sum(abs(a)))
  # SumTestR.All <-apply(TCGA_CNV.df, 1, function(a)sum(abs(a)))
  
  TCGA_CNV.df["Sum",] <- apply(TCGA_CNV.df, 2, function(a)sum(abs(a)))
  TCGA_CNV.df[,"Sum"] <- apply(TCGA_CNV.df, 1, function(a)sum(abs(a)))
  TCGA_CNV_Top.df <- TCGA_CNV.df %>% arrange(desc(Sum)) 
  TCGA_CNV_Top.df <- TCGA_CNV_Top.df[2:51,1:(ncol(TCGA_CNV_Top.df)-1)]
    
##### Primary Heatmap #####  
  Heatmap(TCGA_CNV_Top.df, use_raster=T)
  
  library(circlize)
  col_fun = colorRamp2(c(min(TCGA_CNV_Top.df), 0, max(TCGA_CNV_Top.df)), 
                       c( "#2776e6","white", "#db3784"))
  col_fun = colorRamp2(c(min(TCGA_CNV_Top.df), 0, max(TCGA_CNV_Top.df)), 
                       c( "#2776e6","#0c1829", "#db3784"))
  col_fun = colorRamp2(c(min(TCGA_CNV_Top.df), 0, max(TCGA_CNV_Top.df)), 
                       c( "#02994d","#0c1829", "#e81c4b"))
 
  col_fun(seq(-3, 3))
  
  Heatmap(TCGA_CNV_Top.df, name = "Num", col = col_fun, 
          show_column_names = F)
  
##### Group Samples by RNA expression #####  
  
  

  
  
  
    