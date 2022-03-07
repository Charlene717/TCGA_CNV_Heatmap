## Ref: https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  ##### 
  #install.packages("tidyverse")
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(eoffice) # Export plot to PPT
  
##### Function setting #####
  source("FUN_Beautify_ggplot.R")
  source("FUN_ggPlot_vline.R")
  source("FUN_TOP_CNV.R")
  source("FUN_GeneExp_Group.R")
  
##### Files setting and import * ##### 
  ## File setting*
  RNAFileName <- "TCGA_PAAD_HiSeqV2"
  CNVFileName <- "TCGA_Gistic2_CopyNumber_Gistic2_all_data_by_genes"
  
  
##### Conditions setting* ##### 
  Target_gene_name <- "TOP2A"
  Mode_Group <- list(Mode="Mean",SD=1) # Mode_Group <- list(Mode="Quartile",Q2="Only")
  
##### Current path and new folder setting* ##### 
  Result_Folder_Name <- paste0(Target_gene_name,"_",Sys.Date()) ## Generate output folder automatically
  dir.create(Result_Folder_Name)
  
##### Load and prepossess CNV data ##### 
  CNV.df <- read.delim(CNVFileName,
                      header = F,sep = "\t")
  
  colnames(CNV.df) <- CNV.df[1,]
  row.names(CNV.df) <- CNV.df[,1]
  CNV.df <- CNV.df[-1,-1]
  Colname <- colnames(CNV.df)
  Rowname <- row.names(CNV.df)
  
  CNV.df <- lapply(CNV.df,as.numeric) %>% as.data.frame()
  colnames(CNV.df) <- Colname
  row.names(CNV.df) <- Rowname

  rm(Colname,Rowname)
  
##### Classify the TOP CNV data #####       
  CNV_Top_Sum.lt <- list()
  
  # Total
  CNV_Top.lt <- TOP_CNV(CNV.df,CNVmode="Total")
  CNV_Top.df <- CNV.df[CNV_Top.lt[["Gene"]],CNV_Top.lt[["Sample"]]]
  
  # Dup
  CNV_Top_Dup.lt <- TOP_CNV(CNV.df,CNVmode="Dup")
  CNV_Top_Dup.df <- CNV.df[CNV_Top_Dup.lt[["Gene"]],CNV_Top_Dup.lt[["Sample"]]]
  
  # Del
  CNV_Top_Del.lt <- TOP_CNV(CNV.df,CNVmode="Del")
  CNV_Top_Del.df <- CNV.df[CNV_Top_Del.lt[["Gene"]],CNV_Top_Del.lt[["Sample"]]]
  
  # TOP Dup+Del
  CNV_Top_2D.df <- rbind(CNV_Top_Dup.df,CNV_Top_Del.df)

  
##### Primary CNV Heatmap #####  
  Heatmap(CNV_Top.df, use_raster=T)

  library(circlize)
  col_fun = colorRamp2(c(min(CNV_Top.df), 0, max(CNV_Top.df)), 
                       c( "#2776e6","white", "#db3784"))
  col_fun = colorRamp2(c(min(CNV_Top.df), 0, max(CNV_Top.df)), 
                       c( "#2776e6","#0c1829", "#db3784"))
  
  col_fun = colorRamp2(c(min(CNV_Top.df), 0, max(CNV_Top.df)), 
                       c( "#02994d","#0c1829", "#e81c4b"))
  col_fun(seq(-3, 3))
  
  Heatmap(CNV_Top.df, name = "Num", col = col_fun, 
          show_column_names = F)
  Heatmap(CNV_Top_Dup.df, name = "Num", col = col_fun, 
          show_column_names = F)
  Heatmap(CNV_Top_Del.df, name = "Num", col = col_fun, 
          show_column_names = F)
  
  col_fun = colorRamp2(c(min(CNV_Top_2D.df), 0, max(CNV_Top_2D.df)), 
                       c("#02994d","#0c1829", "#e81c4b"))
  col_fun(seq(-3, 3))
  Heatmap(CNV_Top_2D.df, name = "Num", col = col_fun, 
          show_column_names = F)
  
  
  
##### Group Samples by RNA expression #####  
  ##### Import genetic data file #####
    GeneExp.df <- read.table(RNAFileName, 
                             header=T, row.names = 1, sep="\t")
    GeneExp_colname <- read.table(RNAFileName, 
                       header=F, row.names = 1, sep="\t")
    colnames(GeneExp.df) <-  GeneExp_colname[1,]
    rm(GeneExp_colname)
    
  ##### Extract Target gene and Statistics ####
    # Extract data with Target_gene_name
    Target_gene_Mean <- GeneExp.df[Target_gene_name,] %>%
                        as.numeric() %>% mean()
    
    #rowMeans(data.matrix(Target_gene))
    Target_gene_SD <- GeneExp.df[Target_gene_name,] %>%
                      as.numeric() %>% sd()
    
    # Quartile
    Target_gene_Q <- GeneExp.df[Target_gene_name,] %>%
                     as.numeric() %>% quantile()

    GeneExp_Group.lt <- GeneExp_Group(GeneExp.df, Target_gene_name, 
                                      Mode="Mean",SD=1)
    
##### CNV Heatmap of Gene expression group #####  
    # TGH: Target gene high
    # TGL: Target gene low
    # TGM: Target gene medium
    
    # GeneExp_Group.lt[["GeneExp_medium.set"]] <- setdiff(colnames(GeneExp.df),
    #                                             c(GeneExp_Group.lt[["GeneExp_high.set"]],GeneExp_Group.lt[["GeneExp_low.set"]]))
    
    TarGeGroup.set <- c(GeneExp_Group.lt[["GeneExp_high.set"]], GeneExp_Group.lt[["GeneExp_low.set"]], 
                        GeneExp_Group.lt[["GeneExp_medium.set"]])
    Anno_GeneExp.df <- data.frame(Sample = TarGeGroup.set,
                                  TarGene = c(rep("High",length(GeneExp_Group.lt[["GeneExp_high.set"]])),
                                              rep("Low",length(GeneExp_Group.lt[["GeneExp_low.set"]])),
                                              rep("Med",length(GeneExp_Group.lt[["GeneExp_medium.set"]]))))

    
    CNV_TGH.df <- CNV.df[,colnames(CNV.df) %in% TarGeGroup.set] 
    Anno_CNV.df <- data.frame(Sample = colnames(CNV_TGH.df))
    Anno_CNV.df <- left_join(Anno_CNV.df,Anno_GeneExp.df)
      
    # Total
    CNV_Top.lt <- TOP_CNV(CNV.df, CNVmode="Total",TopNGene = 2000)
    CNV_Top.df <- CNV.df[CNV_Top.lt[["Gene"]],CNV_Top.lt[["Sample"]]]
    CNV_TGH_Top.df <- CNV_Top.df[,colnames(CNV.df) %in% TarGeGroup.set]
    CNV_TGH_Top.df <- CNV_Top.df[,Anno_CNV.df$Sample]
    
    ## Check
    colnames(CNV_TGH_Top.df) == Anno_CNV.df$Sample
    
    # Dup
    CNV_Top_Dup.lt <- TOP_CNV(CNV.df, CNVmode="Dup",TopNGene = 2000)
    CNV_Top_Dup.df <- CNV.df[CNV_Top_Dup.lt[["Gene"]],CNV_Top_Dup.lt[["Sample"]]]
    CNV_TGH_Top_Dup.df <- CNV_Top_Dup.df[,colnames(CNV.df) %in% TarGeGroup.set]
    
    
    # Del
    CNV_Top_Del.lt <- TOP_CNV(CNV.df, CNVmode="Del",TopNGene = 2000)
    CNV_Top_Del.df <- CNV.df[CNV_Top_Del.lt[["Gene"]],CNV_Top_Del.lt[["Sample"]]]
    CNV_TGH_Top_Del.df <- CNV_Top_Del.df[,colnames(CNV.df) %in% TarGeGroup.set]
    
    # TOP Dup+Del
    CNV_TGH_Top_2D.df <- rbind(CNV_TGH_Top_Dup.df, CNV_TGH_Top_Del.df)
    
    ##### CNV Heatmap #####
      col_fun = colorRamp2(c(min(CNV_TGH_Top.df), 0, max(CNV_TGH_Top.df)), 
                           c( "#02994d","#0c1829", "#e81c4b"))
      col_fun = colorRamp2(c(min(CNV_TGH_Top_2D.df), 0, max(CNV_TGH_Top_2D.df)), 
                           c( "#02994d","#0c1829", "#e81c4b"))
      col_fun(seq(-3, 3))
    
      Heatmap(CNV_TGH_Top.df, name = "Num", col = col_fun, 
              show_column_names = F)
      Heatmap(CNV_TGH_Top_Dup.df, name = "Num", col = col_fun, 
              show_column_names = F)
      Heatmap(CNV_TGH_Top_Del.df, name = "Num", col = col_fun, 
              show_column_names = F)
      Heatmap(CNV_TGH_Top_2D.df, name = "Num", col = col_fun, 
              show_column_names = F)
      
      col_fun = colorRamp2(c(min(CNV_TGH_Top.df), 0, max(CNV_TGH_Top.df)), 
                           c( "#02994d","#0c1829", "#e81c4b"))
      column_ha = HeatmapAnnotation(TarGene = Anno_CNV.df$TarGene,
                  col = list(TarGene = c("High" = "#e04f70", "Low" = "#4474db" ,"Med" ="#adadad")))
      Heatmap(CNV_TGH_Top.df , name = "Num", col = col_fun, 
              show_column_names = F,show_row_names = F, 
              cluster_columns = F, top_annotation = column_ha)
      
      # CNV_TGH_Top_2D.df
      col_fun = colorRamp2(c(min(CNV_TGH_Top_2D.df), 0, max(CNV_TGH_Top_2D.df)), 
                           c( "#02994d","#0c1829", "#e81c4b"))
      column_ha = HeatmapAnnotation(TarGene = Anno_CNV.df$TarGene,
                                    col = list(TarGene = c("High" = "#e04f70", "Low" = "#4474db" ,"Med" ="#adadad")))
      Heatmap(CNV_TGH_Top_2D.df , name = "Num", col = col_fun, 
              show_column_names = F,show_row_names = F, 
              cluster_columns = T, top_annotation = column_ha)
      
      
    # ## Annotation
    # 
    # set.seed(123)
    # mat = matrix(rnorm(100), 10)
    # rownames(mat) = paste0("R", 1:10)
    # colnames(mat) = paste0("C", 1:10)
    # column_ha = HeatmapAnnotation(foo1 = runif(10), bar1 = anno_barplot(runif(10)))
    # row_ha = rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(runif(10)))
    # Heatmap(mat, name = "mat", top_annotation = column_ha, right_annotation = row_ha)
    
    