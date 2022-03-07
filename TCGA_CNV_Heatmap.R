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

##### Files setting and import * ##### 
  ## File setting*
  RNAFileName <- "TCGA_PAAD_HiSeqV2"
  CNVFileName <-"TCGA_Gistic2_CopyNumber_Gistic2_all_data_by_genes"
  
  
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
  TOP_CNV <- function(CNV.df, CNVmode="Total", 
                          TopNGene = 50,TopNSample = ncol(CNV.df)) { 
    # CNVmode = c("Total","Dup","Del")
    CNV_Sum.df <- CNV.df
    
    if(CNVmode=="Dup"){
      CNV_Sum.df["Sum",] <- apply(CNV_Sum.df, 2, function(a)sum(abs(which(a>0))))
      CNV_Sum.df[,"Sum"] <- apply(CNV_Sum.df, 1, function(a)sum(abs(which(a>0))))
      
    }else if(CNVmode=="Del"){
      CNV_Sum.df["Sum",] <- apply(CNV_Sum.df, 2, function(a)sum(abs(which(a<0))))
      CNV_Sum.df[,"Sum"] <- apply(CNV_Sum.df, 1, function(a)sum(abs(which(a<0))))
      
    }else{
      CNV_Sum.df["Sum",] <- apply(CNV_Sum.df, 2, function(a)sum(abs(a)))
      CNV_Sum.df[,"Sum"] <- apply(CNV_Sum.df, 1, function(a)sum(abs(a)))
    }
 
    CNV_Top.df <- CNV_Sum.df %>% arrange(desc(Sum)) 
    
    CNV_Top.df <- CNV_Top.df[,!colnames(CNV_Top.df)==c("Sum")]
    CNV_Top.df <- CNV_Top.df[!row.names(CNV_Top.df)==c("Sum"),]
    
    CNV_Top.df <- CNV_Top.df[1:TopNGene,1:TopNSample]
    
    CNV_Top_Gene.set <- row.names(CNV_Top.df)
    CNV_Top_Sample.set <- colnames(CNV_Top.df)
    
    CNV_Top.lt <- list()
    CNV_Top.lt[["Gene"]] <- CNV_Top_Gene.set
    CNV_Top.lt[["Sample"]] <- CNV_Top_Sample.set
    return(CNV_Top.lt)
  }
  
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
  CNV_Top_All.df <- rbind(CNV_Top_Dup.df,CNV_Top_Del.df)
  
  # ## Old Version ## 
  # # CNV.df["Sum",] <- apply(CNV.df, 2, function(a)sum(abs(a)))
  # # CNV.df[,"Sum"] <- apply(CNV.df, 1, function(a)sum(abs(a)))
  # # CNV_Top.df <- CNV.df %>% arrange(desc(Sum)) 
  # # CNV_Top.df <- CNV_Top.df[2:51,1:(ncol(CNV_Top.df)-1)]
  # 
  # # ## Check ##  
  # # SumTestC1 <- sum(abs(CNV.df[,1]))
  # # SumTestR1 <- sum(abs(CNV.df[1,]))
  # # SumTestC.All <-apply(CNV.df, 2, function(a)sum(abs(a)))
  # # SumTestR.All <-apply(CNV.df, 1, function(a)sum(abs(a)))
  
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
  Heatmap(CNV_Top_All.df, name = "Num", col = col_fun, 
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

  ##### Group the expression matrix according to the expression level of Target gene ####  
    if(Mode_Group$Mode=="Mean"){
      if(Mode_Group$SD==0){
        GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD*(Mode_Group$SD)]
        GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] < Target_gene_Mean-Target_gene_SD*(Mode_Group$SD)]
      }else{
        GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD*(Mode_Group$SD)]
        GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] <= Target_gene_Mean-Target_gene_SD*(Mode_Group$SD)]
      }
      #rm(Target_gene_Mean, Target_gene_SD)
      
    }else{
      if(Mode_Group$Q2=="Only"){ # Mode="Quartile"
        GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Q[3]]
        GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] < Target_gene_Q[3]]
      }else{
        GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Q[4]]
        GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] <= Target_gene_Q[2]]
      }
      #rm(Target_gene_Q)
      
    }
    
    GeneExp_medium.set <- setdiff(colnames(GeneExp.df),
                                  c(GeneExp_high.set,GeneExp_low.set))
    
  ##### Visualization #####  
    ## https://www.jianshu.com/p/9e5b7ffcf80f
    
    data <- reshape2::melt(GeneExp.df[Target_gene_name,]%>%
                             as.numeric())
    TGeneDen.p <- ggplot(data,aes(value,fill=value, color=value)) + 
      xlab("Expression level") + 
      geom_density(alpha = 0.6, fill = "lightgray") + 
      geom_rug() + theme_bw()
    
    ## Set the color
    Mean_SD.clr <- list(rect="#ecbdfc", line="#994db3",text="#6a3b7a" )
    Mean_Q.clr <- list(rect="#abede1", line="#12705f",text="#12705f" )
    
    ## Plot Mean and SD
    TGeneDen.p 
    TGeneDen_SD.p <- ggPlot_vline(TGeneDen.p,data)
    TGeneDen_SD.p+ labs(title= Target_gene_name,
                        x = "Expression level", y = "Density")
    
    ## Plot Quartiles 
    TGeneDen_Q.p <- ggPlot_vline(TGeneDen.p,data,
                                 Line.clr = Mean_Q.clr,
                                 Line1 = Target_gene_Q[2],
                                 Line2 = Target_gene_Q[3],
                                 Line3 = Target_gene_Q[4],
                                 Text.set = c("Q1","Q2","Q3"),
                                 rectP = list(xWidth=0.015, yminP=0.45, ymaxP=0.55,alpha=0.8) 
    )
    
    TGeneDen_Q.p + labs(title= Target_gene_name,
                        x = "Expression level", y = "Density")
    
    ## Plot Quartiles & Mean and SD
    TGeneDen_SD_Q.p <- ggPlot_vline(TGeneDen_SD.p,data,
                                    Line.clr = Mean_Q.clr,
                                    Line1 = Target_gene_Q[2],
                                    Line2 = Target_gene_Q[3],
                                    Line3 = Target_gene_Q[4],
                                    Text.set = c("Q1","Q2","Q3"),
                                    Text.yPos = 0.35,
                                    rectP = list(xWidth=0.015, yminP=0.3, ymaxP=0.4,alpha=0.8) 
    )
    
    TGeneDen_SD_Q.p + labs(title= Target_gene_name,
                           x = "Expression level", y = "Density")
    
    
    pdf(
      file = paste0(Result_Folder_Name,"/",RNAFileName,"_",Target_gene_name,"_DensityPlot.pdf"),
      width = 10,  height = 8
    )
    
    TGeneDen_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) + 
      labs(title= Target_gene_name,
           x ="Expression level", y = "Density")
    TGeneDen_SD.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) + 
      labs(title= Target_gene_name,
           x ="Expression level", y = "Density")
    TGeneDen_SD_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) + 
      labs(title= Target_gene_name,
           x ="Expression level", y = "Density")
    
    dev.off()
    
    # Export PPT
    TGeneDen_SD_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7,
                                        OL_Thick = 1.5) + 
      labs(title= Target_gene_name,
           x ="Expression level", y = "Density") -> TGeneDen_SD_Q2.p
    
    topptx(TGeneDen_SD_Q2.p,paste0(Result_Folder_Name,"/",RNAFileName,"_",Target_gene_name,"_DensityPlot.pptx"))
    
    rm(TGeneDen_SD_Q2.p)
    rm(TGeneDen_Q.p, TGeneDen_SD_Q.p, TGeneDen_SD.p, TGeneDen.p)
    
##### CNV Heatmap of Gene expression group #####  
    # TGH: Target gene high
    # TGL: Target gene low
    # TGM: Target gene medium
    # GeneExp_medium.set <- setdiff(colnames(GeneExp.df),
    #                               c(GeneExp_high.set,GeneExp_low.set))
    
    GeneExp_Anno.df <- data.frame(Sample = c(GeneExp_high.set, GeneExp_low.set, GeneExp_medium.set),
                                  TarGene = c(rep("High",length(GeneExp_high.set)),
                                              rep("Low",length(GeneExp_low.set)),
                                              rep("Med",length(GeneExp_medium.set))))

    
    CNV_TGH.df <- CNV.df[,colnames(CNV.df) %in% 
                           c(GeneExp_high.set,GeneExp_low.set, GeneExp_medium.set)] 
    CNV_Anno.df <- data.frame(Sample = colnames(CNV_TGH.df))
    CNV_Anno.df <- left_join(CNV_Anno.df,GeneExp_Anno.df)
      
    # Total
    CNV_TOP_TGH.lt <- TOP_CNV(CNV_TGH.df, CNVmode="Total")
    CNV_TOP_TGH.df <- CNV.df[CNV_TOP_TGH.lt[["Gene"]],CNV_TOP_TGH.lt[["Sample"]]]
    
    # Dup
    CNV_TOP_Dup_TGH.lt <- TOP_CNV(CNV_TGH.df, CNVmode="Dup")
    CNV_TOP_Dup_TGH.df <- CNV.df[CNV_TOP_Dup_TGH.lt[["Gene"]],
                                 CNV_TOP_Dup_TGH.lt[["Sample"]]]
    
    # Del
    CNV_TOP_Del_TGH.lt <- TOP_CNV(CNV_TGH.df, CNVmode="Del")
    CNV_TOP_Del_TGH.df <- CNV.df[CNV_TOP_Del_TGH.lt[["Gene"]],
                                 CNV_TOP_Del_TGH.lt[["Sample"]]]
    
    # TOP Dup+Del
    CNV_Top_All_TGH.df <- rbind(CNV_TOP_Dup_TGH.df, CNV_TOP_Del_TGH.df)
    
    ##### CNV Heatmap #####
      col_fun = colorRamp2(c(min(CNV_TOP_TGH.df), 0, max(CNV_TOP_TGH.df)), 
                           c( "#02994d","#0c1829", "#e81c4b"))
      col_fun = colorRamp2(c(min(CNV_Top_All_TGH.df), 0, max(CNV_Top_All_TGH.df)), 
                           c( "#02994d","#0c1829", "#e81c4b"))
      col_fun(seq(-3, 3))
    
      Heatmap(CNV_TOP_TGH.df, name = "Num", col = col_fun, 
              show_column_names = F)
      Heatmap(CNV_TOP_Dup_TGH.df, name = "Num", col = col_fun, 
              show_column_names = F)
      Heatmap(CNV_TOP_Del_TGH.df, name = "Num", col = col_fun, 
              show_column_names = F)
      Heatmap(CNV_Top_All_TGH.df, name = "Num", col = col_fun, 
              show_column_names = F)
      
      col_fun = colorRamp2(c(min(CNV_TOP_TGH.df), 0, max(CNV_TOP_TGH.df)), 
                           c( "#02994d","#0c1829", "#e81c4b"))
      column_ha = HeatmapAnnotation(TarGene = CNV_Anno.df$TarGene,
                  col = list(TarGene = c("High" = "#e04f70", "Low" = "#4474db" ,"Med" ="#adadad")))
      Heatmap(CNV_TOP_TGH.df , name = "Num", col = col_fun, 
              show_column_names = F, top_annotation = column_ha)
      
      
      
    ## Annotation
    
    set.seed(123)
    mat = matrix(rnorm(100), 10)
    rownames(mat) = paste0("R", 1:10)
    colnames(mat) = paste0("C", 1:10)
    column_ha = HeatmapAnnotation(foo1 = runif(10), bar1 = anno_barplot(runif(10)))
    row_ha = rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(runif(10)))
    Heatmap(mat, name = "mat", top_annotation = column_ha, right_annotation = row_ha)
    
    