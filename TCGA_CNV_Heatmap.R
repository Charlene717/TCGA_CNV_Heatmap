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
##### Classify the CNV data #####       
  ClassifyCNV <- function(CNV.df, CNVmode="Total", 
                          TopNGene = 50,TopNSample = ncol(CNV.df)) { 
    # mode = c("Total","Dup","Del")
    if(CNVmode=="Total"){
      CNV_Sum.df <- CNV.df
    }
    
    CNV_Sum.df["Sum",] <- apply(CNV_Sum.df, 2, function(a)sum(abs(a)))
    CNV_Sum.df[,"Sum"] <- apply(CNV_Sum.df, 1, function(a)sum(abs(a)))
    
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
  CNV_Top.lt <- ClassifyCNV(CNV.df)
  CNV_Top.df <- CNV.df[CNV_Top.lt[["Gene"]],CNV_Top.lt[["Sample"]]]
  
  ## Old Version ## 
  # CNV.df["Sum",] <- apply(CNV.df, 2, function(a)sum(abs(a)))
  # CNV.df[,"Sum"] <- apply(CNV.df, 1, function(a)sum(abs(a)))
  # CNV_Top.df <- CNV.df %>% arrange(desc(Sum)) 
  # CNV_Top.df <- CNV_Top.df[2:51,1:(ncol(CNV_Top.df)-1)]
  
  # ## Check ##  
  # SumTestC1 <- sum(abs(CNV.df[,1]))
  # SumTestR1 <- sum(abs(CNV.df[1,]))
  # SumTestC.All <-apply(CNV.df, 2, function(a)sum(abs(a)))
  # SumTestR.All <-apply(CNV.df, 1, function(a)sum(abs(a)))
  
##### Primary Heatmap #####  
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

  
  ##### Group Samples by RNA expression #####  
    ##### Import genetic data file #####
      GeneExp.df <- read.table(RNAFileName, 
                               header=T, row.names = 1, sep="\t")
  
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
      
      
    
    
    
      