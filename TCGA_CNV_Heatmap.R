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
  
##### Conditions setting* ##### 
  Target_gene_name <- "TOP2A"
  Mode_Group <- list(Mode="Mean",SD=1) # Mode_Group <- list(Mode="Quartile",Q2="Only")
  HeatmapTopGene =  1000
  
##### Current path and new folder setting* ##### 
  Result_Folder_Name <- paste0(Target_gene_name,"_",Sys.Date()) ## Generate output folder automatically
  dir.create(Result_Folder_Name)
  
##### Classify the TOP CNV data #####       
  CNV_Top_Sum.lt <- list()
  
  ### Total
    CNV_Top_Sum.lt[["CNV_Top.lt"]] <- TOP_CNV(CNV.df,CNVmode="Total",TopNGene = HeatmapTopGene)
    CNV_Top_Sum.lt[["CNV_Top.df"]] <- CNV.df[CNV_Top_Sum.lt[["CNV_Top.lt"]][["Gene"]],
                                             CNV_Top_Sum.lt[["CNV_Top.lt"]][["Sample"]]]

  ### Dup
    CNV_Top_Sum.lt[["CNV_Top_Dup.lt"]] <- TOP_CNV(CNV.df,CNVmode="Dup",TopNGene = HeatmapTopGene)
    CNV_Top_Sum.lt[["CNV_Top_Dup.df"]] <- CNV.df[CNV_Top_Sum.lt[["CNV_Top_Dup.lt"]][["Gene"]],
                                                 CNV_Top_Sum.lt[["CNV_Top_Dup.lt"]][["Sample"]]]
  
  ### Del
    CNV_Top_Sum.lt[["CNV_Top_Del.lt"]] <- TOP_CNV(CNV.df,CNVmode="Del",TopNGene = HeatmapTopGene)
    CNV_Top_Sum.lt[["CNV_Top_Del.df"]] <- CNV.df[CNV_Top_Sum.lt[["CNV_Top_Del.lt"]][["Gene"]],
                                                 CNV_Top_Sum.lt[["CNV_Top_Del.lt"]][["Sample"]]]
  
  ### TOP Dup+Del
    CNV_Top_Sum.lt[["CNV_Top_2D.df"]] <- rbind(CNV_Top_Sum.lt[["CNV_Top_Dup.df"]],
                                               CNV_Top_Sum.lt[["CNV_Top_Del.df"]])

  
##### Primary CNV Heatmap #####  
  # Heatmap(CNV_Top_Sum.lt[["CNV_Top.df"]], use_raster=T)

  library(circlize)
  col_fun = colorRamp2(c(min(CNV.df), 0, max(CNV.df)), 
                       c( "#2776e6","white", "#db3784"))
  col_fun = colorRamp2(c(min(CNV.df), 0, max(CNV.df)), 
                       c( "#2776e6","#0c1829", "#db3784"))
  
  col_fun = colorRamp2(c(min(CNV.df), 0, max(CNV.df)), 
                       c( "#02994d","#0c1829", "#e81c4b"))
  col_fun(seq(-3, 3))
  
  pdf(
    file = paste0(Result_Folder_Name,"/CNV_Heatmap_Overal_",Target_gene_name,
                  "_Top",HeatmapTopGene,".pdf"),
    width = 10,  height = 8
  )
  
    Heatmap(CNV_Top_Sum.lt[["CNV_Top.df"]], name = "Num", col = col_fun, 
            show_column_names = F, show_row_names = F) -> Heatmap.p
    
    draw(Heatmap.p, column_title = paste0("Overal CNV heatmaps (TOP",HeatmapTopGene,")"),
         column_title_gp = gpar(fontsize = 16)) #%>% print()
    rm(Heatmap.p)
    # Heatmap(CNV_Top_Sum.lt[["CNV_Top_Dup.df"]], name = "Num", col = col_fun, 
    #         show_column_names = F, show_row_names = F)
    # Heatmap(CNV_Top_Sum.lt[["CNV_Top_Del.df"]], name = "Num", col = col_fun, 
    #         show_column_names = F, show_row_names = F)
    
    Heatmap(CNV_Top_Sum.lt[["CNV_Top_2D.df"]], name = "Num", col = col_fun, 
            show_column_names = F, show_row_names = F) -> Heatmap.p
    draw(Heatmap.p, column_title = paste0("Overal CNV heatmaps (±TOP",HeatmapTopGene,")"),
         column_title_gp = gpar(fontsize = 16)) #%>% print()
    rm(Heatmap.p)
  dev.off()
  
##### Group Samples by RNA expression #####  
  ##### Import genetic data file #####
    GeneExp.df <- read.table(RNAFileName, 
                             header=T, row.names = 1, sep="\t")
    GeneExp_colname <- read.table(RNAFileName, 
                       header=F, row.names = 1, sep="\t")
    colnames(GeneExp.df) <-  GeneExp_colname[1,]
    rm(GeneExp_colname)
    
  ##### Extract Target gene and Statistics ####
    GeneExp_Group.lt <- GeneExp_Group(GeneExp.df, Target_gene_name, 
                                      Mode="Mean",SD=1)
    
    ### Open when the extreme groups(High or Low) are too small
    # GeneExp_Group.lt[["GeneExp_medium.set"]] <- setdiff(colnames(GeneExp.df),
    #                                             c(GeneExp_Group.lt[["GeneExp_high.set"]],GeneExp_Group.lt[["GeneExp_low.set"]]))
    
##### CNV Heatmap of Gene expression group #####  
    # TGH: Target gene high
    # TGL: Target gene low
    # TGM: Target gene medium
    # TGS: Target gene Sum
    
    ##### TGS #####
      TarGeGroup.set <- c(GeneExp_Group.lt[["GeneExp_high.set"]], GeneExp_Group.lt[["GeneExp_low.set"]], 
                          GeneExp_Group.lt[["GeneExp_medium.set"]])
      Anno_GeneExp.df <- data.frame(Sample = TarGeGroup.set,
                                    TarGene = c(rep("High",length(GeneExp_Group.lt[["GeneExp_high.set"]])),
                                                rep("Low",length(GeneExp_Group.lt[["GeneExp_low.set"]])),
                                                rep("Med",length(GeneExp_Group.lt[["GeneExp_medium.set"]]))))
      
      CNV_TGS_Top_Sum.lt <- list()
      
      CNV_TGS_Top_Sum.lt[["CNV_TGS.df"]] <- CNV.df[,colnames(CNV.df) %in% TarGeGroup.set] 
      Anno_CNV.df <- data.frame(Sample = colnames(CNV_TGS_Top_Sum.lt[["CNV_TGS.df"]]))
      Anno_CNV.df <- left_join(Anno_CNV.df,Anno_GeneExp.df) %>% arrange(TarGene)
      

    ### Total
      CNV_TGS_Top_Sum.lt[["CNV_TGS_Top.df"]] <- CNV_Top_Sum.lt[["CNV_Top.df"]][,Anno_CNV.df$Sample]
        ## Check
        colnames(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top.df"]]) == Anno_CNV.df$Sample
  
    ### Dup
      CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_Dup.df"]] <- CNV_Top_Sum.lt[["CNV_Top_Dup.df"]][,Anno_CNV.df$Sample]
      
    ### Del
      CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_Del.df"]] <- CNV_Top_Sum.lt[["CNV_Top_Del.df"]][,Anno_CNV.df$Sample]

    
    ### TOP Dup+Del
      CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_2D.df"]] <- rbind(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_Dup.df"]], 
                                                CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_Del.df"]])
        ## Check
        colnames(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_2D.df"]]) == Anno_CNV.df$Sample
    
    
    ##### CNV Heatmap #####
    ## Ref: https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/ComplexHeatmap/inst/doc/s9.examples.html
      
      # CNV_TGS_Top.df
      column_ha = HeatmapAnnotation(TarGene = Anno_CNV.df$TarGene,
                  col = list(TarGene = c("High" = "#e04f70", "Low" = "#4474db" ,"Med" ="#adadad")))
        ## Modify the label name
        #column_ha@anno_list[["TarGene"]]@name <- Target_gene_name
        column_ha@anno_list[["TarGene"]]@label <- Target_gene_name
        #column_ha@anno_list[["TarGene"]]@name_param[["label"]] <- Target_gene_name
        column_ha@anno_list[["TarGene"]]@color_mapping@name <- Target_gene_name
       
      pdf(
        file = paste0(Result_Folder_Name,"/CNV_Heatmap_TGS_",Target_gene_name,
                      "_Top",HeatmapTopGene,".pdf"),
        width = 10,  height = 8
      )   
      # CNV_TGS_Top.df
        Heatmap(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top.df"]] , name = "Num", col = col_fun, 
                show_column_names = F,show_row_names = F, 
                cluster_columns = T, top_annotation = column_ha)-> Heatmap.p
        draw(Heatmap.p, column_title = paste0("TGS CNV heatmaps (TOP",HeatmapTopGene,")"),
             column_title_gp = gpar(fontsize = 16)) #%>% print()
        rm(Heatmap.p)
        
        Heatmap(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top.df"]] , name = "Num", col = col_fun, 
                show_column_names = F,show_row_names = F, 
                cluster_columns = F, top_annotation = column_ha)-> Heatmap.p
        draw(Heatmap.p, column_title = paste0("TGS CNV heatmaps (TOP",HeatmapTopGene,")"),
             column_title_gp = gpar(fontsize = 16)) #%>% print()
        rm(Heatmap.p)
        
      # CNV_TGS_Top_2D.df
        Heatmap(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_2D.df"]] , name = "Num", col = col_fun, 
                show_column_names = F,show_row_names = F, 
                cluster_columns = T, top_annotation = column_ha)-> Heatmap.p
        draw(Heatmap.p, column_title = paste0("TGS CNV heatmaps (±TOP",HeatmapTopGene,")"),
             column_title_gp = gpar(fontsize = 16)) #%>% print()
        rm(Heatmap.p)
        
        Heatmap(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_2D.df"]] , name = "Num", col = col_fun, 
                show_column_names = F,show_row_names = F, 
                cluster_columns = F, top_annotation = column_ha)-> Heatmap.p
        draw(Heatmap.p, column_title = paste0("TGS CNV heatmaps (±TOP",HeatmapTopGene,")"),
             column_title_gp = gpar(fontsize = 16)) #%>% print()
        rm(Heatmap.p)
        
        rm(column_ha)
      dev.off()
    ##### TGH #####
      ##### TGS #####
      # TarGeGroup.set <- c(GeneExp_Group.lt[["GeneExp_high.set"]], GeneExp_Group.lt[["GeneExp_low.set"]], 
      #                     GeneExp_Group.lt[["GeneExp_medium.set"]])
      # Anno_GeneExp.df <- data.frame(Sample = TarGeGroup.set,
      #                               TarGene = c(rep("High",length(GeneExp_Group.lt[["GeneExp_high.set"]])),
      #                                           rep("Low",length(GeneExp_Group.lt[["GeneExp_low.set"]])),
      #                                           rep("Med",length(GeneExp_Group.lt[["GeneExp_medium.set"]]))))
      # 
      # CNV_TGS_Top_Sum.lt <- list()
      # 
      # CNV_TGS_Top_Sum.lt[["CNV_TGS.df"]] <- CNV.df[,colnames(CNV.df) %in% TarGeGroup.set] 
      # Anno_CNV.df <- data.frame(Sample = colnames(CNV_TGS_Top_Sum.lt[["CNV_TGS.df"]]))
      # Anno_CNV.df <- left_join(Anno_CNV.df,Anno_GeneExp.df) %>% arrange(TarGene)

      ##### Classify the TOP CNV data #####       
        CNV_TGH_Top_Sum.lt <- list()
        CNV_TGH.df <- CNV.df[,GeneExp_Group.lt[["GeneExp_high.set"]]]
      
      ### Total
        CNV_TGH_Top_Sum.lt[["CNV_Top.lt"]] <- TOP_CNV(CNV_TGH.df,CNVmode="Total",TopNGene = HeatmapTopGene)
        CNV_TGH_Top_Sum.lt[["CNV_Top.df"]] <- CNV.df[CNV_TGH_Top_Sum.lt[["CNV_Top.lt"]][["Gene"]],
                                                     colnames(CNV.df) %in% TarGeGroup.set]
        
      ### Dup
        CNV_TGH_Top_Sum.lt[["CNV_Top_Dup.lt"]] <- TOP_CNV(CNV_TGH.df,CNVmode="Dup",TopNGene = HeatmapTopGene)
        CNV_TGH_Top_Sum.lt[["CNV_Top_Dup.df"]] <- CNV.df[CNV_TGH_Top_Sum.lt[["CNV_Top_Dup.lt"]][["Gene"]],
                                                         colnames(CNV.df) %in% TarGeGroup.set]
      
      ### Del
        CNV_TGH_Top_Sum.lt[["CNV_Top_Del.lt"]] <- TOP_CNV(CNV_TGH.df,CNVmode="Del",TopNGene = HeatmapTopGene)
        CNV_TGH_Top_Sum.lt[["CNV_Top_Del.df"]] <- CNV.df[CNV_TGH_Top_Sum.lt[["CNV_Top_Del.lt"]][["Gene"]],
                                                         colnames(CNV.df) %in% TarGeGroup.set]
      
      ### TOP Dup+Del
        CNV_TGH_Top_Sum.lt[["CNV_Top_2D.df"]] <- rbind(CNV_TGH_Top_Sum.lt[["CNV_Top_Dup.df"]],
                                                   CNV_TGH_Top_Sum.lt[["CNV_Top_Del.df"]])
      
      
      # ### Total
      #   CNV_TGS_Top_Sum.lt[["CNV_TGS_Top.df"]] <- CNV_TGH_Top_Sum.lt[["CNV_Top.df"]][,Anno_CNV.df$Sample]
      # ## Check
      #  colnames(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top.df"]]) == Anno_CNV.df$Sample
      # 
      # ### Dup
      #   CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_Dup.df"]] <- CNV_TGH_Top_Sum.lt[["CNV_Top_Dup.df"]][,Anno_CNV.df$Sample]
      # 
      # ### Del
      #  CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_Del.df"]] <- CNV_TGH_Top_Sum.lt[["CNV_Top_Del.df"]][,Anno_CNV.df$Sample]
      # 
      # 
      # ### TOP Dup+Del
      #   CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_2D.df"]] <- rbind(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_Dup.df"]], 
      #                                                    CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_Del.df"]])
      # ## Check
      #   colnames(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_2D.df"]]) == Anno_CNV.df$Sample
      
      
      ##### CNV Heatmap #####
      ## Ref: https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/ComplexHeatmap/inst/doc/s9.examples.html
      
      pdf(
        file = paste0(Result_Folder_Name,"/CNV_Heatmap_TGH_",Target_gene_name,
                      "_Top",HeatmapTopGene,".pdf"),
        width = 10,  height = 8
      )   
        
        # CNV_TGH_Top.df
        column_ha = HeatmapAnnotation(TarGene = Anno_CNV.df$TarGene,
                                      col = list(TarGene = c("High" = "#e04f70", "Low" = "#4474db" ,"Med" ="#adadad")))
        Heatmap(CNV_TGH_Top_Sum.lt[["CNV_Top.df"]] , name = "Num", col = col_fun, 
                show_column_names = F,show_row_names = F, 
                cluster_columns = T, top_annotation = column_ha)
        Heatmap(CNV_TGH_Top_Sum.lt[["CNV_Top.df"]] , name = "Num", col = col_fun, 
                show_column_names = F,show_row_names = F, 
                cluster_columns = F, top_annotation = column_ha)
        
        # CNV_TGH_Top_2D.df
        Heatmap(CNV_TGH_Top_Sum.lt[["CNV_Top_2D.df"]] , name = "Num", col = col_fun, 
                show_column_names = F,show_row_names = F, 
                cluster_columns = T, top_annotation = column_ha)
        Heatmap(CNV_TGH_Top_Sum.lt[["CNV_Top_2D.df"]] , name = "Num", col = col_fun, 
                show_column_names = F,show_row_names = F, 
                cluster_columns = F, top_annotation = column_ha)
        
        rm(column_ha)
        
      dev.off()
      
       
    