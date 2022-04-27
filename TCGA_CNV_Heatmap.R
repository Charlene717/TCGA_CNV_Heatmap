## Ref: https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html

##### To Do List ######
  # - [ ] Consider Phenotype (Normal,Primary Tumor)
  # - [ ] Design method of clustering (Gene A high expression group and low expression group)
  # - [ ] Extract the gene list by the greatest difference in the comparison of the groups classified by gene expression
  # - [ ] Find nearby genes
  # - [ ] Gene name conversion

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  ##### 
  # install.packages("tidyverse")
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
  RNAFileName <- "RNA_TCGA-PAAD" # RNAFileName <- "TCGA_PAAD_HiSeqV2"
  CNVFileName <- "CNV_TCGA-PAAD.gistic.tsv"
  PhenoFileName <- "TCGA-PAAD.GDC_phenotype.tsv"

##### Conditions setting* ##### 
  Target_gene_name <- "ENSG00000131747" # ENSG00000131747 TOP2A
  Mode_Group <- list(Mode="Mean",SD=1) # Mode_Group <- list(Mode="Quartile",Q2="Only")
  HeatmapTopGene =  1000
  
  # Figure note
  if(Mode_Group[["Mode"]]=="Mean"){
    Note = paste0("Mean_SD",Mode_Group[["SD"]])
  }else{
    Note = paste0("Quarter")
  }
  
  # IntGene.set <- c("TOP2A","NSUN2","MB21D1","TP53","PTK2",
  #                  "PIK3CA","EGFR","MYC","KRAS","BRAF") ## MB21D1= CGAS
  IntGene.set <- c("ENSG00000131747","ENSG00000037474","ENSG00000164430","ENSG00000141510","ENSG00000169398",
                   "ENSG00000121879","ENSG00000146648","ENSG00000136997","ENSG00000133703","ENSG00000157764") ## MB21D1= CGAS
  
##### Current path and new folder setting* ##### 
  Result_Folder_Name <- paste0(Target_gene_name,"_",Sys.Date()) ## Generate output folder automatically
  dir.create(Result_Folder_Name)
  
  
##### Load and prepossess CNV data ##### 
  CNV.df <- read.delim(CNVFileName,
                      header = F,sep = "\t")
  
  Colname <- CNV.df[1,]
  Rowname <- CNV.df[,1]
  CNV.df <- CNV.df[-1,-1]
  
  CNV.df <- lapply(CNV.df,as.numeric) %>% as.data.frame()
  colnames(CNV.df) <- Colname[-1]
  row.names(CNV.df) <- Rowname[-1]

  rm(Colname,Rowname)
  
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

  ## Heatmap color setting
  col_fun = colorRamp2(c(min(CNV.df), 0, max(CNV.df)), 
                       c( "#02994d","#0c1829", "#e81c4b"))
  col_fun = colorRamp2(c(-2, 0, 2), 
                       c( "#02994d","#0c1829", "#e81c4b"))
  
  col_fun(seq(-3, 3))
  
  pdf(
    file = paste0(Result_Folder_Name,"/CNV_Heatmap_Overal_",Target_gene_name,
                  "_Top",HeatmapTopGene,".pdf"),
    width = 10,  height = 8
  )
  
    Heatmap.p <- Heatmap(CNV_Top_Sum.lt[["CNV_Top.df"]], name = "Num", col = col_fun, 
              show_column_names = F, show_row_names = F) 
    
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
  ##### Import RNA expression data #####
    GeneExp.df <- read.table(RNAFileName, 
                             header=T, row.names = 1, sep="\t")
    GeneExp_colname <- read.table(RNAFileName, 
                       header=F, row.names = 1, sep="\t")
    colnames(GeneExp.df) <-  GeneExp_colname[1,]
    rm(GeneExp_colname)
    
    # ##### Covert Gene name #####
    # # https://stackoverflow.com/questions/28543517/how-can-i-convert-ensembl-id-to-gene-symbol-in-r
    # # https://www.rdocumentation.org/packages/biomaRt/versions/2.28.0/topics/useDataset
    # library('biomaRt')
    # mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    # genes <- row.names(GeneExp.df)
    # G_list <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"),values=genes,mart= mart)
    # 
    # ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
    # 
    # my_ensembl_gene_id <- getBM(attributes = 'ensembl_gene_id',
    #                             filters = 'chromosome_name',
    #                             values = genes,
    #                             mart = ensembl
    # )
  ##### Extract Target gene and Statistics ####
    #ENSG# Target_gene_name <- row.names(GeneExp.df)[grepl(Target_gene_name, row.names(GeneExp.df))]
    GeneExp_Group.lt <- GeneExp_Group(GeneExp.df, Target_gene_name, 
                                      Mode="Mean",SD=1,Result_Folder_Name=Result_Folder_Name)
    
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
    
    
      ##### TGS CNV Heatmap #####
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
          draw(Heatmap.p, column_title = paste0("TGS CNV heatmaps (TOP",HeatmapTopGene,")","\n",Note),
               column_title_gp = gpar(fontsize = 16)) #%>% print()
          rm(Heatmap.p)
          
          Heatmap(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top.df"]] , name = "Num", col = col_fun, 
                  show_column_names = F,show_row_names = F, 
                  cluster_columns = F, top_annotation = column_ha)-> Heatmap.p
          draw(Heatmap.p, column_title = paste0("TGS CNV heatmaps (TOP",HeatmapTopGene,")","\n",Note),
               column_title_gp = gpar(fontsize = 16)) #%>% print()
          rm(Heatmap.p)
          
        # CNV_TGS_Top_2D.df
          Heatmap(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_2D.df"]] , name = "Num", col = col_fun, 
                  show_column_names = F,show_row_names = F, 
                  cluster_columns = T, top_annotation = column_ha)-> Heatmap.p
          draw(Heatmap.p, column_title = paste0("TGS CNV heatmaps (±TOP",HeatmapTopGene,")","\n",Note),
               column_title_gp = gpar(fontsize = 16)) #%>% print()
          rm(Heatmap.p)
          
          Heatmap(CNV_TGS_Top_Sum.lt[["CNV_TGS_Top_2D.df"]] , name = "Num", col = col_fun, 
                  show_column_names = F,show_row_names = F, 
                  cluster_columns = F, top_annotation = column_ha)-> Heatmap.p
          draw(Heatmap.p, column_title = paste0("TGS CNV heatmaps (±TOP",HeatmapTopGene,")","\n",Note),
               column_title_gp = gpar(fontsize = 16)) #%>% print()
          rm(Heatmap.p)
          
        dev.off()
        
        
        ## Search IntGene.set
        IntGene.set
        IntGene.set2 <- list()
        for (i in 1:length(IntGene.set)) {
          Gene <- row.names(GeneExp.df)[grepl(IntGene.set[i], row.names(GeneExp.df))]
          IntGene.set2[i] <-Gene
        }
        rm(i)
        IntGene.set <- IntGene.set2 %>% unlist()
          
        row.names(GeneExp.df)[grepl(IntGene.set, row.names(GeneExp.df))]
        
        pdf(
          file = paste0(Result_Folder_Name,"/CNV_Heatmap_IntGene_",Target_gene_name,".pdf"),
          width = 10,  height = 8
        )  
          CNV_IntGene_Top.df <- CNV.df[IntGene.set,Anno_CNV.df$Sample]
          Heatmap(CNV_IntGene_Top.df , name = "Num", col = col_fun, 
                  show_column_names = F,show_row_names = T, 
                  cluster_columns = F,cluster_rows = F ,top_annotation = column_ha)-> Heatmap.p
          draw(Heatmap.p, column_title = paste0("TGS CNV heatmaps IntGene","\n",Note),
               column_title_gp = gpar(fontsize = 16)) #%>% print()
          rm(Heatmap.p)
        dev.off() 
        
        #rm(column_ha) 
        
    ##### TGH #####
      ##### Classify the TOP CNV data #####       
        CNV_TGH_Top_Sum.lt <- list()
        CNV_TGH_Top_Sum.lt[["CNV_TGH.df"]] <- CNV.df[,GeneExp_Group.lt[["GeneExp_high.set"]]]
        CNV_TGH.df <- CNV_TGH_Top_Sum.lt[["CNV_TGH.df"]] 
      
      ### Total
        CNV_TGH_Top_Sum.lt[["CNV_THG_Top.lt"]] <- TOP_CNV(CNV_TGH.df,CNVmode="Total",TopNGene = HeatmapTopGene)
        CNV_TGH_Top_Sum.lt[["CNV_THG_Top.df"]] <- CNV.df[CNV_TGH_Top_Sum.lt[["CNV_THG_Top.lt"]][["Gene"]],
                                                     colnames(CNV.df) %in% TarGeGroup.set]
        
      ### Dup
        CNV_TGH_Top_Sum.lt[["CNV_THG_Top_Dup.lt"]] <- TOP_CNV(CNV_TGH.df,CNVmode="Dup",TopNGene = HeatmapTopGene)
        CNV_TGH_Top_Sum.lt[["CNV_THG_Top_Dup.df"]] <- CNV.df[CNV_TGH_Top_Sum.lt[["CNV_THG_Top_Dup.lt"]][["Gene"]],
                                                         colnames(CNV.df) %in% TarGeGroup.set]
      
      ### Del
        CNV_TGH_Top_Sum.lt[["CNV_THG_Top_Del.lt"]] <- TOP_CNV(CNV_TGH.df,CNVmode="Del",TopNGene = HeatmapTopGene)
        CNV_TGH_Top_Sum.lt[["CNV_THG_Top_Del.df"]] <- CNV.df[CNV_TGH_Top_Sum.lt[["CNV_THG_Top_Del.lt"]][["Gene"]],
                                                         colnames(CNV.df) %in% TarGeGroup.set]
      
      ### TOP Dup+Del
        CNV_TGH_Top_Sum.lt[["CNV_THG_Top_2D.df"]] <- rbind(CNV_TGH_Top_Sum.lt[["CNV_THG_Top_Dup.df"]],
                                                   CNV_TGH_Top_Sum.lt[["CNV_THG_Top_Del.df"]])
      
      
      ##### CNV Heatmap #####
      ## Ref: https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/ComplexHeatmap/inst/doc/s9.examples.html
      
      pdf(
        file = paste0(Result_Folder_Name,"/CNV_Heatmap_TGH_",Target_gene_name,
                      "_Top",HeatmapTopGene,".pdf"),
        width = 10,  height = 8
      )   
        
        # CNV_TGH_Top.df
        Heatmap(CNV_TGH_Top_Sum.lt[["CNV_THG_Top.df"]] , name = "Num", col = col_fun, 
                show_column_names = F,show_row_names = F, 
                cluster_columns = T, top_annotation = column_ha)-> Heatmap.p
        draw(Heatmap.p, column_title = paste0("TGH CNV heatmaps (TOP",HeatmapTopGene,")","\n",Note),
             column_title_gp = gpar(fontsize = 16)) #%>% print()
        rm(Heatmap.p)
        
        Heatmap(CNV_TGH_Top_Sum.lt[["CNV_THG_Top.df"]] , name = "Num", col = col_fun, 
                show_column_names = F,show_row_names = F, 
                cluster_columns = F, top_annotation = column_ha)-> Heatmap.p
        draw(Heatmap.p, column_title = paste0("TGH CNV heatmaps (TOP",HeatmapTopGene,")","\n",Note),
             column_title_gp = gpar(fontsize = 16)) #%>% print()
        rm(Heatmap.p)
        
        # CNV_TGH_Top_2D.df
        Heatmap(CNV_TGH_Top_Sum.lt[["CNV_THG_Top_2D.df"]] , name = "Num", col = col_fun, 
                show_column_names = F,show_row_names = F, 
                cluster_columns = T, top_annotation = column_ha)-> Heatmap.p
        draw(Heatmap.p, column_title = paste0("TGH CNV heatmaps (TOP",HeatmapTopGene,")","\n",Note),
             column_title_gp = gpar(fontsize = 16)) #%>% print()
        rm(Heatmap.p)
        
        Heatmap(CNV_TGH_Top_Sum.lt[["CNV_THG_Top_2D.df"]] , name = "Num", col = col_fun, 
                show_column_names = F,show_row_names = F, 
                cluster_columns = F, top_annotation = column_ha)-> Heatmap.p
        draw(Heatmap.p, column_title = paste0("TGH CNV heatmaps (TOP",HeatmapTopGene,")","\n",Note),
             column_title_gp = gpar(fontsize = 16)) #%>% print()
        rm(Heatmap.p)
        
        rm(column_ha)
        
      dev.off()
      
      rm(CNV_TGH.df) 
    
##### Pheno Type #####
      Pheno.df <- read.delim(PhenoFileName, 
                               header=T,  sep="\t") # row.names = 1,
      Pheno_KLC.df <- Pheno.df[,c("submitter_id.samples","sample_type.samples")]
      GeneExp_T.df <- t(GeneExp.df) # %>% as.data.frame()
      GeneExp_T.df <- data.frame(submitter_id.samples = rownames(GeneExp_T.df),GeneExp_T.df)
      GeneExp_T.df <- GeneExp_T.df[,c("submitter_id.samples",row.names(GeneExp.df)[grepl(Target_gene_name, row.names(GeneExp.df))])]
      Pheno_KLC.df <- dplyr::left_join(Pheno_KLC.df,GeneExp_T.df) %>% na.omit()
       
      # Bar  
      Bar.p <- ggplot(data = Pheno_KLC.df, aes(x = submitter_id.samples, y = ENSG00000131747.13, fill = sample_type.samples)) +
        geom_bar(stat = "identity")
      Bar.p
      
      
      ##### Density plot #####
      # https://www.omicsclass.com/article/1555
      library(plyr)
      mu <- ddply(Pheno_KLC.df, "sample_type.samples", summarise, grp.mean=mean(ENSG00000131747.13))
      head(mu)
      
      ggplot(Pheno_KLC.df, aes(x=ENSG00000131747.13, color=sample_type.samples)) +
        geom_density()

      # Add vline of the average
      TGeneDen.p <- ggplot(Pheno_KLC.df,aes(ENSG00000131747.13,fill=sample_type.samples, color=sample_type.samples)) + 
        xlab("Expression level") + 
        geom_density(alpha = 0.6, fill = "lightgray") + 
         geom_vline(data=mu, aes(xintercept=grp.mean, color=sample_type.samples),
                   linetype="dashed")
      TGeneDen.p %>% BeautifyggPlot(LegPos = c(0.85, 0.85),AxisTitleSize=1.7,
                                   OL_Thick = 1.5) + 
        labs(title= Target_gene_name,
             x ="Expression level", y = "Density")
      
      ## Count the sample number of different phenotype
      nrow(Pheno_KLC.df[Pheno_KLC.df[,"sample_type.samples"] %in% "Primary Tumor",])
      nrow(Pheno_KLC.df[Pheno_KLC.df[,"sample_type.samples"] %in% "Solid Tissue Normal",])
      nrow(Pheno_KLC.df[Pheno_KLC.df[,"sample_type.samples"] %in% "Metastatic",])
      
      