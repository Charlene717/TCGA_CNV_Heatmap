GeneExp_Group <- function(GeneExp.df, Target_gene_name, Mode="Mean",SD=1) { 
  
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
      GeneExp_TarGene.df <- GeneExp.df[Target_gene_name,] 
      GeneExp_TarGene.df <- GeneExp_TarGene.df[,order(GeneExp_TarGene.df[1,])]
      
      if(Mode_Group$Mode=="Mean"){
        if(Mode_Group$SD==0){
          GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD*(Mode_Group$SD)]
          GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] < Target_gene_Mean-Target_gene_SD*(Mode_Group$SD)]
        }else{
          GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD*(Mode_Group$SD)]
          GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] <= Target_gene_Mean-Target_gene_SD*(Mode_Group$SD)]
        }
        #rm(Target_gene_Mean, Target_gene_SD)
        
        ## GeneExp_medium.set
        # # Test
        # min(abs(GeneExp_TarGene.df[1,]-Target_gene_Mean))
        # which(GeneExp_TarGene.df==10.5429) 
        
        Len <- length(c(GeneExp_high.set,GeneExp_low.set))
        Mean <- min(abs(GeneExp_TarGene.df[1,]-Target_gene_Mean)) 
        Mean.loc <- which(GeneExp_TarGene.df[1,]-Target_gene_Mean == Mean)
        GeneExp_medium.set <- colnames(GeneExp_TarGene.df)[(Mean.loc-Len/2):(Mean.loc+Len/2)]
        rm(Len,Mean,Mean.loc)
        
      }else{
        if(Mode_Group$Q2=="Only"){ # Mode="Quartile"
          GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Q[3]]
          GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] < Target_gene_Q[3]]
        }else{
          GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Q[4]]
          GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] <= Target_gene_Q[2]]
        }
        #rm(Target_gene_Q)
        
        Len <- length(c(GeneExp_high.set,GeneExp_low.set))
        Q2 <- min(abs(GeneExp_TarGene.df[1,]-Target_gene_Q[2])) 
        Q2.loc <- which(GeneExp_TarGene.df[1,]-Target_gene_Q[2] == Q2)
        GeneExp_medium.set <- colnames(GeneExp_TarGene.df)[(Q2.loc-Len/2):(Q2.loc+Len/2)]
        rm(Len,Q2,Q2.loc)
      }
      
      # GeneExp_medium.set <- setdiff(colnames(GeneExp.df),
      #                               c(GeneExp_high.set,GeneExp_low.set))
      
      rm(GeneExp_TarGene.df)
    
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
    TGeneDen_SD.p <- ggPlot_vline(TGeneDen.p,data,                         
                                  Line1 = Target_gene_Mean+Target_gene_SD,
                                  Line2 = Target_gene_Mean,
                                  Line3 = Target_gene_Mean-Target_gene_SD)
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
  
  ## Set export
    GeneExp_Group.lt <- list()
    GeneExp_Group.lt[["GeneExp_high.set"]] <- GeneExp_high.set
    GeneExp_Group.lt[["GeneExp_low.set"]] <- GeneExp_low.set
    GeneExp_Group.lt[["GeneExp_medium.set"]] <- GeneExp_medium.set
    GeneExp_Group.lt[["Target_gene_Mean"]] <- Target_gene_Mean
    GeneExp_Group.lt[["Target_gene_SD"]] <- Target_gene_SD
    GeneExp_Group.lt[["Target_gene_Q"]] <- Target_gene_Q
    
  return(GeneExp_Group.lt)
}