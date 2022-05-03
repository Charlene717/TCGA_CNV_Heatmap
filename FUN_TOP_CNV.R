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
  
  CNV_Top.df <- tibble::rownames_to_column(CNV_Sum.df) %>% arrange(desc(Sum)) 
  row.names(CNV_Top.df) <- CNV_Top.df[,1]
  CNV_Top.df <- CNV_Top.df[,-1]
  
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