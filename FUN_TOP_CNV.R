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