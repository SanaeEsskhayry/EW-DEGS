DEGs<-function(countData, colData, design, FoldChange, condition){
  
  #the columns of the "data" and the rows of the "pheno" (information about samples) are in the same order.
  colData<-colData[colnames(countData),]
  
  #The deseq2 package require the count data values to be integers 
  #save the gene names in a variable
  
  genes=row.names(countData)
  
  #convert the data values to integers
  countData=apply(countData,2,as.integer) #(2=>column, 1 =>Row)
  
  
  #rename the rows of the data
  row.names(countData)=genes
  
  ################ DO the differential EXP analysis using DeSeq2 #########################
  
  #specify the conditions that you want to compare according to the phenotypic table
  
  
  conds<-names(table(colData[condition]))
  
  
  
  #create a deseq dataset object
  
  colData["Cell"] <- colData[condition]
  
  Cell <- as.formula(colData["Cell"])
  
  dds= DESeqDataSetFromMatrix( countData = countData , colData = colData, design = ~ Cell )
  
  #run the deseq2 worflow
  dds.run = DESeq(dds)
  
  #specifying teh contrast (to make a res object based on two specific conditions)
  res=results(dds.run)

  # remove nulls
  res=as.data.frame(res[complete.cases(res), ])
  
  #chose the statstical significant differentaily expressed genes (DEGs) based
  #on the p adjusted value less than 0.05 and biological significance  based
  #on the fold change more than 2
  deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)>as.numeric(FoldChange),]
  
  return(deseq.deg)
}