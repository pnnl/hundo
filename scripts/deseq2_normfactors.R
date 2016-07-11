UqNormDESeq <- function(data, q){
   # data is a p x n data.frame, where n is the number of samples and p is the number of OTUs
   # q is the desired quantile to normalize with (between 0 and 1)#
   #

   
   col.q = apply(data, 2, function(x) quantile(x[x != 0],probs = q))
   g.q = quantile(data[data != 0], probs = q)
   sizes = col.q/g.q
   return(sizes)
 }


 scalar_vec <- UqNormDESeq(raw_data, 0.75)
 
 # Store the input values, intermediate calculations and results of an analysis of differential expression
  DESeq_data_set_object <- DESeqDataSetFromMatrix(countData = raw_count_data, colData = colData, design = ~Group)
  
  #Include normalization factors (sizeFactors)
  sizeFactors(DESeq_data_set_object) = scalar_vec
