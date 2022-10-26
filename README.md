# ssgea-single-sample-geneset-enrichment-analysis

Single Sample Geneset Enrichment Analysis, extension of GSEA method, working at the level of a single sample rather than a sample population as in the original GSEA application.
ssGSEA reflects the degree of input gene up- or downregulated within a sample.

**calling library**

library(matrixStats)

library(circlize)

**read the data**

data1= readRDS("log.RDS")
gene_set= read.csv("C:/Users/shiva/Desktop/shivani college/rgitbt/cancer genomics/genes/markers2Sep.csv")
View(gene_set)
gene_sets= as.list(as.data.frame(gene_set))

# create function

ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  **Calculate enrichment score (es) for each sample (column)**
  
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    **Calc es for each gene set**
    
    es_sample = sapply(gene_sets, function(gene_set_idx) {
    
      **pos: match (within the gene set)**
      
      **neg: non-match (outside the gene set)**
      
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      **Normalize by gene number**
      
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
       **Use ssGSEA or not**
       
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  **Normalize by absolute diff between max and min**
  
  if (norm) es = es / diff(range(es))
  
   **Prepare output**
   
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}


res= ssgsea(data1,gene_sets,norm= F,scale=T)

View(res) #value of gene in that particular sample

# Transpose data
res1 = t(res)

head(res1)

mat= as.matrix(res1)

for(i in 1:nrow(mat)){
  vec= as.numeric(mat[i,])
  mat[i, 1:ncol(mat)] = (vec-mean(vec))/sd(vec)
}

View(mat)

# calling library for heatmap

library(ComplexHeatmap)

library(circlize)

# HEATMAP

Heatmap(t(mat),col= colorRamp2(c(-2,0,2),c("green","white","yellow")))  

![image](https://user-images.githubusercontent.com/66779651/197950042-106ceb9b-3391-45d3-b271-90ebe6e5e0b8.png)

            
