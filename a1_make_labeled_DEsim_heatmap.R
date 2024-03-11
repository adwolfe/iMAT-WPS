library(pheatmap)
library(stringr)


# reproduce the heatmap
input = read.csv('./../../MetabolicLibrary/2_DE/output/cosineSimilarity_FC_denoised_stationery_metabolic.csv',row.names = 1)
N_DETbl = read.csv('./../../MetabolicLibrary/2_DE/output/RNAi_condition_metaInfo.csv',row.names = 1)
WBID = read.csv('./../../MetabolicLibrary/input_data/otherTbls/WBIDtbl.txt',sep = '\t')
ann2 = data.frame(conID = colnames(input))
ann2$log2NDE = log2(N_DETbl$N_DE_targetExcluded[match(ann2$conID, str_replace(N_DETbl$RNAiID,' ','_'))])
rownames(ann2) = ann2$conID
wormPath = read.csv('./../../MetabolicLibrary/input_data/otherTbls/WormPaths_Tables/wormPathTable.csv')
wormpathLevels = read.csv('./../../MetabolicLibrary/input_data/otherTbls/WormPaths_Tables/LEVEL1_GENES.csv',header = F)
ann2$RNAiName = N_DETbl$RNAi_geneName[match(ann2$conID, str_replace(N_DETbl$RNAiID,' ','_'))]
ann2$WBID = N_DETbl$RNAi_WBID[match(ann2$conID, str_replace(N_DETbl$RNAiID,' ','_'))]
ann2$wormPath = wormPath$LEVEL.1[match(ann2$WBID, wormPath$WormBase.ID)]
mylevels = wormpathLevels$V1
ann_col = list()
for (i in 1:length(mylevels)){
  if (any(str_detect(ann2$wormPath, paste(mylevels[i],';',sep = '')),na.rm = T) | 
      any(str_detect(ann2$wormPath, paste(mylevels[i],'$',sep = '')),na.rm = T)){
    ann2[,mylevels[i]] = 'No'
    ann2[which(str_detect(ann2$wormPath, paste(mylevels[i],';',sep = ''))),mylevels[i]] = 'Yes'
    ann2[which(str_detect(ann2$wormPath, paste(mylevels[i],'$',sep = ''))),mylevels[i]] = 'Yes'
    cols = c('red','white')
    names(cols) = c('Yes','No')     
    ann_col[[mylevels[i]]] = cols
  }
  
}



# we directly use the cosine distance (1-cosine similarity) to cluster
pdf('figures/DE_similarity_used_pair_heatmap.pdf',width = 10,height = 12)
upperl = 0.5
seq1 = seq(-upperl,upperl,length.out = 100)
cl_cosine = fastcluster::hclust(as.dist(1-as.matrix(input)), method = 'complete') # for reproducibility; the original hclust can produce the same tree but with same-height branches reordered!
color_palette <- rev(colorRampPalette(c("red", "grey100", "blue"))(100))
pheatmap(input, color = color_palette,breaks = seq1
         ,show_rownames = F,show_colnames = F
         ,annotation_col = ann2[,-which(colnames(ann2) %in% c('conID','WBID','RNAiName','wormPath')),drop = F]
         ,annotation_colors = ann_col
         ,fontsize = 7,annotation_legend = F
         ,cellheight = 1,cellwidth = 1
         ,cluster_rows = cl_cosine,cluster_cols = cl_cosine
)


# produce the labeling heatmap 
df = read.csv('input/WPS/DEsim_table_integration_summary.csv')
df$pair_gene1 = str_replace_all(df$pair_gene1,'\\.','_')
df$pair_gene2 = str_replace_all(df$pair_gene2,'\\.','_')
mat = matrix(data = 0, nrow = nrow(input), ncol = ncol(input))
mat = as.data.frame(mat)
colnames(mat) = str_replace_all(colnames(input),'\\.','_')
rownames(mat) = str_replace_all(rownames(input),'\\.','_')
for (i in 1:nrow(df)){
  mat[df$pair_gene1[i], df$pair_gene2[i]] = 1
}
pheatmap(mat, color = color_palette,breaks = seq1
         ,show_rownames = F,show_colnames = F
         ,fontsize = 7,annotation_legend = F
         ,cellheight = 1,cellwidth = 1
         ,cluster_rows = cl_cosine,cluster_cols = cl_cosine
)

dev.off()


# show ahcy1 and cbs-1
pdf('figures/DE_similarity_example_pair_heatmap.pdf',width = 10,height = 12)
upperl = 0.5
seq1 = seq(-upperl,upperl,length.out = 100)
cl_cosine = fastcluster::hclust(as.dist(1-as.matrix(input)), method = 'complete') # for reproducibility; the original hclust can produce the same tree but with same-height branches reordered!
color_palette <- rev(colorRampPalette(c("red", "grey100", "blue"))(100))
pheatmap(input, color = color_palette,breaks = seq1
         ,show_rownames = F,show_colnames = F
         ,annotation_col = ann2[,-which(colnames(ann2) %in% c('conID','WBID','RNAiName','wormPath')),drop = F]
         ,annotation_colors = ann_col
         ,fontsize = 7,annotation_legend = F
         ,cellheight = 1,cellwidth = 1
         ,cluster_rows = cl_cosine,cluster_cols = cl_cosine
)


# produce the labeling heatmap 
df = read.csv('input/WPS/DEsim_table_integration_summary.csv')
df = df[df$mets == 'hcys-L[c]',]
df$pair_gene1 = str_replace_all(df$pair_gene1,'\\.','_')
df$pair_gene2 = str_replace_all(df$pair_gene2,'\\.','_')
mat = matrix(data = 0, nrow = nrow(input), ncol = ncol(input))
mat = as.data.frame(mat)
colnames(mat) = str_replace_all(colnames(input),'\\.','_')
rownames(mat) = str_replace_all(rownames(input),'\\.','_')
for (i in 1:nrow(df)){
  mat[df$pair_gene1[i], df$pair_gene2[i]] = 1
}
pheatmap(mat, color = color_palette,breaks = seq1
         ,show_rownames = F,show_colnames = F
         ,fontsize = 7,annotation_legend = F
         ,cellheight = 1,cellwidth = 1
         ,cluster_rows = cl_cosine,cluster_cols = cl_cosine
)

dev.off()
