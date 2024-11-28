####### pathway analysis of solution space ######
library(stringr)
library(matrixStats)
library(reshape2)
library(pheatmap)

# load the directionally constrained reactions
fluxTbl = read.csv('./../output/fluxTable.csv')
dirConsRxns = fluxTbl$rxns[fluxTbl$OFD_bounded_exp_resp_simi == 1]

# # convert reaction sets to gene sets so that we can run enrichment analysis 
# rxnGeneMat = read.csv('./../../../MetabolicLibrary/4_networkAnalysis/input/rxnGeneMat_iCEL1314.csv',row.names = 1)
# # get relevant genes 
# dirConsGenes = colnames(rxnGeneMat)[colAnys(rxnGeneMat[rownames(rxnGeneMat) %in% dirConsRxns,]==1)]
# # these genes are associated with 2000 reactions, indicating the specificity is lost in the conversion. Give up this enrichment method 
# a = rownames(rxnGeneMat)[rowAnys(rxnGeneMat[,colnames(rxnGeneMat) %in% dirConsGenes]==1)]

# since we dont have reaction-centered WormPath, let's change to conventional subsystems in the model for this analysis 
model = readxl::read_excel(path = './../input/model/makeWormModel/iCEL1314_with_uptakes.xlsx')
colnames(model)[1] = 'rxnID'
model$PATHWAY[model$PATHWAY == 'NA'] = 'Uptake reactions'

# we can check the overrepresentation among either all reactions, or only internal reactions (then constrained reactions also internal only)
# play with this by define diff universe 
library(clusterProfiler)
universe = model$ID # or just internal

# sig enrich threshold
p_sig = 0.01

# model pathway table
TERM2GENE=data.frame(term = model$PATHWAY, gene = model$rxnID)
# keep everything in the universe
myrxns = intersect(dirConsRxns, universe)
TERM2GENE = TERM2GENE[TERM2GENE$gene %in% universe, ]

# first is enriched pathways for constrained reactions
constrained_enriched = enricher(myrxns, TERM2GENE=TERM2GENE, universe = universe, qvalueCutoff = 1, pvalueCutoff = p_sig,
             minGSSize = 1, maxGSSize = 20000)
constrained_enriched = as.data.frame(constrained_enriched)
constrained_enriched = constrained_enriched[order(constrained_enriched$p.adjust,decreasing = F),]
# add a few annotations
constrained_enriched$set_size = as.numeric(unlist(lapply(strsplit(constrained_enriched$BgRatio,'/'),function(x){x[1]})))
constrained_enriched$set_recall = constrained_enriched$Count / constrained_enriched$set_size
# the enriched reactions covers around half of the total reactions bounded.
sum(constrained_enriched$Count) / length(myrxns)


# then is enriched pathways for unconstrained reactions
unconstrained_enriched = enricher(setdiff(universe, myrxns), TERM2GENE=TERM2GENE, universe = universe, qvalueCutoff = 1, pvalueCutoff = p_sig,
                                minGSSize = 1, maxGSSize = 20000)
unconstrained_enriched = as.data.frame(unconstrained_enriched)
unconstrained_enriched = unconstrained_enriched[order(unconstrained_enriched$p.adjust,decreasing = F),]
# add a few annotations
unconstrained_enriched$set_size = as.numeric(unlist(lapply(strsplit(unconstrained_enriched$BgRatio,'/'),function(x){x[1]})))
unconstrained_enriched$set_recall = unconstrained_enriched$Count / unconstrained_enriched$set_size
# the enriched reactions covers around half of the total reactions bounded.
sum(unconstrained_enriched$Count) / length(setdiff(universe, myrxns))


# enrichment looks missing the many pathways that are not enriched for either - lets do the barplot with enrichment annotated by color 
# barplot of constrained reaction proportion of each subsystems 
plot_pathway_distribution = function(model, dirConsRxns,constrained_enriched, unconstrained_enriched){
  library('reshape')
  library('ggsci')
  library(data.table)

  catTbl1 = data.frame(row.names = unique(model$PATHWAY))
  catTbl1$constrained = 0
  catTbl1$unconstrained = 0

  for (i in 1:nrow(catTbl1)){
    rxns = model$ID[model$PATHWAY == rownames(catTbl1)[i]]
    catTbl1$constrained[i] = sum(rxns %in% dirConsRxns)
    catTbl1$unconstrained[i] = sum(!(rxns %in% dirConsRxns))
  }
  catTbl1 = as.data.frame(t(catTbl1))
  catTbl1$catName = rownames(catTbl1)
  catTbl1 = catTbl1[,c(ncol(catTbl1),1:(ncol(catTbl1)-1))]
  catTbl1_mlt = melt(as.data.table(catTbl1),id.vars = 'catName')
  catTbl1_mlt$catName = factor(catTbl1_mlt$catName,levels = c('unconstrained','constrained'))
  # rearrange the levels
  N_rxn = c()
  N_rxn_abs = c()
  for (i in 2:ncol(catTbl1)){
    rxns = model$ID[model$PATHWAY == colnames(catTbl1)[i]]
    N_rxn = c(N_rxn, sum(rxns %in% dirConsRxns)/length(rxns))
    N_rxn_abs = c(N_rxn_abs, sum(rxns %in% dirConsRxns))
  }
  names(N_rxn) = colnames(catTbl1)[2:ncol(catTbl1)]
  names(N_rxn_abs) = colnames(catTbl1)[2:ncol(catTbl1)]
  
  # first rank by absolute
  b = colnames(catTbl1)[2:ncol(catTbl1)]
  b[rank(N_rxn_abs,ties.method ='first')] = b
  N_rxn_abs = N_rxn_abs[b]
  N_rxn = N_rxn[b]
  # then rank by relative
  b[rank(-N_rxn,ties.method ='first')] = b
  
  # reorder
  N_rxn_abs = N_rxn_abs[b]
  catTbl1_mlt$variable = factor(catTbl1_mlt$variable,levels = b)

  N_rxn_abs_acc = c()
  for(i in 1:length(N_rxn_abs)){N_rxn_abs_acc = c(N_rxn_abs_acc, sum(N_rxn_abs[1:i]))}
  # plot
  myFill = c('orange','grey')
  names(myFill) = c('constrained','unconstrained')
  
  # color the x name based on sig 
  red_names <- constrained_enriched$ID
  blue_names <- unconstrained_enriched$ID
  black_names <- setdiff(levels(catTbl1_mlt$variable), c(red_names, blue_names))
  x_labels <- levels(catTbl1_mlt$variable)
  x_label_colors <- ifelse(x_labels %in% red_names, "red", ifelse(x_labels %in% blue_names, "blue", "black"))
  
  p1 <- ggplot(catTbl1_mlt, aes(x=variable, y = value,fill=catName)) +
    geom_bar(stat="identity", position="fill") + scale_y_continuous(labels = scales::percent,expand = c(0,0))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1, colour = x_label_colors),
          legend.title = element_blank(),
          legend.text = element_text(colour ='black'))+ 
    #scale_fill_discrete(myFill)+
    scale_fill_manual(values = myFill)+ #, labels = names(myFill)
    labs(y = 'Proportion')+
    geom_text(aes(label = value), position = position_fill(vjust = 0.5))
    
    
  p1 <- p1 + geom_line(data = data.frame(variable = names(N_rxn_abs), N_rxn_abs_acc = N_rxn_abs_acc / max(N_rxn_abs_acc)),
                       aes(x = variable, y = N_rxn_abs_acc, group = 1, fill = NULL), color = "black", linetype = "dashed") +
    scale_y_continuous(
      sec.axis = sec_axis(~ . * max(N_rxn_abs_acc), name = "Cumulative number of directionally constrained reactions"),
      labels = scales::percent, expand = c(0, 0)
    )  
  return(p1)
  #p1+p2
}

p = plot_pathway_distribution(model, dirConsRxns,constrained_enriched, unconstrained_enriched)


pdf(paste('figures/directionally_constrained_reactions_pathway_distribution.pdf',sep = ''),width = 15,height = 8)
p
dev.off()



####### annotate the difference in solution space ######
# per request of the reviewer, we will label reactions whose flux is significantly differentially predicted based on the perturbation constraints
fluxTbl = read.csv('./../output/fluxTable.csv')
# we arbitrarily define the differential reactions as those directionally constrained in the full integration (PFD level, highest confidence)
# but shows very different flux  in dual integration
diff_reactions = fluxTbl$rxns[
  fluxTbl$PFD_bounded_exp_resp_simi &
                    abs(fluxTbl$normalized_OFD_exp_only - fluxTbl$normalized_OFD_exp_resp_simi)/abs(fluxTbl$normalized_OFD_exp_resp_simi) > 0.8
  ]
# we will label these reaction with a star 
fluxTbl$rxns[fluxTbl$rxns %in% diff_reactions] = paste(fluxTbl$rxns[fluxTbl$rxns %in% diff_reactions],'*',sep = '')
write.csv(fluxTbl,'./output/fluxTable_with_star.csv',row.names = F)


####### annotate constraining of solution space for diff. types of rxns ######
library(stringr)
library(ggplot2)

fluxTbl = read.csv('./../output/fluxTable.csv')

# OFD
OFD_N_total = sum(fluxTbl$OFD_bounded_exp_resp_simi)/nrow(fluxTbl)
OFD_N_R = sum(fluxTbl$OFD_bounded_exp_resp_simi & str_detect(fluxTbl$rxns,'^R|DGR|BIO'))/sum(str_detect(fluxTbl$rxns,'^R|DGR|BIO'))
OFD_N_TP = sum(fluxTbl$OFD_bounded_exp_resp_simi & str_detect(fluxTbl$rxns,'^TC'))/sum(str_detect(fluxTbl$rxns,'^TC'))
OFD_N_EX = sum(fluxTbl$OFD_bounded_exp_resp_simi & str_detect(fluxTbl$rxns,'^DM|SNK|EX|UP'))/sum(str_detect(fluxTbl$rxns,'^DM|SNK|EX|UP'))
df = data.frame(type = c('all','internal','transport','exchange'),
                count = c(OFD_N_total,OFD_N_R,OFD_N_TP,OFD_N_EX),
                text = c(paste(sum(fluxTbl$OFD_bounded_exp_resp_simi),'/',nrow(fluxTbl),sep = ''),
                         paste(sum(fluxTbl$OFD_bounded_exp_resp_simi & str_detect(fluxTbl$rxns,'^R|DGR|BIO')),'/',sum(str_detect(fluxTbl$rxns,'^R|DGR|BIO')),sep = ''),
                         paste(sum(fluxTbl$OFD_bounded_exp_resp_simi & str_detect(fluxTbl$rxns,'^TC')),'/',sum(str_detect(fluxTbl$rxns,'^TC')),sep = ''),
                         paste(sum(fluxTbl$OFD_bounded_exp_resp_simi & str_detect(fluxTbl$rxns,'^DM|SNK|EX|UP')),'/',sum(str_detect(fluxTbl$rxns,'^DM|SNK|EX|UP')),sep = '')))
df_merged = df
df_merged$model = 'OFM'
# PFD
PFD_N_total = sum(fluxTbl$PFD_bounded_exp_resp_simi)/nrow(fluxTbl)
PFD_N_R = sum(fluxTbl$PFD_bounded_exp_resp_simi & str_detect(fluxTbl$rxns,'^R|DGR|BIO'))/sum(str_detect(fluxTbl$rxns,'^R|DGR|BIO'))
PFD_N_TP = sum(fluxTbl$PFD_bounded_exp_resp_simi & str_detect(fluxTbl$rxns,'^TC'))/sum(str_detect(fluxTbl$rxns,'^TC'))
PFD_N_EX = sum(fluxTbl$PFD_bounded_exp_resp_simi & str_detect(fluxTbl$rxns,'^DM|SNK|EX|UP'))/sum(str_detect(fluxTbl$rxns,'^DM|SNK|EX|UP'))
df = data.frame(type = c('all','internal','transport','exchange'),
                count = c(PFD_N_total,PFD_N_R,PFD_N_TP,PFD_N_EX),
                text = c(paste(sum(fluxTbl$PFD_bounded_exp_resp_simi),'/',nrow(fluxTbl),sep = ''),
                         paste(sum(fluxTbl$PFD_bounded_exp_resp_simi & str_detect(fluxTbl$rxns,'^R|DGR|BIO')),'/',sum(str_detect(fluxTbl$rxns,'^R|DGR|BIO')),sep = ''),
                         paste(sum(fluxTbl$PFD_bounded_exp_resp_simi & str_detect(fluxTbl$rxns,'^TC')),'/',sum(str_detect(fluxTbl$rxns,'^TC')),sep = ''),
                         paste(sum(fluxTbl$PFD_bounded_exp_resp_simi & str_detect(fluxTbl$rxns,'^DM|SNK|EX|UP')),'/',sum(str_detect(fluxTbl$rxns,'^DM|SNK|EX|UP')),sep = '')))
df$model = 'PFM'
df_merged = rbind(df_merged, df)

df_merged$type = factor(df_merged$type, levels = c('all','internal','transport','exchange'))
barplot = ggplot(df_merged, aes(x = type, y = count, fill = model)) +
  geom_bar(stat = "identity", position = 'dodge', color = "black", width = 0.8) +
  geom_text(aes(label = text), position = position_dodge(width  = 1), size = 4, color = "black") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(size = 14, color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5)
  ) +
  labs(
    x = "",  # Replace with your x-axis label
    y = "Coverage (%)",      # Replace with your y-axis label
  )

pdf('figures/proportion_directionally_constrained_reactions_by_reaction_type.pdf',width = 3.5,height = 4)
print(barplot)
dev.off()

