library(eulerr)
library(matrixStats)
library(stringr)


# load flux table

myTbl = read.csv('output/fluxTable.csv')

# some supplementary figures

tbl = myTbl[,c("OFD_bounded_exp_resp_simi",
               "OFD_bounded_exp_resp",
               "OFD_bounded_exp_simi",
               "OFD_bounded_exp_only",
               "OFD_bounded_no_data")]
colnames(tbl) = str_remove(colnames(tbl),'OFD_bounded_')
fit1 <- euler(tbl, 
              shape = "circle",
              #loss = 'region',
              #loss_aggregator = 'max', 
              control = list(extraopt = TRUE, extraopt_threshold = 0)
)
# double check if the fit of euler plot area is good
fit1
totalN = sum(rowAnys(as.matrix(myTbl[,c("OFD_bounded_exp_resp_simi","OFD_bounded_exp_resp","OFD_bounded_exp_simi","OFD_bounded_exp_only","OFD_bounded_no_data")])))
# total count should equal
totalN == sum(fit1$original.values)

pdf(paste('figures/euler_plot_bounded_OFD_rxns.pdf',sep = ''),width = 7,height = 7)
plot(fit1, quantities = TRUE, fill = c('#FFA07A','#D8BFD8','#F0E68C','#E0FFFF','#87CEEB'),
     main = paste('make sure total is ',totalN))
dev.off()

tbl = myTbl[,c("PFD_bounded_exp_resp_simi",
               "PFD_bounded_exp_resp",
               "PFD_bounded_exp_simi",
               "PFD_bounded_exp_only",
               "PFD_bounded_no_data")]
colnames(tbl) = str_remove(colnames(tbl),'PFD_bounded_')
fit2 <- euler(tbl, 
              shape = "circle",
              #loss = 'region',
              #loss_aggregator = 'max', 
              control = list(extraopt = TRUE, extraopt_threshold = 0)
)
fit2
totalN = sum(rowAnys(as.matrix(myTbl[,c("PFD_bounded_exp_resp_simi","PFD_bounded_exp_resp","PFD_bounded_exp_simi","PFD_bounded_exp_only","PFD_bounded_no_data")])))
# total count should equal
totalN == sum(fit2$original.values)

pdf(paste('figures/euler_plot_bounded_PFD_rxns.pdf',sep = ''),width = 7,height = 7)
plot(fit2, quantities = TRUE, fill = c('#FFA07A','#D8BFD8','#F0E68C','#E0FFFF','#87CEEB'),
     main = paste('make sure total is ',totalN))
dev.off()

