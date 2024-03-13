library(ggplot2)
library(RColorBrewer)
library(data.table)

df = read.csv('input/WPS/gene_category_summary.csv',row.names = 1)
df$undetermined = df$total - df$responsive - df$nonresponsive

# the first layer
data2 = data.frame(category = rownames(df), value = df$total)
data2 = data2[data2$category %in% c('zero','low','dynamic','high'),]
data2$layer=1
data2$value1=data2$value/sum(data2$value) # the position in the circle
# align the order of the table 
data2 = rbind(data2[data2$category == 'high',],
                data2[data2$category == 'dynamic',],
                data2[data2$category == 'low',],
                data2[data2$category == 'zero',])

# the second layer
tmp = as.data.table(df[,2:4])
tmp$category = rownames(df)
data1= as.data.frame(melt(tmp,id.vars = 'category'))
colnames(data1) = c('category','responsiveness','value')
data1$layer=2
data1 = rbind(data1[data1$category == 'high',],
              data1[data1$category == 'dynamic',],
              data1[data1$category == 'low',],
              data1[data1$category == 'zero',])
data1$value1=data1$value/sum(data1$value)
# make sure the row orders of data1 and data2 are aligned and layed out in the order of display


colors = c( 'undetermined' = 'white',
            'responsive' = '#D95319',
            'nonresponsive' = 'grey',
            'high' = '#77AC30',
            'dynamic' = '#EDB120',
            'low' = '#7E2F8E',
            'zero' = '#A2142F')

p=ggplot() +
  geom_col(aes(x = layer, y = value1, fill = category,group=layer),
           data =data2,width = 1, color = 'black')+
  geom_text(aes(label = category, x= layer, y = value1),
            data = data2,position=position_fill(vjust=0.5))+
  geom_col(aes(x = layer, y = value1, fill = responsiveness,group=layer), 
           data = data1,width = 1, color = 'black')+
  geom_text(aes(label = value, x= layer, y = value1),
            data = data1,position=position_fill(vjust=0.5))+
  scale_fill_manual(values = colors)+
  xlim(c(-0.25, 3))+
  coord_polar(theta = "y") +
  coord_polar(theta = "y") + theme_classic()

#

p

ggsave(file="figures/stacked_pie_gene_category_summary.pdf", plot=p,width=8, height=8)


# plot the DEG distribution for iCEL genes only 
conditionInfo = read.csv('input/WPS/RNAi_condition_information.csv')
conditionInfo = conditionInfo[conditionInfo$isICEL, ]
N_DE = conditionInfo$N_DE_targetExcluded
sum(N_DE == 0)
dev.off()
pdf('figures/N_DE_distribution_iCEL_only.pdf',width = 5,height = 5)
hist(log2(N_DE+0.1),xlab = 'log2(number of DEGs)',breaks = 50)
abline(v = log2(5+0.1))
dev.off()





