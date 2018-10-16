source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/plot_utils.R')
source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')

counts <- read.delim('Experimental_group_TB7Sc_batch_combine_counts.csv',sep=',',row.names = 1,header = T)
colnames(counts) <- gsub('^X','',colnames(counts))
group_table <- read.csv('group_table_demo_batch_TB7Sc_45_SOL_3PCA.csv',sep=',',as.is = T,row.names = 1)

cis_genes=read.table('TB7Sc_cis_genes.txt',header=F,sep="\t",as.is=T)$V1
plot_vol_edgeR(counts,group_table[1:7,],gene_keep = cis_genes)

