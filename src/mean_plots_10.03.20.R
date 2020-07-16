library(data.table)
library(stringr)
coldata <- fread('/Volumes/Seagate Backup /hd_in/09.19_hd/CategoricalAnnotation.txt', data.table = F)
rownames(coldata) <- coldata$Sample
coldata$Sample <- as.character(coldata$Sample)

protein <- fread('/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/20200309_HDProteomics_Log2_SubtractMedianColumn_AverageSample_+20.txt', data.table=F)
protein <- protein[,c('Gene', colnames(protein)[which(colnames(protein) %in% colnames(rna))])]

protein_fb <- protein[, c(colnames(protein)[which(colnames(protein) %in% subset(coldata, coldata$CellType == 'FB')$Sample)], 'Gene')]
colnames(protein_fb)[-ncol(protein_fb)] <- paste(colnames(protein_fb)[-ncol(protein_fb)], 'protein', sep='_')
colnames(protein_fb)[ncol(protein_fb)] <- 'gene_name'
protein_in <- protein[, c(colnames(protein)[which(colnames(protein) %in% subset(coldata, coldata$CellType == 'IN')$Sample)], 'Gene')]
colnames(protein_in)[-ncol(protein_in)] <- paste(colnames(protein_in)[-ncol(protein_in)], 'protein', sep='_')
colnames(protein_in)[ncol(protein_in)] <- 'gene_name'

rna <- fread('/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/1_readcounts/gene_count_matrix_2.csv', data.table=F)

gene_transcript <- fread('/Volumes/Seagate Backup /annotation/human/gencode/gencode.v30.annotation.transcr.gene.tab', data.table=F, header=F)
colnames(gene_transcript) <- c('transcript_id', 'gene_id', 'gene_name', 'gene_type') 

rna <- merge(rna, unique(gene_transcript[,c(2,3)]), by.x='Geneid', by.y='gene_id')
rownames(rna) <- rna$Geneid

rna <- subset(rna, rna$gene_name %in% protein$Gene)
rna_aggr <- aggregate(rna[,which(colnames(rna) %in% coldata$Sample)], by=list(rna$gene_name), FUN=sum)
rownames(rna_aggr) <- rna_aggr$Group.1
rna_aggr <- rna_aggr[,-1]

for (i in 1:ncol(rna_aggr)){
  rna_aggr[,i] <- ifelse(rna_aggr[,i] == 0, NaN, rna_aggr[,i])
}

rna_aggr_log2 <- log2(rna_aggr)
normalizeMedianValues(rna_aggr_log2)
rna_aggr_log2 <- rna_aggr_log2 + 20

rna_fb <- rna_aggr_log2[, c(colnames(rna_aggr_log2)[which(colnames(rna_aggr_log2) %in% subset(coldata, coldata$CellType == 'FB')$Sample)])]
colnames(rna_fb)[-ncol(rna_fb)] <- paste(colnames(rna_fb)[-ncol(rna_fb)], 'rna', sep='_')

rna_in <- rna_aggr_log2[, c(colnames(rna_aggr_log2)[which(colnames(rna_aggr_log2) %in% subset(coldata, coldata$CellType == 'IN')$Sample)])]
colnames(rna_in)[-ncol(rna_in)] <- paste(colnames(rna_in)[-ncol(rna_in)], 'rna', sep='_')

fb <- merge(protein_fb, rna_fb, by.x='gene_name', by.y='row.names')
rownames(fb) <- make.unique(fb$gene_name)
fb <- fb[,-1]
fb_hd <- fb[,c(colnames(fb)[startsWith(colnames(fb), 'HD')])]

fb_hd_plot <- data.frame(fb_hd_yaxis = rowMeans(fb_hd[,endsWith(colnames(fb_hd), 'protein')]),
                         fb_hd_xaxis = rowMeans(fb_hd[,endsWith(colnames(fb_hd), 'rna')]))

ggplotRegression <- function (fit, title, xlab, ylab) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste(title, paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                                    "Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)), sep='\n'),
         x=xlab, y=ylab) + theme_classic()
}
fb_hd_lm <- lm(data=fb_hd_plot, fb_hd_yaxis~fb_hd_xaxis)

fb_ctrl <- fb[,c(colnames(fb)[startsWith(colnames(fb), 'C')])]

fb_ctrl_plot <- data.frame(fb_ctrl_yaxis = rowMeans(fb_ctrl[,endsWith(colnames(fb_ctrl), 'protein')]),
                           fb_ctrl_xaxis = rowMeans(fb_ctrl[,endsWith(colnames(fb_ctrl), 'rna')]))

fb_ctrl_lm <- lm(data=fb_ctrl_plot, fb_ctrl_yaxis~fb_ctrl_xaxis)

iN <- merge(protein_in, rna_in, by.x='gene_name', by.y='row.names')
rownames(iN) <- make.unique(iN$gene_name)
iN <- iN[,-1]
in_hd <- iN[,c(colnames(iN)[startsWith(colnames(iN), 'HD')])]

in_hd_plot <- data.frame(in_hd_yaxis = rowMeans(in_hd[,endsWith(colnames(in_hd), 'protein')]),
                         in_hd_xaxis = rowMeans(in_hd[,endsWith(colnames(in_hd), 'rna')]))

in_hd_lm <- lm(data=in_hd_plot, in_hd_yaxis~in_hd_xaxis)

in_ctrl <- iN[,c(colnames(iN)[startsWith(colnames(iN), 'C')])]

in_ctrl_plot <- data.frame(in_ctrl_yaxis = rowMeans(in_ctrl[,endsWith(colnames(in_ctrl), 'protein')]),
                           in_ctrl_xaxis = rowMeans(in_ctrl[,endsWith(colnames(in_ctrl), 'rna')]))

in_ctrl_lm <- lm(data=in_ctrl_plot, in_ctrl_yaxis~in_ctrl_xaxis)

library(limma)
ggarrange(ggplotRegression(in_ctrl_lm, 'iN Control', xlab='normalized mean RNA expression', ylab='normalized mean protein expression'), 
          ggplotRegression(in_hd_lm, 'iN HD', xlab='normalized mean RNA expression', ylab='normalized mean protein expression'), 
          ggplotRegression(fb_ctrl_lm, 'Fibroblasts Control', xlab='normalized mean RNA expression', ylab='normalized mean protein expression'), 
          ggplotRegression(fb_hd_lm, 'Fibroblasts HD', xlab='normalized mean RNA expression', ylab='normalized mean protein expression'))


in_protein_plot <- data.frame(in_ctrl_yaxis = rowMeans(in_ctrl[,endsWith(colnames(in_ctrl), 'protein')]),
                           in_hd_xaxis = rowMeans(in_hd[,endsWith(colnames(in_hd), 'protein')]))
t.test(in_protein_plot)
in_protein_lm <- lm(data=in_protein_plot, in_ctrl_yaxis~in_hd_xaxis)

in_rna_plot <- data.frame(in_ctrl_yaxis = rowMeans(in_ctrl[,endsWith(colnames(in_ctrl), 'rna')]),
                              in_hd_xaxis = rowMeans(in_hd[,endsWith(colnames(in_hd), 'rna')]))

in_rna_lm <- lm(data=in_rna_plot, in_ctrl_yaxis~in_hd_xaxis)

fb_protein_plot <- data.frame(fb_ctrl_yaxis = rowMeans(fb_ctrl[,endsWith(colnames(fb_ctrl), 'protein')]),
                              fb_hd_xaxis = rowMeans(fb_hd[,endsWith(colnames(fb_hd), 'protein')]))

fb_protein_lm <- lm(data=fb_protein_plot, fb_ctrl_yaxis~fb_hd_xaxis)

fb_rna_plot <- data.frame(fb_ctrl_yaxis = rowMeans(fb_ctrl[,endsWith(colnames(fb_ctrl), 'rna')]),
                              fb_hd_xaxis = rowMeans(fb_hd[,endsWith(colnames(fb_hd), 'rna')]))

fb_rna_lm <- lm(data=fb_rna_plot, fb_ctrl_yaxis~fb_hd_xaxis)

ggarrange(ggplotRegression(in_protein_lm, 'iN Protein (control vs HD)', xlab='normalized mean control expression', ylab='normalized mean HD expression'), 
          ggplotRegression(in_rna_lm, 'iN RNA (control vs HD)', xlab='normalized mean control expression', ylab='normalized mean HD expression'), 
          ggplotRegression(fb_protein_lm, 'Fb Protein (control vs HD)', xlab='normalized mean control expression', ylab='normalized mean HD expression'), 
          ggplotRegression(fb_rna_lm, 'Fb RNA (control vs HD)', xlab='normalized mean control expression', ylab='normalized mean HD expression'))

