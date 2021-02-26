setwd('/Volumes/My Passport/hd_in/24.02.20/')

library(ggplot2)
library(RColorBrewer)
library(data.table)
library(DESeq2)
library(stringr)
library(rjson)
library(xlsx)
library(Hmisc)

more_10 <- function(row){
  return(length(which(row > 10)))
}

gene_transcript <- fread('/Volumes/My Passport/annotation/human/gencode/v30/gencode.v30.annotation.transcr.gene.tab', data.table=F, header=F)
colnames(gene_transcript) <- c('transcript_id', 'gene_id', 'gene_name', 'gene_type') 

coldata <- fread('/Volumes/My Passport/hd_in/09.19_hd/CategoricalAnnotation.txt', data.table = F)
rownames(coldata) <- coldata$Sample
coldata$Sample <- as.character(coldata$Sample)

rna <- fread('/Volumes/My Passport/hd_in/09.19_hd/3_stdmapping/1_readcounts/gene_count_matrix_2.csv', data.table=F)
rownames(rna) <- rna$Geneid

colnames(rna) <- ifelse(colnames(rna) == 'C52_IN', 'CBS_IN',
                        ifelse(colnames(rna) == 'C52_FB', 'CBS_FB', 
                               ifelse(colnames(rna) == 'C1_IN', 'CKP_IN', colnames(rna))))

# Filter out genes with less than 10 reads in one sample
rna_expressed <- rna[which(apply(rna[,coldata$Sample], 1, more_10) > 0), coldata$Sample]

dds_rna <- DESeqDataSetFromMatrix(rna_expressed[,coldata$Sample], coldata, design= ~Stage)
dds_rna$Stage <- relevel(dds_rna$Stage, 'CTRL')
dds_rna <- DESeq(dds_rna)
rna_expressed_norm <- as.data.frame(counts(dds_rna, normalize=T))
dds_rna$sizeFactor

# RNAseq ----
rna_expressed_norm$gene_id <- rownames(rna_expressed_norm)

rna_norm_test <- apply(rna_expressed_norm, 1, FUN=function(row){
  row <- as.data.frame(row)
  colnames(row) <- row['gene_id',]
  
  row <- row[coldata$Sample,,drop=F]
  row <- merge(row, coldata[,4,drop=F], by='row.names')
  row[,2] <- as.numeric(as.character(row[,2]))
  
  mean_in <- mean(subset(row, row$CellType == 'IN')[,2])
  mean_fb <- mean(subset(row, row$CellType == 'FB')[,2])
  
  gene_id <- colnames(row)[2]
  colnames(row) <- c('sample', 'expression', 'condition')
  row$expression <- as.numeric(as.character(row$expression))
  
  if(sum(row$expression) > 0){
    if(sum(subset(row, row$condition == 'IN')$expression) > 0){
      normality_in <- shapiro.test(subset(row, row$condition == 'IN')$expression)
    }
    else{
      normality_in <- data.frame(p.value=NA)
    }
    if(sum(subset(row, row$condition == 'FB')$expression) > 0){
      normality_fb <- shapiro.test(subset(row, row$condition == 'FB')$expression)
    }
    else{
      normality_fb <- data.frame(p.value=NA)
    }
    
    fb_expression <- subset(row, row$condition == 'FB')$expression
    in_expression <- subset(row, row$condition == 'IN')$expression
    test <- t.test(fb_expression, in_expression, paired=F, exact=F)
    
    results <- c(gene_id, mean_in, mean_fb, log2((mean_fb/mean_in)), normality_in$p.value, normality_fb$p.value, test$statistic, test$p.value, test$conf.int)
    
  }
  else{
    results <- c(gene_id, mean_in, mean_fb, log2((mean_fb/mean_in)), NA, NA, NA, NA, NA, NA, NA)
  }
  
  
  names(results) <- c("Gene id", "Mean iN", "Mean FB", "Log2FC",  'Shapiro Pvalue iN', 'Shapiro Pvalue FB', 'T', 'Pvalue', 'Low conf int', 'High conf int')
  
  return(results)
})

rna_norm_test <- as.data.frame(t(rna_norm_test))
rna_norm_test$`Gene id` <- as.character(rna_norm_test$`Gene id`)
rna_norm_test$`Mean iN` <- as.numeric(as.character(rna_norm_test$`Mean iN`))
rna_norm_test$`Mean FB` <- as.numeric(as.character(rna_norm_test$`Mean FB`))
rna_norm_test$`Log2FC` <- as.numeric(as.character(rna_norm_test$`Log2FC`))
rna_norm_test$`Shapiro Pvalue iN` <- as.numeric(as.character(rna_norm_test$`Shapiro Pvalue iN`))
rna_norm_test$`Shapiro Pvalue FB` <- as.numeric(as.character(rna_norm_test$`Shapiro Pvalue FB`))
rna_norm_test$T <- as.numeric(as.character(rna_norm_test$T))
rna_norm_test$Pvalue <- as.numeric(as.character(rna_norm_test$Pvalue))
rna_norm_test$`Low conf int` <- as.numeric(as.character(rna_norm_test$`Low conf int`))
rna_norm_test$`High conf int` <- as.numeric(as.character(rna_norm_test$`High conf int`))
# Filter out genes that are not normally distributed
rna_norm_test <- rna_norm_test[which(rna_norm_test$`Shapiro Pvalue iN` > 0.05 & rna_norm_test$`Shapiro Pvalue FB` > 0.05),]
# Get the counts of the ones that are normally distributed
rna_expressed_norm <- subset(rna_expressed_norm, rna_expressed_norm$gene_id %in% rna_norm_test$`Gene id`)

# Which one of those is significantly different?
rna_signdiff <- rna_norm_test[which(rna_norm_test$Pvalue < 0.05),]
# Tag them in one
rna_norm_test$type <- ifelse(rna_norm_test$Pvalue < 0.05 & rna_norm_test$Log2FC < 0, 'Downregulated', 
                             ifelse(rna_norm_test$Pvalue < 0.05 & rna_norm_test$Log2FC > 0, 'Upregulated', 'Not significant'))

rna_norm_test$colours <- ifelse(rna_norm_test$type == 'Downregulated', 'steelblue', 
                                ifelse(rna_norm_test$type == 'Upregulated', 'tomato3', 'black'))

pdf(paste(getwd(), '/plots/rna_fb.in_meanplot.pdf', sep=''))
plot(log2(rna_norm_test$`Mean iN`+0.5), 
     log2(rna_norm_test$`Mean FB`+0.5), 
     col=rna_norm_test$colours, 
     cex=0.5, 
     pch=16, 
     xlab='log2(mean iN)', 
     ylab='log2(mean FB)',
     main='RNA FB vs iN (p-value < 0.05; |log2FC| > 0)')

legend("bottomright", legend = c(paste("up (",as.numeric(table(rna_norm_test$type)["Upregulated"]),")",sep=""),
                                 paste("down (",as.numeric(table(rna_norm_test$type)["Downregulated"]),")",sep = ""),
                                 paste("not significant (",as.numeric(table(rna_norm_test$type)["Not significant"]),")",sep = "")),
       pch=16,col=c("firebrick3","steelblue4","black"),cex=1)

dev.off()

svg(paste(getwd(), '/plots/rna_fb.in_meanplot.svg', sep=''))
plot(log2(rna_norm_test$`Mean iN`+0.5), 
     log2(rna_norm_test$`Mean FB`+0.5), 
     col=rna_norm_test$colours, 
     cex=0.5, 
     pch=16, 
     xlab='log2(mean iN)', 
     ylab='log2(mean FB)',
     main='RNA FB vs iN (p-value < 0.05; |log2FC| > 0)')

legend("bottomright", legend = c(paste("up (",as.numeric(table(rna_norm_test$type)["Upregulated"]),")",sep=""),
                                 paste("down (",as.numeric(table(rna_norm_test$type)["Downregulated"]),")",sep = ""),
                                 paste("not significant (",as.numeric(table(rna_norm_test$type)["Not significant"]),")",sep = "")),
       pch=16,col=c("firebrick3","steelblue4","black"),cex=1)
dev.off()

# Significantly different
rna_norm_test_sign <- subset(rna_norm_test, rna_norm_test$type != 'Not significant')

# Up and down regulated for PANTHER
rna_norm_test_upreg <- merge(rna_norm_test[which(rna_norm_test$type == 'Upregulated'),], unique(gene_transcript[,c(2,3)]), by.x='row.names', by.y='gene_id')
write.table(rna_norm_test_upreg$gene_name, quote=F, col.names = F, row.names = F,
            file=paste(getwd(), '/3_stdmapping/GO_analysis/upregulated/rna_upreg_fb.in.tab', sep=''))

rna_norm_test_dwnreg <- merge(rna_norm_test[which(rna_norm_test$type == 'Downregulated'),], unique(gene_transcript[,c(2,3)]), by.x='row.names', by.y='gene_id')
write.table(rna_norm_test_dwnreg$gene_name, quote=F, col.names = F, row.names = F,
            file=paste(getwd(), '/3_stdmapping/GO_analysis/downregulated/rna_dwnreg_fb.in.tab', sep=''))

rna_norm_test <- merge(rna_norm_test, unique(gene_transcript[,c('gene_id', 'gene_name')]), by.x='Gene id', by.y='gene_id')
write.table(rna_norm_test$gene_name, quote=F, col.names = F, row.names = F,
            file=paste(getwd(), '/3_stdmapping/GO_analysis/rna_expressed.tab', sep=''))


to_file <- c("Gene", "50% quartile iN", "50% quartile FB", "Log2FC", "Shapiro Pvalue - iN",
             "Shapiro Pvalue - FB", "T", "Pvalue", "Low conf int", "High conf int", "type")
rna_norm_test_upreg_toxlsx <- rna_norm_test_upreg
colnames(rna_norm_test_upreg_toxlsx) <- c("gene_id","Gene id","50% quartile iN","50% quartile FB",
                                   "Log2FC","Shapiro Pvalue - iN",
                                   "Shapiro Pvalue - FB","T","Pvalue","Low conf int","High conf int","type",
                                   "colours","Gene")
rna_norm_in_test_dwnreg_toxlsx <- rna_norm_test_dwnreg
colnames(rna_norm_in_test_dwnreg_toxlsx) <- c("gene_id","Gene id","50% quartile iN","50% quartile FB",
                                       "Log2FC","Shapiro Pvalue - iN",
                                       "Shapiro Pvalue - FB","T","Pvalue","Low conf int","High conf int","type",
                                       "colours","Gene")

write.xlsx(x = rna_norm_test_upreg_toxlsx[,to_file], row.names = F,
           file = paste(getwd(), "/tables/rna_fb.in_signdiff_genes.xlsx", sep=''),
           sheetName = "Upreg in FB | Downreg in iN")
write.xlsx(x = rna_norm_in_test_dwnreg_toxlsx[,to_file],  row.names = F,
           file = paste(getwd(), "/tables/rna_fb.in_signdiff_genes.xlsx", sep=''),
           sheetName = "Downreg in FB | Upreg in iN", append=TRUE)

save.image('RNA_in_vs_fb.RData')
