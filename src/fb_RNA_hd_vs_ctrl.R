setwd('/Volumes/My Passport/hd_in/')
load('RNA_in_vs_fb.RData')
# RNAseq FBs ----
coldata_fb <- subset(coldata, coldata$CellType == 'FB')

# Filter out genes with less than 10 reads in one sample
rna_expressed_fb <- rna[which(apply(rna[,coldata_fb$Sample], 1, more_10) > 0), coldata_fb$Sample]

# DESeq DEA and normalizing with median of ratios
dds_fb <- DESeqDataSetFromMatrix(rna_expressed_fb[,coldata_fb$Sample], coldata_fb, design= ~Stage)
dds_fb$Stage <- relevel(dds_fb$Stage, 'CTRL')
dds_fb <- DESeq(dds_fb)
res_fb <- results(dds_fb)

rna_norm_fb <- as.data.frame(counts(dds_fb, normalize=T))
rna_norm_fb$gene_id <- rownames(rna_norm_fb)
# rna_log2norm_fb <- cbind(log2(rna_norm_fb[,-ncol(rna_norm_fb)]+0.5), gene_id=rna_norm_fb[,ncol(rna_norm_fb)])
# For each gene (normalized and filtered by expression) 
# 1) Shapiro test
# 2) T test
rna_norm_fb_test <- apply(rna_norm_fb, 1, FUN=function(row){
  row <- as.data.frame(row)
  colnames(row) <- row['gene_id',]
  
  row <- row[coldata_fb$Sample,,drop=F]
  row <- merge(row, coldata_fb[,2,drop=F], by='row.names')
  row[,2] <- as.numeric(as.character(row[,2]))
  
  mean_ctrl <- mean(subset(row, row$Stage == 'CTRL')[,2])
  mean_hd <- mean(subset(row, row$Stage == 'HD')[,2])
  
  gene_id <- colnames(row)[2]
  colnames(row) <- c('sample', 'expression', 'condition')
  row$expression <- as.numeric(as.character(row$expression))
  # print(row)
  if(sum(row$expression) > 0){
    if(sum(subset(row, row$condition == 'CTRL')$expression) > 0){
      normality_ctrl <- shapiro.test(subset(row, row$condition == 'CTRL')$expression)
    }
    else{
      normality_ctrl <- data.frame(p.value=NA)
    }
    if(sum(subset(row, row$condition == 'HD')$expression) > 0){
      normality_hd <- shapiro.test(subset(row, row$condition == 'HD')$expression)
    }
    else{
      normality_hd <- data.frame(p.value=NA)
    }
    # normality <- shapiro.test(row$expression)
    test <- t.test(expression~condition, data=row, paired=F, exact=F)
    
    results <- c(gene_id, mean_ctrl, mean_hd, log2((mean_hd/mean_ctrl)), normality_ctrl$p.value, normality_hd$p.value, test$statistic, test$p.value, test$conf.int)
    
  }
  else{
    results <- c(gene_id, mean_ctrl, mean_hd, log2((mean_hd/mean_ctrl)), NA, NA, NA,NA, NA, NA)
  }
  
  names(results) <- c("Gene id", "Mean Control", "Mean HD", "Log2FC",  'Shapiro Pvalue - Control', 'Shapiro Pvalue - HD', 'T', 'Pvalue', 'Low conf int', 'High conf int')
  # names(results) <- c("Gene id", "Mean Control", "Mean HD", "Log2FC", 'T', 'Pvalue')
  
  return(results)
})

rna_norm_fb_test <- as.data.frame(t(rna_norm_fb_test))
rna_norm_fb_test$`Gene id` <- as.character(rna_norm_fb_test$`Gene id`)
rna_norm_fb_test$`Mean Control` <- as.numeric(as.character(rna_norm_fb_test$`Mean Control`))
rna_norm_fb_test$`Mean HD` <- as.numeric(as.character(rna_norm_fb_test$`Mean HD`))
rna_norm_fb_test$`Log2FC` <- as.numeric(as.character(rna_norm_fb_test$`Log2FC`))
rna_norm_fb_test$`Shapiro Pvalue - HD` <- as.numeric(as.character(rna_norm_fb_test$`Shapiro Pvalue - HD`))
rna_norm_fb_test$`Shapiro Pvalue - Control` <- as.numeric(as.character(rna_norm_fb_test$`Shapiro Pvalue - Control`))
rna_norm_fb_test$T <- as.numeric(as.character(rna_norm_fb_test$T))
rna_norm_fb_test$Pvalue <- as.numeric(as.character(rna_norm_fb_test$Pvalue))
rna_norm_fb_test$`Low conf int` <- as.numeric(as.character(rna_norm_fb_test$`Low conf int`))
rna_norm_fb_test$`High conf int` <- as.numeric(as.character(rna_norm_fb_test$`High conf int`))
# Filter out genes that are not normally distributed
# rna_norm_fb_test <- rna_norm_fb_test[which(rna_norm_fb_test$`Shapiro Pvalue - HD` > 0.05 & rna_norm_fb_test$`Shapiro Pvalue - Control` > 0.05),]
rna_norm_fb_test <- rna_norm_fb_test[which(rna_norm_fb_test$`Shapiro Pvalue - HD` > 0.05 & rna_norm_fb_test$`Shapiro Pvalue - Control` > 0.05),]
# Get the counts of the ones that are normally distributed
rna_norm_fb <- subset(rna_norm_fb, rna_norm_fb$gene_id %in% rna_norm_fb_test$`Gene id`)

# Which one of those is significantly different?
rna_signdiff_fb <- rna_norm_fb_test[which(rna_norm_fb_test$Pvalue < 0.05),]

# Tag them in one
rna_norm_fb_test$type <- ifelse(rna_norm_fb_test$Pvalue < 0.05 & rna_norm_fb_test$Log2FC < 0, 'Downregulated', 
                                ifelse(rna_norm_fb_test$Pvalue < 0.05 & rna_norm_fb_test$Log2FC > 0, 'Upregulated', 'Not significant'))

table(rna_norm_fb_test$type)
# Put them colorrsss
rna_norm_fb_test$colours <- ifelse(rna_norm_fb_test$type == 'Downregulated', 'steelblue', 
                                   ifelse(rna_norm_fb_test$type == 'Upregulated', 'tomato3', 'black'))

# Size of the points for the mean plot
rna_norm_fb_test$cexs <- ifelse(rownames(rna_norm_fb_test) %in% rownames(rna_signdiff_fb) , 1, 0.5)

# png('/Volumes/My Passport/hd_in/07.20_hd/plots/rna_fb_meanplot.png', width = 15, height = 15, units = "cm", res=200)
plot(log2(rna_norm_fb_test$`Mean Control`+0.5), 
     log2(rna_norm_fb_test$`Mean HD`+0.5), 
     col=rna_norm_fb_test$colours, 
     cex=rna_norm_fb_test$cexs, 
     pch=16, 
     xlab='log2(mean Control)', 
     ylab='log2(mean HD)',
     main='RNA FB - HD vs Control (p-value < 0.05; |log2FC| > 0)')

legend("bottomright", legend = c(paste("up (",as.numeric(table(rna_norm_fb_test$type)["Upregulated"]),")",sep=""),
                                 paste("down (",as.numeric(table(rna_norm_fb_test$type)["Downregulated"]),")",sep = ""),
                                 paste("not significant (",as.numeric(table(rna_norm_fb_test$type)["Not significant"]),")",sep = "")),
       pch=16,col=c("firebrick3","steelblue4","black"),cex=1)

# dev.off()

# Significantly different 
rna_norm_fb_test_signdiff <- subset(rna_norm_fb_test, rna_norm_fb_test$type != 'Not significant')

# Up and down regulated for PANTHER
rna_norm_fb_test_upreg <- merge(rna_norm_fb_test_signdiff[which(rna_norm_fb_test_signdiff$type == 'Upregulated'),], unique(gene_transcript[,c(2,3)]), by.x='row.names', by.y='gene_id')
# write.table(rna_norm_fb_test_upreg$gene_name, quote=F, col.names = F, row.names = F, 
# file='/Volumes/My Passport/hd_in/07.20_hd/3_stdmapping/GO_analysis/upregulated/fb_rna_upreg.tab')

rna_norm_fb_test_dwnreg <- merge(rna_norm_fb_test_signdiff[which(rna_norm_fb_test_signdiff$type == 'Downregulated'),], unique(gene_transcript[,c(2,3)]), by.x='row.names', by.y='gene_id')
# write.table(rna_norm_fb_test_dwnreg$gene_name, quote=F, col.names = F, row.names = F, 
# file='/Volumes/My Passport/hd_in/07.20_hd/3_stdmapping/GO_analysis/downregulated/fb_rna_dwnreg.tab')

rna_norm_fb_test <- merge(rna_norm_fb_test, unique(gene_transcript[,c('gene_id', 'gene_name')]), by.x='Gene id', by.y='gene_id')
# write.table(rna_norm_fb_test$gene_name, quote=F, col.names = F, row.names = F, 
# file='/Volumes/My Passport/hd_in/07.20_hd/3_stdmapping/GO_analysis/fb_rna_expressed.tab')

save.image('fb_RNA_hd_vs_ctrl.RData')

