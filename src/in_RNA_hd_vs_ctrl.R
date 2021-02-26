setwd('/Volumes/My Passport/hd_in/24.02.20/')
load('fb_RNA_hd_vs_ctrl.RData')

coldata_in <- subset(coldata, coldata$CellType == 'IN')

# DESeq DEA and normalizing with median of ratios
dds_in <- DESeqDataSetFromMatrix(rna[,coldata_in$Sample], coldata_in, design= ~Stage)
dds_in$Stage <- relevel(dds_in$Stage, 'CTRL')
dds_in <- DESeq(dds_in)
res_in <- results(dds_in)

rna_norm_in <- as.data.frame(counts(dds_in, normalize=T))
rna_norm_in$gene_id <- rownames(rna_norm_in)

rna_norm_in_test <- apply(rna_norm_in, 1, FUN=function(row){
  row <- as.data.frame(row)
  colnames(row) <- row['gene_id',]
  # print(row)
  row <- row[coldata_in$Sample,,drop=F]
  row <- merge(row, coldata_in[,2,drop=F], by='row.names')
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
    
    hd_expression <- subset(row, row$condition == 'HD')$expression
    ctrl_expression <- subset(row, row$condition == 'CTRL')$expression
    test <- t.test(hd_expression, ctrl_expression, paired=F, exact=F)
    
    results <- c(gene_id, mean_ctrl, mean_hd, log2((mean_hd/mean_ctrl)), normality_ctrl$p.value, normality_hd$p.value, test$statistic, test$p.value, test$conf.int)
    
  }
  else{
    results <- c(gene_id, mean_ctrl, mean_hd, log2((mean_hd/mean_ctrl)), NA, NA, NA, NA, NA, NA)
  }
  
  names(results) <- c("Gene id", "Mean Control", "Mean HD", "Log2FC", 'Shapiro Pvalue - Control', 'Shapiro Pvalue - HD', 'T', 'Pvalue', 'Low conf int', 'High conf int')
  # names(results) <- c("Gene id", "Mean Control", "Mean HD", "Log2FC", 'T', 'Pvalue')
  
  return(results)
})

rna_norm_in_test <- as.data.frame(t(rna_norm_in_test))
rna_norm_in_test$`Gene id` <- as.character(rna_norm_in_test$`Gene id`)
rna_norm_in_test$`Mean Control` <- as.numeric(as.character(rna_norm_in_test$`Mean Control`))
rna_norm_in_test$`Mean HD` <- as.numeric(as.character(rna_norm_in_test$`Mean HD`))
rna_norm_in_test$`Log2FC` <- as.numeric(as.character(rna_norm_in_test$`Log2FC`))
rna_norm_in_test$`Shapiro Pvalue - HD` <- as.numeric(as.character(rna_norm_in_test$`Shapiro Pvalue - HD`))
rna_norm_in_test$`Shapiro Pvalue - Control` <- as.numeric(as.character(rna_norm_in_test$`Shapiro Pvalue - Control`))
rna_norm_in_test$T <- as.numeric(as.character(rna_norm_in_test$T))
rna_norm_in_test$Pvalue <- as.numeric(as.character(rna_norm_in_test$Pvalue))
rna_norm_in_test$`Low conf int` <- as.numeric(as.character(rna_norm_in_test$`Low conf int`))
rna_norm_in_test$`High conf int` <- as.numeric(as.character(rna_norm_in_test$`High conf int`))

# Filter out genes that are not normally distributed
rna_norm_in_test <- rna_norm_in_test[which(rna_norm_in_test$`Shapiro Pvalue - HD` > 0.05 & rna_norm_in_test$`Shapiro Pvalue - Control` > 0.05),]
# Get the counts of the ones that are normally distributed
rna_norm_in <- subset(rna_norm_in, rna_norm_in$gene_id %in% rna_norm_in_test$`Gene id`)
# Calculate the log2(mean per condition + 0.5)

# Which one of those is significantly different?
# Tag them in one
rna_norm_in_test$type <- ifelse(rna_norm_in_test$Pvalue < 0.05 & rna_norm_in_test$Log2FC < 0, 'Downregulated', 
                                ifelse(rna_norm_in_test$Pvalue < 0.05 & rna_norm_in_test$Log2FC > 0, 'Upregulated', 'Not significant'))

rna_norm_in_test$colours <- ifelse(rna_norm_in_test$type == 'Downregulated', 'steelblue', 
                                   ifelse(rna_norm_in_test$type == 'Upregulated', 'tomato3', 'black'))

# Size of the points for the mean plot
rna_norm_in_test$cexs <- ifelse(rna_norm_in_test$Pvalue < 0.05 , 1, 0.5)

# png(paste(getwd(), '/plots/rna_in_meanplot.png', sep=''), width = 15, height = 15, units = "cm", res=200)
pdf(paste(getwd(), '/plots/rna_in_meanplot.pdf', sep=''))
plot(log2(rna_norm_in_test$`Mean Control` + 0.5), 
     log2(rna_norm_in_test$`Mean HD` + 0.5), 
     col=rna_norm_in_test$colours, 
     cex=rna_norm_in_test$cexs, 
     pch=16, 
     xlab='log2(mean Control)', 
     ylab='log2(mean HD)',
     main='RNA iN - HD vs Control (p-value < 0.05; |log2FC| > 0)')

legend("bottomright", legend = c(paste("up (",as.numeric(table(rna_norm_in_test$type)["Upregulated"]),")",sep=""),
                                 paste("down (",as.numeric(table(rna_norm_in_test$type)["Downregulated"]),")",sep = ""),
                                 paste("not significant (",as.numeric(table(rna_norm_in_test$type)["Not significant"]),")",sep = "")),
       pch=16,col=c("firebrick3","steelblue4","black"),cex=1)

dev.off()

svg(paste(getwd(), '/plots/rna_in_meanplot.svg', sep=''))
plot(log2(rna_norm_in_test$`Mean Control` + 0.5), 
     log2(rna_norm_in_test$`Mean HD` + 0.5), 
     col=rna_norm_in_test$colours, 
     cex=rna_norm_in_test$cexs, 
     pch=16, 
     xlab='log2(mean Control)', 
     ylab='log2(mean HD)',
     main='RNA iN - HD vs Control (p-value < 0.05; |log2FC| > 0)')

legend("bottomright", legend = c(paste("up (",as.numeric(table(rna_norm_in_test$type)["Upregulated"]),")",sep=""),
                                 paste("down (",as.numeric(table(rna_norm_in_test$type)["Downregulated"]),")",sep = ""),
                                 paste("not significant (",as.numeric(table(rna_norm_in_test$type)["Not significant"]),")",sep = "")),
       pch=16,col=c("firebrick3","steelblue4","black"),cex=1)

dev.off()

# Significantly different 
rna_norm_in_test_sign <- subset(rna_norm_in_test, rna_norm_in_test$type != 'Not significant')

# Up and down regulated for PANTHER
rna_norm_in_test_upreg <- merge(rna_norm_in_test[which(rna_norm_in_test$type == 'Upregulated'),], unique(gene_transcript[,c(2,3)]), by.x='row.names', by.y='gene_id')
write.table(rna_norm_in_test_upreg$gene_name, quote=F, col.names = F, row.names = F,
            file=paste(getwd(), '/3_stdmapping/GO_analysis/upregulated/in_rna_hd.ctrl_upreg.tab', sep=''))

rna_norm_in_test_dwnreg <- merge(rna_norm_in_test[which(rna_norm_in_test$type == 'Downregulated'),], unique(gene_transcript[,c(2,3)]), by.x='row.names', by.y='gene_id')
write.table(rna_norm_in_test_dwnreg$gene_name, quote=F, col.names = F, row.names = F,
            file=paste(getwd(), '/3_stdmapping/GO_analysis/downregulated/in_rna_hd.ctrl_dwnreg.tab', sep=''))

rna_norm_in_test <- merge(rna_norm_in_test, unique(gene_transcript[,c('gene_id', 'gene_name')]), by.x='Gene id', by.y='gene_id')
write.table(rna_norm_in_test$gene_name, quote=F, col.names = F, row.names = F,
            file=paste(getwd(), '/3_stdmapping/GO_analysis/in_rna_hd.ctrl_expressed.tab', sep=''))

save.image('in_RNA_hd_vs_ctrl.RData')

