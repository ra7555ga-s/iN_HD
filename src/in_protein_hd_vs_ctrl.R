setwd('/Volumes/My Passport/hd_in/24.02.20/')
load('fb_protein_hd_vs_ctrl.RData')

# Protein ----
protein <- fread('/Volumes/My Passport/hd_in/09.19_hd/5_proteomics/20200309_HDProteomics_Log2_SubtractMedianColumn_AverageSample_+20.txt', data.table=F)
protein <- protein[,c('Gene', colnames(protein)[which(colnames(protein) %in% colnames(rna))])]
protein$gene_unique <- make.unique(protein$Gene)

protein_fb <- protein[, c(colnames(protein)[which(colnames(protein) %in% subset(coldata, coldata$CellType == 'FB')$Sample)], 'Gene', 'gene_unique')]
colnames(protein_fb)[(ncol(protein_fb)-1):ncol(protein_fb)] <- c('gene_name', 'gene_unique')
protein_in <- protein[, c(colnames(protein)[which(colnames(protein) %in% subset(coldata, coldata$CellType == 'IN')$Sample)], 'Gene', 'gene_unique')]
colnames(protein_in)[(ncol(protein_in)-1):ncol(protein_in)] <- c('gene_name', 'gene_unique')

protein_in[is.na(protein_in)] <- 0
protein_fb[is.na(protein_fb)] <- 0

autophagy_protein_in_control <- protein_in[which(protein_in$gene_name %in% autophagy), ]
rownames(autophagy_protein_in_control) <- autophagy_protein_in_control$gene_unique
autophagy_protein_in_control <- autophagy_protein_in_control[,subset(coldata, coldata$CellType == 'IN' & coldata$Stage == "CTRL")$Sample]
library(pheatmap)
autophagy_protein_in_control_colannot <- coldata[which(coldata$CellType == 'IN' & coldata$Stage == "CTRL"), c("AgeContinuous", "RealName"), drop=F]
rownames(autophagy_protein_in_control_colannot) <- autophagy_protein_in_control_colannot$RealName

colnames(autophagy_protein_in_control) <- coldata[colnames(autophagy_protein_in_control), "RealName"]

pdf('/Volumes/My Passport/hd_in/24.02.20/plots/protein_in_control_autophagy.pdf', height = 50)
pheatmap(autophagy_protein_in_control, 
         scale='row', 
         annotation_col = autophagy_protein_in_control_colannot[,"AgeContinuous", drop=F], cluster_rows = T,
         fontsize = 7)
dev.off()

# Protein iN ----
protein_in_test <- apply(protein_in, 1, FUN=function(row){
  row <- as.data.frame(row)
  colnames(row) <- row['gene_name',]
  
  row <- row[coldata_in$Sample,,drop=F]
  row <- merge(row, coldata_in[,2,drop=F], by='row.names')
  row[,2] <- as.numeric(as.character(row[,2]))
  
  quantile_ctrl <- quantile(subset(row, row$Stage == 'CTRL')[,2])
  quantile_hd <- quantile(subset(row, row$Stage == 'HD')[,2])
  
  gene_name <- colnames(row)[2]
  colnames(row) <- c('sample', 'expression', 'condition')
  row$expression <- as.numeric(as.character(row$expression))
  
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
    # Rest ----
    hd_expression <- subset(row, row$condition == 'HD')$expression
    ctrl_expression <- subset(row, row$condition == 'CTRL')$expression
    test <- t.test(hd_expression, ctrl_expression, paired=F, exact=F)
    
    results <- c(gene_name, quantile_ctrl, quantile_hd, log2((quantile_hd["50%"]/quantile_ctrl["50%"])), normality_ctrl$p.value, normality_hd$p.value, test$statistic, test$p.value, test$conf.int)
    
  }
  else{
    results <- c(gene_name, quantile_ctrl, quantile_hd, log2((quantile_hd["50%"]/quantile_ctrl["50%"])), NA, NA, NA, NA, NA, NA)
  }
  
  names(results) <- c("Gene name", "0% quartile Control", "25% quartile Control", "50% quartile Control", "75% quartile Control", "100% quartile Control",
                      "0% quartile HD", "25% quartile HD", "50% quartile HD", "75% quartile HD", "100% quartile HD",
                      "Log2FC", 'Shapiro Pvalue - Control', 'Shapiro Pvalue - HD', 'T', 'Pvalue', 'Low conf int', 'High conf int')
  
  return(results)
})

protein_in_test <- as.data.frame(t(protein_in_test))
protein_in_test$`Gene name` <- as.character(protein_in_test$`Gene name`)
protein_in_test$`0% quartile Control` <- as.numeric(as.character(protein_in_test$`0% quartile Control`))
protein_in_test$`25% quartile Control` <- as.numeric(as.character(protein_in_test$`25% quartile Control`))
protein_in_test$`50% quartile Control` <- as.numeric(as.character(protein_in_test$`50% quartile Control`))
protein_in_test$`75% quartile Control` <- as.numeric(as.character(protein_in_test$`75% quartile Control`))
protein_in_test$`100% quartile Control` <- as.numeric(as.character(protein_in_test$`100% quartile Control`))
protein_in_test$`0% quartile HD` <- as.numeric(as.character(protein_in_test$`0% quartile HD`))
protein_in_test$`25% quartile HD` <- as.numeric(as.character(protein_in_test$`25% quartile HD`))
protein_in_test$`50% quartile HD` <- as.numeric(as.character(protein_in_test$`50% quartile HD`))
protein_in_test$`75% quartile HD` <- as.numeric(as.character(protein_in_test$`75% quartile HD`))
protein_in_test$`100% quartile HD` <- as.numeric(as.character(protein_in_test$`100% quartile HD`))
protein_in_test$`Log2FC` <- as.numeric(as.character(protein_in_test$`Log2FC`))
protein_in_test$`Shapiro Pvalue - HD` <- as.numeric(as.character(protein_in_test$`Shapiro Pvalue - HD`))
protein_in_test$`Shapiro Pvalue - Control` <- as.numeric(as.character(protein_in_test$`Shapiro Pvalue - Control`))
protein_in_test$T <- as.numeric(as.character(protein_in_test$T)) # T= Mean difference 
protein_in_test$Pvalue <- as.numeric(as.character(protein_in_test$Pvalue))
protein_in_test$`Low conf int` <- as.numeric(as.character(protein_in_test$`Low conf int`))
protein_in_test$`High conf int` <- as.numeric(as.character(protein_in_test$`High conf int`))

# Filter out genes that are not normally distributed
protein_in_test <- protein_in_test[which(protein_in_test$`Shapiro Pvalue - HD` > 0.05 & protein_in_test$`Shapiro Pvalue - Control` > 0.05),]

# Which one of those is significantly different?
protein_signdiff_in <- protein_in_test[which(protein_in_test$Pvalue < 0.05 ),]

# Tag them in one
protein_in_test$type <- ifelse(protein_in_test$Pvalue < 0.05 & protein_in_test$Log2FC < 0, 'Downregulated', 
                               ifelse(protein_in_test$Pvalue < 0.05 & protein_in_test$Log2FC > 0, 'Upregulated', 'Not significant'))

# Put them colorrsss
protein_in_test$colours <- ifelse(protein_in_test$type == 'Downregulated', 'steelblue', 
                                  ifelse(protein_in_test$type == 'Upregulated', 'tomato3', 'black'))

# Size of the points for the mean plot
protein_in_test$cexs <- ifelse(protein_in_test$type != "Not significant" , 1, 0.5)

# png(paste(getwd(), '/plots/in_protein_hd.ctrl_meanplot.png', sep=''), width = 15, height = 15, units = "cm", res=200)
pdf(paste(getwd(), '/plots/in_protein_hd.ctrl_meanplot.pdf', sep=''))
plot(log2(protein_in_test$`50% quartile Control`), 
     log2(protein_in_test$`50% quartile HD`), 
     col=protein_in_test$colours, 
     cex=protein_in_test$cexs, 
     pch=16, 
     xlab='log2(mean Control)', 
     ylab='log2(mean HD)',
     main='Protein iN - HD vs Control (p-value < 0.05; |log2FC| > 0)',
     ylim=c(3.7, 5), xlim=c(3.8, 5))#, ylim=c(14,30), xlim=c(14,30))

legend("bottomright", legend = c(paste("up (",as.numeric(table(protein_in_test$type)["Upregulated"]),")",sep=""),
                                 paste("down (",as.numeric(table(protein_in_test$type)["Downregulated"]),")",sep = ""),
                                 paste("not significant (",as.numeric(table(protein_in_test$type)["Not significant"]),")",sep = "")),
       pch=16,col=c("firebrick3","steelblue4","black"),cex=1)

dev.off()

svg(paste(getwd(), '/plots/in_protein_hd.ctrl_meanplot.svg', sep=''))
plot(log2(protein_in_test$`50% quartile Control`), 
     log2(protein_in_test$`50% quartile HD`), 
     col=protein_in_test$colours, 
     cex=protein_in_test$cexs, 
     pch=16, 
     xlab='log2(mean Control)', 
     ylab='log2(mean HD)',
     main='Protein iN - HD vs Control (p-value < 0.05; |log2FC| > 0)',
     ylim=c(3.7, 5), xlim=c(3.8, 5))#, ylim=c(14,30), xlim=c(14,30))

legend("bottomright", legend = c(paste("up (",as.numeric(table(protein_in_test$type)["Upregulated"]),")",sep=""),
                                 paste("down (",as.numeric(table(protein_in_test$type)["Downregulated"]),")",sep = ""),
                                 paste("not significant (",as.numeric(table(protein_in_test$type)["Not significant"]),")",sep = "")),
       pch=16,col=c("firebrick3","steelblue4","black"),cex=1)

dev.off()

# protein_in_test_upreg <- protein_in_test[which(protein_in_test$type == 'Upregulated'),]
# protein_in_test_upreg <- merge(protein[,c('Gene', 'gene_unique')], protein_in_test_upreg, by.y='Gene name', by.x='gene_unique')
# write.table(protein_in_test_upreg$`Gene`, quote=F, col.names = F, row.names = F,
#             file=paste(getwd(), '/5_proteomics/GO_analysis/upregulated/in_protein_hd.ctrl_upreg.tab', sep=''))
# 
# protein_in_test_dwnreg <- protein_in_test[which(protein_in_test$type == 'Downregulated'),]
# protein_in_test_dwnreg <- merge(protein[,c('Gene', 'gene_unique')], protein_in_test_dwnreg, by.y='Gene name', by.x='gene_unique')
# write.table(protein_in_test_dwnreg$`Gene`, quote=F, col.names = F, row.names = F,
#             file=paste(getwd(), '/5_proteomics/GO_analysis/downregulated/in_protein_hd.ctrl_dwnreg.tab', sep=''))
# 
# protein_in_test <- merge(protein[,c('Gene', 'gene_unique')], protein_in_test, by.y='Gene name', by.x='gene_unique')
# write.table(protein_in_test$`Gene`, quote=F, col.names = F, row.names = F,
#             file=paste(getwd(), '/5_proteomics/GO_analysis/in_protein_hd.ctrl_expressed.tab', sep=''))

save.image('in_protein_hd_vs_ctrl.RData')

