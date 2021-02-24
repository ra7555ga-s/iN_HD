setwd('/Volumes/My Passport/hd_in/24.02.20/')
load('RNA_GOanalysis.RData')

protein <- fread('/Volumes/My Passport/hd_in/09.19_hd/5_proteomics/20200309_HDProteomics_Log2_SubtractMedianColumn_AverageSample_+20.txt', data.table=F)
protein <- protein[,c('Gene', colnames(protein)[which(colnames(protein) %in% colnames(rna))])]
protein$gene_unique <- make.unique(protein$Gene)

protein_fb <- protein[, c(colnames(protein)[which(colnames(protein) %in% subset(coldata, coldata$CellType == 'FB')$Sample)], 'Gene', 'gene_unique')]
colnames(protein_fb)[(ncol(protein_fb)-1):ncol(protein_fb)] <- c('gene_name', 'gene_unique')
protein_in <- protein[, c(colnames(protein)[which(colnames(protein) %in% subset(coldata, coldata$CellType == 'IN')$Sample)], 'Gene', 'gene_unique')]
colnames(protein_in)[(ncol(protein_in)-1):ncol(protein_in)] <- c('gene_name', 'gene_unique')

protein_in[is.na(protein_in)] <- 0
protein_fb[is.na(protein_fb)] <- 0

colnames(protein)[1] <- 'gene_name'
protein[is.na(protein)] <- 0

protein_test <- apply(protein, 1, FUN=function(row){
  row <- as.data.frame(row)
  colnames(row) <- row['gene_name',]
  
  row <- row[coldata$Sample,,drop=F]
  row <- merge(row, coldata[,4,drop=F], by='row.names')
  row[,2] <- as.numeric(as.character(row[,2]))
  
  mean_in <- mean(subset(row, row$CellType == 'IN')[,2])
  mean_fb <- mean(subset(row, row$CellType == 'FB')[,2])
  
  gene_name <- colnames(row)[2]
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
    
    results <- c(gene_name, mean_in, mean_fb, log2((mean_fb/mean_in)), normality_in$p.value, normality_fb$p.value, test$statistic, test$p.value, test$conf.int)
    
  }
  else{
    results <- c(gene_name, mean_in, mean_fb, log2((mean_fb/mean_in)), NA, NA, NA, NA, NA, NA)
  }
  
  names(results) <- c("Gene name", "Mean iN", "Mean FB", "Log2FC", 'Shapiro Pvalue - iN', 'Shapiro Pvalue - FB', 'T', 'Pvalue', 'Low conf int', 'High conf int')
  
  return(results)
})
protein_test <- as.data.frame(t(protein_test))
protein_test$`Gene name` <- as.character(protein_test$`Gene name`)
protein_test$`Mean iN` <- as.numeric(as.character(protein_test$`Mean iN`))
protein_test$`Mean FB` <- as.numeric(as.character(protein_test$`Mean FB`))
protein_test$`Log2FC` <- as.numeric(as.character(protein_test$`Log2FC`))
protein_test$`Shapiro Pvalue - iN` <- as.numeric(as.character(protein_test$`Shapiro Pvalue - iN`))
protein_test$`Shapiro Pvalue - FB` <- as.numeric(as.character(protein_test$`Shapiro Pvalue - FB`))
protein_test$T <- as.numeric(as.character(protein_test$T))
protein_test$Pvalue <- as.numeric(as.character(protein_test$Pvalue))
protein_test$`Low conf int` <- as.numeric(as.character(protein_test$`Low conf int`))
protein_test$`High conf int` <- as.numeric(as.character(protein_test$`High conf int`))
# Filter out genes that are not normally distributed
protein_test <- protein_test[which(protein_test$`Shapiro Pvalue - iN` > 0.05 & protein_test$`Shapiro Pvalue - FB` > 0.05),]
# Get the counts of the ones that are normally distributed
protein <- subset(protein, protein$gene_id %in% protein_test$`Gene name`)

# Which one of those is significantly different?
protein_signdiff <- protein_test[which(protein_test$Pvalue < 0.05),]

# Tag them in one
protein_test$type <- ifelse(protein_test$Pvalue < 0.05 & protein_test$Log2FC < 0, 'Downregulated', 
                            ifelse(protein_test$Pvalue < 0.05 & protein_test$Log2FC > 0, 'Upregulated', 'Not significant'))

# Put them colorrsss
protein_test$colours <- ifelse(protein_test$type == 'Downregulated', 'steelblue', 
                               ifelse(protein_test$type == 'Upregulated', 'tomato3', 'black'))

# Size of the points for the mean plot
protein_test$cexs <- ifelse(rownames(protein_test) %in% rownames(protein_signdiff) , 1, 0.5)

png(paste(getwd(), '/plots/protein_fb.in_meanplot.png', sep=''), res = 200, height = 15, width = 15, units = "cm")
plot(log2(protein_test$`Mean iN`+0.5), 
     log2(protein_test$`Mean FB`+0.5), 
     col=protein_test$colours, 
     cex=0.5, 
     pch=16, 
     xlab='log2(mean iN)', 
     ylab='log2(mean FB)', 
     main='Protein FB vs iN (p-value < 0.05; |log2FC| > 0)')

legend("bottomright", legend = c(paste("up (",as.numeric(table(protein_test$type)["Upregulated"]),")",sep=""),
                                 paste("down (",as.numeric(table(protein_test$type)["Downregulated"]),")",sep = ""),
                                 paste("not significant (",as.numeric(table(protein_test$type)["Not significant"]),")",sep = "")),
       pch=16,col=c("firebrick3","steelblue4","black"),cex=1)

dev.off()

save.image('protein_in_vs_fb.RData')
