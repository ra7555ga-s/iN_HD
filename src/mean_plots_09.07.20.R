# Settings ----
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(DESeq2)
library(stringr)
library(rjson)
library(xlsx)
library(Hmisc)

id_to_name <- function(input_gene_ids, transcript_gene){
  tmp <- as.data.frame(str_split(input_gene_ids, ', ')[[1]])
  colnames(tmp) <- 'gene_id_panther'
  tmp2 <- merge(tmp, unique(transcript_gene[,c(3,5)]), by='gene_id_panther')
  tmp2 <- tmp2[match(tmp$gene_id_panther, tmp2$gene_id_panther),]
  colnames(tmp2) <- c('gene_id', 'gene_name')
  output_gene_names <- paste(tmp2$gene_name, collapse=', ')
  return(output_gene_names)
}

GO_json <- function(file_name, prefix, botlog2fc=NULL, toplog2fc=NULL, morethan=NULL, ttl='', pdf_width = NULL, pdf_heigth = NULL){
  json <- fromJSON(file=file_name)
  json$overrepresentation$group[[1]]$result[[1]]$input_list$number_in_list
  results <- vector("list",  length(json$overrepresentation$group))
  for(j in 1:length(results)){
    
    for(i in 1:length(json$overrepresentation$group[[j]]$result)){
      
      if(length(json$overrepresentation$group[[j]]$result[[i]]) > 1){
        if(!is.null(json$overrepresentation$group[[j]]$result[[i]]$input_list$number_in_list)){
          if(json$overrepresentation$group[[j]]$result[[i]]$input_list$number_in_list > 0 ){
            tmp <- data.frame(GO_label= json$overrepresentation$group[[j]]$result[[i]]$term$label,
                              GO_id = json$overrepresentation$group[[j]]$result[[i]]$term$id,
                              gene_id=paste(json$overrepresentation$group[[j]]$result[[i]]$input_list$mapped_id_list$mapped_id, collapse=', '),
                              fold_enrichment=json$overrepresentation$group[[j]]$result[[i]]$input_list$fold_enrichment,
                              sign=json$overrepresentation$group[[j]]$result[[i]]$input_list$plus_minus,
                              expected=json$overrepresentation$group[[j]]$result[[i]]$input_list$expected,
                              fdr=json$overrepresentation$group[[j]]$result[[i]]$input_list$fdr,
                              pval = json$overrepresentation$group[[j]]$result[[i]]$input_list$pValue)
            
            
            results[[j]] <- rbind(results[[j]], tmp)
          }
        }
      }
    }
    
  }
  
  results <- results[-which(sapply(results, is.null))]
  results <- results[which(sapply(results, nrow) != 0)]
  results_df <- do.call(rbind, results)
  
  results_df$log2fold_enrichment <- log2(results_df$fold_enrichment+0.5)
  
  results_df$gene_id <- as.character(results_df$gene_id)
  
  if(startsWith(results_df$gene_id[1], 'ENS')){
    transcript_gene <- fread('/Volumes/Seagate Backup /annotation/human/gencode/gencode.v30.annotation.transcr.gene.tab', data.table = FALSE, header=FALSE)
    colnames(transcript_gene) <- c('transcript_id', 'gene_id', 'gene_name', 'gene_type')
    
    transcript_gene$gene_id_panther <- unlist(lapply(str_split(transcript_gene$gene_id, '[.]'), `[[` , 1))
    
    results_df$gene_name <- unlist(lapply(results_df$gene_id, id_to_name, transcript_gene=transcript_gene))
  }
  
  write.xlsx(x=results_df[order(results_df$fdr), ], file=paste(prefix, '.xlsx', sep = ''), col.names = T, row.names = F)
  
  results_df_filtered <- results_df
  
  if(!is.null(toplog2fc)){
    if(!is.null(botlog2fc)){
      results_df_filtered <- subset(results_df, results_df$log2fold_enrichment > toplog2fc | results_df$log2fold_enrichment < botlog2fc)
    }
    else{
      results_df_filtered <- subset(results_df, results_df$log2fold_enrichment > toplog2fc)
    }
  }else{
    if(!is.null(botlog2fc)){
      results_df_filtered <- subset(results_df, results_df$log2fold_enrichment < botlog2fc)
    }
  }
  if(!is.null(morethan)){
    results_df_filtered <- results_df_filtered[which(unlist(lapply(str_split(results_df_filtered$gene_id, ', '), length)) > morethan),]
  }
  # results_df_filtered <- subset(results_df_filtered, results_df_filtered$pval < pvalue)
  
  results_df_filtered$GO_label <- as.character(results_df_filtered$GO_label)
  for(i in 1:nrow(results_df_filtered)){
    if(length(unlist(str_split(as.character(results_df_filtered[i,'GO_label']), " "))) > 10){
      newline <- ceiling(length(unlist(str_split(as.character(results_df_filtered[i,'GO_label']), " ")))/2)
      no_words <- floor(length(unlist(str_split(as.character(results_df_filtered[i,'GO_label']), " ")))) 
      results_df_filtered[i,'GO_label'] <- paste(c(unlist(str_split(as.character(results_df_filtered[i,'GO_label']), " "))[1:newline], "\n", unlist(str_split(as.character(results_df_filtered[i,'GO_label']), " "))[(newline+2):no_words]), collapse = " ")
    }  
  }
  results_df_filtered$fdr_log10 <- -log10(results_df_filtered$fdr)
  results_df_filtered$num_of_genes <- unlist(lapply(str_split(results_df_filtered$gene_id, ', '), FUN=length))
  
  results_for_plot <- results_df_filtered 
  if("TRUE" %in% names(table(results_df_filtered$fdr < 0.05)) & nrow(results_df_filtered) > 100){
    results_for_plot <- results_df_filtered[which(results_df_filtered$fdr < 0.05),]
    if(nrow(results_for_plot) > 30)
    {
      plot_title <- "Top 30 significant GO terms"
      results_for_plot <- results_for_plot[1:30,]
    }
    else
    {
      plot_title <- "All significant GO terms"
      if(nrow(results_for_plot) < 3)
      {
        results_for_plot <- results_df_filtered[which(results_df_filtered$fdr < 0.05 | results_df_filtered$log2fold_enrichment > 1),]
        results_for_plot <- results_for_plot[order(results_for_plot$log2fold_enrichment),]
        results_for_plot <- results_for_plot[1:30,]
        plot_title <- "Top 30 GO terms"
      }
    }
  }else{
    if(nrow(results_for_plot) > 30)
    {
      results_for_plot <- results_for_plot[order(results_for_plot$fdr),]
      results_for_plot <- results_for_plot[1:30,]
      plot_title <- "Top 30 GO terms"
    }else{
      plot_title <- "All GO terms"
    }
  }
  
  steps <- ifelse(max(results_for_plot$num_of_genes) > 1, ceiling((max(results_for_plot$num_of_genes) - min(results_for_plot$num_of_genes))/5), 1)
  if(steps > 1){
    results_for_plot$num_of_genes <- cut2(results_for_plot$num_of_genes, seq(min(results_for_plot$num_of_genes), max(results_for_plot$num_of_genes), steps), digits=0, oneval=FALSE)
    results_for_plot$num_of_genes <- sub(',', '-', gsub('^.|.$', '', as.character(results_for_plot$num_of_genes)))
  }
  
  library(RColorBrewer)
  getPalette = brewer.pal(length(unique(results_for_plot$num_of_genes))+1, "Blues")[-1]
  results_for_plot$GO_label <- factor(results_for_plot$GO_label, levels=results_for_plot[order(results_for_plot$fdr, decreasing = T),'GO_label'])
  results_for_plot$num_of_genes <- as.factor(results_for_plot$num_of_genes)
  
  plot <- ggplot(data=results_for_plot, aes(x=GO_label, y=log2fold_enrichment)) +
    geom_bar(stat="identity", aes(fill=num_of_genes))+
    expand_limits(y=c(min(results_for_plot$log2fold_enrichment), (max(results_for_plot$log2fold_enrichment)+0.5)))+
    xlab("GO term") + ylab("\nlog2(Fold Change)")+ theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip() +
    geom_hline(yintercept=0, color='black')+
    ggtitle(plot_title) +
    geom_point(aes(y=fdr_log10), size=3) + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "-log10(padj)"))+
    geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red') + labs(fill="Num. of genes") +
    theme(text = element_text(size=24), axis.text.x = element_text(angle=90, hjust=1)) 
  
  if(length(unique(results_for_plot$num_of_genes)) > 1){
    plot <- plot + scale_fill_manual(values=getPalette)
  }
  
  ggsave(plot, file=paste(prefix, 'fdr.png', sep = ''), width=32, height=length(results_for_plot$GO_label)/1.2, units="cm")
  
  if(!is.null(pdf_width)){
    if(!is.null(pdf_heigth)){
      pdf(paste(prefix, 'fdr.pdf', sep = ''), height=pdf_heigth, width = pdf_width)
      print(plot)
      dev.off()  
    }else{
      pdf(paste(prefix, 'fdr.pdf', sep = ''), height=length(results_for_plot$GO_label)/1.7, width = pdf_width)
      print(plot)
      dev.off()
    }
  }else{
    if(!is.null(pdf_heigth)){
      pdf(paste(prefix, 'fdr.pdf', sep = ''), height=pdf_heigth, width = 25)
      print(plot)
      dev.off()
    }else{
      pdf(paste(prefix, 'fdr.pdf', sep = ''), height=length(results_for_plot$GO_label)/1.7, width = 25)
      print(plot)
      dev.off()
    }
  }
  
  
  return(plot)
}
more_10 <- function(row){
  return(length(which(row > 10)))
}

# Coldata ----
coldata <- fread('/Volumes/Seagate Backup /hd_in/09.19_hd/CategoricalAnnotation.txt', data.table = F)
rownames(coldata) <- coldata$Sample
coldata$Sample <- as.character(coldata$Sample)
# RNAseq ----
rna <- fread('/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/1_readcounts/gene_count_matrix_2.csv', data.table=F)
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

dds_ct_rna <- DESeqDataSetFromMatrix(rna_expressed[,coldata$Sample], coldata, design= ~CellType)
dds_ct_rna$CellType <- relevel(dds_ct_rna$CellType, 'FB')
dds_ct_rna <- DESeq(dds_ct_rna)
vst_ct_rna <- varianceStabilizingTransformation(dds_ct_rna)
plotPCA(vst_ct_rna, intgroup='CellType')
dds_ct_rna$sizeFactor

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

png('/Volumes/Seagate Backup /hd_in/07.20_hd/plots/rna_fb_meanplot.png', width = 15, height = 15, units = "cm", res=200)
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

dev.off()

# Significantly different 
rna_norm_fb_test_signdiff <- subset(rna_norm_fb_test, rna_norm_fb_test$type != 'Not significant')

gene_transcript <- fread('/Volumes/Seagate Backup /annotation/human/gencode/gencode.v30.annotation.transcr.gene.tab', data.table=F, header=F)
colnames(gene_transcript) <- c('transcript_id', 'gene_id', 'gene_name', 'gene_type') 

# Up and down regulated for PANTHER
rna_norm_fb_test_upreg <- merge(rna_norm_fb_test_signdiff[which(rna_norm_fb_test_signdiff$type == 'Upregulated'),], unique(gene_transcript[,c(2,3)]), by.x='row.names', by.y='gene_id')
write.table(rna_norm_fb_test_upreg$gene_name, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/upregulated/fb_rna_upreg.tab')

rna_norm_fb_test_dwnreg <- merge(rna_norm_fb_test_signdiff[which(rna_norm_fb_test_signdiff$type == 'Downregulated'),], unique(gene_transcript[,c(2,3)]), by.x='row.names', by.y='gene_id')
write.table(rna_norm_fb_test_dwnreg$gene_name, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/downregulated/fb_rna_dwnreg.tab')

rna_norm_fb_test <- merge(rna_norm_fb_test, unique(gene_transcript[,c('gene_id', 'gene_name')]), by.x='Gene id', by.y='gene_id')
write.table(rna_norm_fb_test$gene_name, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/fb_rna_expressed.tab')

# RNAseq iNs ----
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

# Put them colorrsss
rna_norm_in_test$colours <- ifelse(rna_norm_in_test$type == 'Downregulated', 'steelblue', 
                                   ifelse(rna_norm_in_test$type == 'Upregulated', 'tomato3', 'black'))

# Size of the points for the mean plot
rna_norm_in_test$cexs <- ifelse(rna_norm_in_test$Pvalue < 0.05 , 1, 0.5)

png('/Volumes/Seagate Backup /hd_in/07.20_hd/plots/rna_in_meanplot.png', width = 15, height = 15, units = "cm", res=200)
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
            file='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/upregulated/in_rna_upreg.tab')

rna_norm_in_test_dwnreg <- merge(rna_norm_in_test[which(rna_norm_in_test$type == 'Downregulated'),], unique(gene_transcript[,c(2,3)]), by.x='row.names', by.y='gene_id')
write.table(rna_norm_in_test_dwnreg$gene_name, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/downregulated/in_rna_dwnreg.tab')

rna_norm_in_test <- merge(rna_norm_in_test, unique(gene_transcript[,c('gene_id', 'gene_name')]), by.x='Gene id', by.y='gene_id')
write.table(rna_norm_in_test$gene_name, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/in_rna_expressed.tab')


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
# rna_norm_test <- rna_norm_test[which(rna_norm_test$`Shapiro Pvalue - HD` > 0.05 & rna_norm_test$`Shapiro Pvalue - Control` > 0.05),]
rna_norm_test <- rna_norm_test[which(rna_norm_test$`Shapiro Pvalue iN` > 0.05 & rna_norm_test$`Shapiro Pvalue FB` > 0.05),]
# Get the counts of the ones that are normally distributed
rna_expressed_norm <- subset(rna_expressed_norm, rna_expressed_norm$gene_id %in% rna_norm_test$`Gene id`)

# Which one of those is significantly different?
rna_signdiff <- rna_norm_test[which(rna_norm_test$Pvalue < 0.05),]
# Tag them in one
rna_norm_test$type <- ifelse(rna_norm_test$Pvalue < 0.05 & rna_norm_test$Log2FC < 0, 'Downregulated', 
                                ifelse(rna_norm_test$Pvalue < 0.05 & rna_norm_test$Log2FC > 0, 'Upregulated', 'Not significant'))

# Put them colorrsss
rna_norm_test$colours <- ifelse(rna_norm_test$type == 'Downregulated', 'steelblue', 
                                   ifelse(rna_norm_test$type == 'Upregulated', 'tomato3', 'black'))

png('/Volumes/Seagate Backup /hd_in/07.20_hd/plots/rna_fb.in_meanplot.png',res=200, width=15, height=15,units = "cm")
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

# GO analysis ----
# Biological process
GO_json(file_name='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/upregulated/in_rna_upreg_biological_process.json', prefix='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/upregulated/in_rna_upreg_biological_process_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 10)
GO_json(file_name='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/downregulated/in_rna_downreg_biological_process.json', prefix='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/downregulated/in_rna_downreg_biological_process_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 15, pdf_heigth= 15)

# Molecular function
GO_json(file_name='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/upregulated/in_rna_upreg_molecular_function.json', prefix='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/upregulated/in_rna_upreg_molecular_function_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 15)
GO_json(file_name='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/downregulated/in_rna_downreg_molecular_function.json', prefix='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/downregulated/in_rna_downreg_molecular_function_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth = 15)

# Biological process
GO_json(file_name='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/upregulated/fb_rna_upreg_biological_process.json', prefix='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/upregulated/fb_rna_upreg_biological_process_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 10)
GO_json(file_name='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/downregulated/fb_rna_downreg_biological_process.json', prefix='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/downregulated/fb_rna_downreg_biological_process_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 15, pdf_heigth= 15)

# Molecular function
GO_json(file_name='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/upregulated/fb_rna_upreg_molecular_function.json', prefix='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/upregulated/fb_rna_upreg_molecular_function_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 15)
GO_json(file_name='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/downregulated/fb_rna_downreg_molecular_function.json', prefix='/Volumes/Seagate Backup /hd_in/07.20_hd/3_stdmapping/GO_analysis/downregulated/fb_rna_downreg_molecular_function_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth = 15)

# Protein ----
protein <- fread('/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/20200309_HDProteomics_Log2_SubtractMedianColumn_AverageSample_+20.txt', data.table=F)
protein <- protein[,c('Gene', colnames(protein)[which(colnames(protein) %in% colnames(rna))])]
protein$gene_unique <- make.unique(protein$Gene)

protein_fb <- protein[, c(colnames(protein)[which(colnames(protein) %in% subset(coldata, coldata$CellType == 'FB')$Sample)], 'Gene', 'gene_unique')]
colnames(protein_fb)[(ncol(protein_fb)-1):ncol(protein_fb)] <- c('gene_name', 'gene_unique')
protein_in <- protein[, c(colnames(protein)[which(colnames(protein) %in% subset(coldata, coldata$CellType == 'IN')$Sample)], 'Gene', 'gene_unique')]
colnames(protein_in)[(ncol(protein_in)-1):ncol(protein_in)] <- c('gene_name', 'gene_unique')

protein_in[is.na(protein_in)] <- 0
protein_fb[is.na(protein_fb)] <- 0

# Protein FB ----
# For each gene (normalized and filtered by expression) 
# 1) Shapiro test
# 2) T test
# rownames(protein_fb) <- make.unique(protein_fb$gene_name)
# row <- t(protein_fb[1,])
# protein_fb$hd_mean <- rowMeans(protein_fb[,startsWith(colnames(protein_fb), 'HD')])
# protein_fb$ctrl_mean <- rowMeans(protein_fb[,!startsWith(colnames(protein_fb), 'HD')])

protein_fb_test <- apply(protein_fb, 1, FUN=function(row){
  row <- as.data.frame(row)
  colnames(row) <- row['gene_unique',]
  
  row <- row[coldata_fb$Sample,,drop=F]
  row <- merge(row, coldata_fb[,2,drop=F], by='row.names')
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
    
    hd_expression <- subset(row, row$condition == 'HD')$expression
    ctrl_expression <- subset(row, row$condition == 'CTRL')$expression
    test <- t.test(hd_expression, ctrl_expression, paired=F, exact=F)
    
    results <- c(gene_name, quantile_ctrl, quantile_hd, log2((quantile_hd["50%"]/quantile_ctrl["50%"])), 
                 normality_ctrl$p.value, 
                 normality_hd$p.value, test$statistic, test$p.value, test$conf.int)
    
  }
  else{
    results <- c(gene_name, quantile_ctrl, quantile_hd, log2((quantile_hd["50%"]/quantile_ctrl["50%"])), NA, NA, NA, NA, NA, NA)
  }
  
  
  names(results) <- c("Gene name", "0% quartile Control", "25% quartile Control", "50% quartile Control", "75% quartile Control", "100% quartile Control",
                      "0% quartile HD", "25% quartile HD", "50% quartile HD", "75% quartile HD", "100% quartile HD",
                      "Log2FC", 'Shapiro Pvalue - Control', 'Shapiro Pvalue - HD', 'T', 'Pvalue', 'Low conf int', 'High conf int')
  
  return(results)
})

protein_fb_test <- as.data.frame(t(protein_fb_test))
protein_fb_test$`Gene name` <- as.character(protein_fb_test$`Gene name`)
protein_fb_test$`0% quartile Control` <- as.numeric(as.character(protein_fb_test$`0% quartile Control`))
protein_fb_test$`25% quartile Control` <- as.numeric(as.character(protein_fb_test$`25% quartile Control`))
protein_fb_test$`50% quartile Control` <- as.numeric(as.character(protein_fb_test$`50% quartile Control`))
protein_fb_test$`75% quartile Control` <- as.numeric(as.character(protein_fb_test$`75% quartile Control`))
protein_fb_test$`100% quartile Control` <- as.numeric(as.character(protein_fb_test$`100% quartile Control`))
protein_fb_test$`0% quartile HD` <- as.numeric(as.character(protein_fb_test$`0% quartile HD`))
protein_fb_test$`25% quartile HD` <- as.numeric(as.character(protein_fb_test$`25% quartile HD`))
protein_fb_test$`50% quartile HD` <- as.numeric(as.character(protein_fb_test$`50% quartile HD`))
protein_fb_test$`75% quartile HD` <- as.numeric(as.character(protein_fb_test$`75% quartile HD`))
protein_fb_test$`100% quartile HD` <- as.numeric(as.character(protein_fb_test$`100% quartile HD`))
protein_fb_test$`Log2FC` <- as.numeric(as.character(protein_fb_test$`Log2FC`))
protein_fb_test$`Shapiro Pvalue - HD` <- as.numeric(as.character(protein_fb_test$`Shapiro Pvalue - HD`))
protein_fb_test$`Shapiro Pvalue - Control` <- as.numeric(as.character(protein_fb_test$`Shapiro Pvalue - Control`))
protein_fb_test$T <- as.numeric(as.character(protein_fb_test$T)) # T= Mean difference 
protein_fb_test$Pvalue <- as.numeric(as.character(protein_fb_test$Pvalue))
protein_fb_test$`Low conf int` <- as.numeric(as.character(protein_fb_test$`Low conf int`))
protein_fb_test$`High conf int` <- as.numeric(as.character(protein_fb_test$`High conf int`))

# Filter out genes that are not normally distributed
protein_fb_test <- protein_fb_test[which(protein_fb_test$`Shapiro Pvalue - HD` > 0.05 & protein_fb_test$`Shapiro Pvalue - Control` > 0.05),]

# Which one of those is significantly different?
# Tag them in one
protein_fb_test$type <- ifelse(protein_fb_test$Pvalue < 0.05 & protein_fb_test$Log2FC < 0, 'Downregulated', 
                               ifelse(protein_fb_test$Pvalue < 0.05 & protein_fb_test$Log2FC > 0, 'Upregulated', 'Not significant'))

# Put them colorrsss
protein_fb_test$colours <- ifelse(protein_fb_test$type == 'Downregulated', 'steelblue', 
                                  ifelse(protein_fb_test$type == 'Upregulated', 'tomato3', 'black'))

# Size of the points for the mean plot
protein_fb_test$cexs <- ifelse(protein_fb_test$Pvalue < 0.05 , 1, 0.5)

protein_fb_test$`log2 Mean Control` <- log2(protein_fb_test$`50% quartile Control` + 0.5) 
protein_fb_test$`log2 Mean HD` <- log2(protein_fb_test$`50% quartile HD` + 0.5)

png('/Volumes/Seagate Backup /hd_in/07.20_hd/plots/protein_fb_meanplot.png', width = 15, height = 15, units = "cm", res=200)
plot(protein_fb_test$`log2 Mean Control`, 
     protein_fb_test$`log2 Mean HD`, 
     col=protein_fb_test$colours, 
     cex=protein_fb_test$cexs, 
     pch=16, 
     xlab='log2(mean Control)', 
     ylab='log2(mean HD)',
     main='Protein FB - HD vs Control (p-value < 0.05; |log2FC| > 0)')

legend("bottomright", legend = c(paste("up (",as.numeric(table(protein_fb_test$type)["Upregulated"]),")",sep=""),
                                 paste("down (",as.numeric(table(protein_fb_test$type)["Downregulated"]),")",sep = ""),
                                 paste("not significant (",as.numeric(table(protein_fb_test$type)["Not significant"]),")",sep = "")),
       pch=16,col=c("firebrick3","steelblue4","black"),cex=1)

dev.off()

# Significantly different 
# Means
protein_fb_test_signdiff <- subset(protein_fb_test, protein_fb_test$type != 'Not significant')

# Up and down regulated for PANTHER
protein_fb_test_upreg <- protein_fb_test_signdiff[which(protein_fb_test_signdiff$type == 'Upregulated'),]
protein_fb_test_upreg <- merge(protein[,c('Gene', 'gene_unique')], protein_fb_test_upreg, by.y='Gene name', by.x='gene_unique')
write.table(protein_fb_test_upreg$Gene, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/upregulated/fb_protein_upreg.tab')

protein_fb_test_dwnreg <- protein_fb_test_signdiff[which(protein_fb_test_signdiff$type == 'Downregulated'),]
protein_fb_test_dwnreg <- merge(protein[,c('Gene', 'gene_unique')], protein_fb_test_dwnreg, by.y='Gene name', by.x='gene_unique')
write.table(protein_fb_test_dwnreg$Gene, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/downregulated/fb_protein_dwnreg.tab')

write.table(protein_fb_test$Gene, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/fb_protein_expressed.tab')

write.xlsx(protein_fb_test_signdiff[,-c((ncol(protein_fb_test_signdiff)-3):ncol(protein_fb_test_signdiff))], 
           file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/fb_protein_signdiff.xlsx', row.names = F)

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

png('/Volumes/Seagate Backup /hd_in/07.20_hd/plots/protein_in_meanplot.png', width = 15, height = 15, units = "cm", res=200)
plot(log2(protein_in_test$`50% quartile Control`), 
     log2(protein_in_test$`50% quartile HD`), 
     col=protein_in_test$colours, 
     cex=protein_in_test$cexs, 
     pch=16, 
     xlab='log2(mean Control)', 
     ylab='log2(mean HD)',
     main='Protein iN - HD vs Control (p-value < 0.05; |log2FC| > 0)')#, ylim=c(14,30), xlim=c(14,30))

legend("bottomright", legend = c(paste("up (",as.numeric(table(protein_in_test$type)["Upregulated"]),")",sep=""),
                                 paste("down (",as.numeric(table(protein_in_test$type)["Downregulated"]),")",sep = ""),
                                 paste("not significant (",as.numeric(table(protein_in_test$type)["Not significant"]),")",sep = "")),
       pch=16,col=c("firebrick3","steelblue4","black"),cex=1)

dev.off()

# Up and down regulated for PANTHER ----
protein_in_test_upreg <- protein_in_test[which(protein_in_test$type == 'Upregulated'),]
protein_in_test_upreg <- merge(protein[,c('Gene', 'gene_unique')], protein_in_test_upreg, by.y='Gene name', by.x='gene_unique')
write.table(protein_in_test_upreg$`Gene`, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/upregulated/in_protein_upreg.tab')

protein_in_test_dwnreg <- protein_in_test[which(protein_in_test$type == 'Downregulated'),]
protein_in_test_dwnreg <- merge(protein[,c('Gene', 'gene_unique')], protein_in_test_dwnreg, by.y='Gene name', by.x='gene_unique')
write.table(protein_in_test_dwnreg$`Gene`, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/downregulated/in_protein_dwnreg.tab')

protein_in_test <- merge(protein[,c('Gene', 'gene_unique')], protein_in_test, by.y='Gene name', by.x='gene_unique')
write.table(protein_in_test$`Gene`, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/in_protein_expressed.tab')


# Significantly diff in INs, expression in all datasets ----
protein_in_test_signdiff <- subset(protein_in_test, protein_in_test$type != 'Not significant')
write.xlsx(protein_in_test_signdiff[,-c((ncol(protein_in_test_signdiff)-3):ncol(protein_in_test_signdiff))], 
           file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/in_protein_signdiff.xlsx', row.names = F)

# protein_in_test_signdiff <- protein_in_test_signdiff[,c(1:19)]
colnames(protein_in_test_signdiff)[3:12] <- paste("Protein iN -", colnames(protein_in_test_signdiff)[3:12])

rownames(protein_fb) <- protein_fb$gene_unique
# Expression of significantly different genes in FB
protein_in_signdiff_fb <- protein_fb[protein_in_test_signdiff$gene_unique,]

# Signficantly different in iNs, quartiles from FB CTRL Proteomics
protein_in_test_signdiff <- merge(protein_in_test_signdiff, 
                                  t(apply(protein_fb[protein_in_signdiff_fb$gene_unique, 
                                                     subset(coldata_fb, coldata_fb$Stage == 'CTRL')$Sample], 1, FUN=quantile)), 
                                  by.x='gene_unique', by.y='row.names', all=T)
colnames(protein_in_test_signdiff)[c((ncol(protein_in_test_signdiff)-4):ncol(protein_in_test_signdiff))] <- c("Protein FB - 0% quartile Control", "Protein FB - 25% quartile Control", "Protein FB - 50% quartile Control", "Protein FB - 75% quartile Control", "Protein FB - 100% quartile Control")

# Signficantly different in iNs, quartiles from FB HD Proteomics
protein_in_test_signdiff <- merge(protein_in_test_signdiff, 
                                  t(apply(protein_fb[protein_in_signdiff_fb$gene_unique, 
                                                     subset(coldata_fb, coldata_fb$Stage == 'HD')$Sample], 1, FUN=quantile)), 
                                  by.x='gene_unique', by.y='row.names', all=T)
colnames(protein_in_test_signdiff)[c((ncol(protein_in_test_signdiff)-4):ncol(protein_in_test_signdiff))] <- c("Protein FB - 0% quartile HD", "Protein FB - 25% quartile HD", "Protein FB - 50% quartile HD", "Protein FB - 75% quartile HD", "Protein FB - 100% quartile HD")

rna_norm_fb <- merge(rna_norm_fb, unique(gene_transcript[,c('gene_id', 'gene_name')]), by='gene_id')
rna_norm_in <- merge(rna_norm_in, unique(gene_transcript[,c('gene_id', 'gene_name')]), by='gene_id')

protein_in_signdiff_rna.fb <- subset(rna_norm_fb, rna_norm_fb$gene_name %in% protein_in_test_signdiff$Gene)
rownames(protein_in_signdiff_rna.fb) <- make.unique(protein_in_signdiff_rna.fb$gene_name)

# Significantly different in iNs, quartiles from FB CTRL RNAseq
protein_in_signdiff_rna.fb <- merge(protein_in_signdiff_rna.fb, 
                                    t(apply(protein_in_signdiff_rna.fb[, subset(coldata_fb, coldata_fb$Stage == 'CTRL')$Sample], 1, FUN=quantile)), 
                                    by='row.names')
rownames(protein_in_signdiff_rna.fb) <- protein_in_signdiff_rna.fb$Row.names
protein_in_signdiff_rna.fb <- protein_in_signdiff_rna.fb[,-c(1,2)]
colnames(protein_in_signdiff_rna.fb)[c((ncol(protein_in_signdiff_rna.fb)-4):ncol(protein_in_signdiff_rna.fb))] <- c("RNA FB - 0% quartile Control", "RNA FB - 25% quartile Control", "RNA FB - 50% quartile Control", "RNA FB - 75% quartile Control", "RNA FB - 100% quartile Control")

# Significantly different in iNs, quartiles from FB HD RNAseq
protein_in_signdiff_rna.fb <- merge(protein_in_signdiff_rna.fb, 
                                    t(apply(protein_in_signdiff_rna.fb[, subset(coldata_fb, coldata_fb$Stage == 'HD')$Sample], 1, FUN=quantile)), 
                                    by='row.names')
rownames(protein_in_signdiff_rna.fb) <- protein_in_signdiff_rna.fb$Row.names
protein_in_signdiff_rna.fb <- protein_in_signdiff_rna.fb[,-c(1)]
colnames(protein_in_signdiff_rna.fb)[c((ncol(protein_in_signdiff_rna.fb)-4):ncol(protein_in_signdiff_rna.fb))] <- c("RNA FB - 0% quartile HD", "RNA FB - 25% quartile HD", "RNA FB - 50% quartile HD", "RNA FB - 75% quartile HD", "RNA FB - 100% quartile HD")

protein_in_test_signdiff <- merge(protein_in_test_signdiff, protein_in_signdiff_rna.fb[,-c(1:14)], by.x='Gene', by.y='gene_name', all=T)
colnames(protein_in_test_signdiff)

protein_in_signdiff_rna.in <- subset(rna_norm_in, rna_norm_in$gene_name %in% protein_in_test_signdiff$Gene)
rownames(protein_in_signdiff_rna.in) <- make.unique(protein_in_signdiff_rna.in$gene_name)

# Significantly different in iNs, quartiles from iN CTRL RNAseq
protein_in_signdiff_rna.in <- merge(protein_in_signdiff_rna.in, 
                                    t(apply(protein_in_signdiff_rna.in[, subset(coldata_in, coldata_in$Stage == 'CTRL')$Sample], 1, FUN=quantile)), 
                                    by='row.names')
rownames(protein_in_signdiff_rna.in) <- protein_in_signdiff_rna.in$Row.names
protein_in_signdiff_rna.in <- protein_in_signdiff_rna.in[,-c(1,2)]
colnames(protein_in_signdiff_rna.in)[c((ncol(protein_in_signdiff_rna.in)-4):ncol(protein_in_signdiff_rna.in))] <- c("RNA iN - 0% quartile Control", "RNA iN - 25% quartile Control", "RNA iN - 50% quartile Control", "RNA iN - 75% quartile Control", "RNA iN - 100% quartile Control")

# Significantly different in iNs, quartiles from iN HD RNAseq
protein_in_signdiff_rna.in <- merge(protein_in_signdiff_rna.in, 
                                    t(apply(protein_in_signdiff_rna.in[, subset(coldata_in, coldata_in$Stage == 'HD')$Sample], 1, FUN=quantile)), 
                                    by='row.names')
rownames(protein_in_signdiff_rna.in) <- protein_in_signdiff_rna.in$Row.names
protein_in_signdiff_rna.in <- protein_in_signdiff_rna.in[,-c(1)]
colnames(protein_in_signdiff_rna.in)[c((ncol(protein_in_signdiff_rna.in)-4):ncol(protein_in_signdiff_rna.in))] <- c("RNA iN - 0% quartile HD", "RNA iN - 25% quartile HD", "RNA iN - 50% quartile HD", "RNA iN - 75% quartile HD", "RNA iN - 100% quartile HD")

protein_in_test_signdiff <- merge(protein_in_test_signdiff, protein_in_signdiff_rna.in[,-c(1:14)], by.x='Gene', by.y='gene_name', all=T)
colnames(protein_in_test_signdiff)

write.xlsx(protein_in_test_signdiff, file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/in_protein_signdiff_alldatasets.xlsx', row.names = F)
write.xlsx(subset(protein_in_test_signdiff, protein_in_test_signdiff$Log2FC < 0), file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/in_protein_signdiff_alldatasets_downreg.xlsx', row.names = F)
write.xlsx(subset(protein_in_test_signdiff, protein_in_test_signdiff$Log2FC > 0), file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/in_protein_signdiff_alldatasets_upreg.xlsx', row.names = F)

length(unique(subset(protein_in_test_signdiff, protein_in_test_signdiff$Log2FC < 0)$Gene))
length(unique(subset(protein_in_test_signdiff, protein_in_test_signdiff$Log2FC > 0)$Gene))

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
IDs <- protein_in_test_dwnreg$Gene
genedesc <- getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = IDs, mart =ensembl)
length(unique(genedesc$external_gene_name))
length(unique(protein_in_test_dwnreg$Gene))
protein_in_test_dwnreg <- merge(genedesc, protein_in_test_dwnreg, by.x='external_gene_name', by.y='Gene name')
colnames(protein_in_test_dwnreg[,-c(1,20:22)])
write.xlsx(protein_in_test_dwnreg[,-c(1,20:22)], col.names = T, row.names = F, 
           file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/downregulated/in_protein_dwnreg.xlsx')

# All proteins, all datasets ----
protein <- fread('/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/20200309_HDProteomics_Log2_SubtractMedianColumn_AverageSample_+20.txt', data.table=F)
protein <- protein[,c('Gene', colnames(protein)[which(colnames(protein) %in% colnames(rna))])]
protein$gene_unique <- make.unique(protein$Gene)
rownames(protein) <- protein$gene_unique
protein[is.na(protein)] <- 0

all_proteins <- t(apply(protein[protein$gene_unique, subset(coldata_in, coldata_in$Stage == 'CTRL')$Sample], 1, FUN=quantile))
colnames(all_proteins) <- paste("Protein iN -", colnames(all_proteins), 'quartile Control')
all_proteins <- merge(all_proteins, protein[,c('Gene', 'gene_unique')], by.x='row.names', by.y='gene_unique')
rownames(all_proteins) <- all_proteins$Row.names
all_proteins <- all_proteins[,-1]

all_proteins <- merge(all_proteins, t(apply(protein[protein$gene_unique, subset(coldata_in, coldata_in$Stage == 'HD')$Sample], 1, FUN=quantile)), by='row.names')
colnames(all_proteins)[c((ncol(all_proteins)-4):ncol(all_proteins))] <- c("Protein iN - 0% quartile HD", "Protein iN - 25% quartile HD", "Protein iN - 50% quartile HD", "Protein iN - 75% quartile HD", "Protein iN - 100% quartile HD")
rownames(all_proteins) <- all_proteins$Row.names
all_proteins <- all_proteins[,-1]

all_proteins <- merge(all_proteins, t(apply(protein[protein$gene_unique, subset(coldata_fb, coldata_fb$Stage == 'CTRL')$Sample], 1, FUN=quantile)), by='row.names')
colnames(all_proteins)[c((ncol(all_proteins)-4):ncol(all_proteins))] <- c("Protein FB - 0% quartile Control", "Protein FB - 25% quartile Control", "Protein FB - 50% quartile Control", "Protein FB - 75% quartile Control", "Protein FB - 100% quartile Control")
rownames(all_proteins) <- all_proteins$Row.names
all_proteins <- all_proteins[,-1]

all_proteins <- merge(all_proteins, t(apply(protein[protein$gene_unique, subset(coldata_fb, coldata_fb$Stage == 'HD')$Sample], 1, FUN=quantile)), by='row.names')
colnames(all_proteins)[c((ncol(all_proteins)-4):ncol(all_proteins))] <- c("Protein FB - 0% quartile HD", "Protein FB - 25% quartile HD", "Protein FB - 50% quartile HD", "Protein FB - 75% quartile HD", "Protein FB - 100% quartile HD")
rownames(all_proteins) <- all_proteins$Row.names
all_proteins <- all_proteins[,-1]

all_proteins_rna_norm_fb <- subset(rna_norm_fb, rna_norm_fb$gene_name %in% all_proteins$Gene)
rownames(all_proteins_rna_norm_fb) <- make.unique(all_proteins_rna_norm_fb$gene_name)

all_proteins_rna_norm_fb <- merge(all_proteins_rna_norm_fb, 
                                    t(apply(all_proteins_rna_norm_fb[, subset(coldata_fb, coldata_fb$Stage == 'CTRL')$Sample], 1, FUN=quantile)), 
                                    by='row.names')
rownames(all_proteins_rna_norm_fb) <- all_proteins_rna_norm_fb$Row.names
all_proteins_rna_norm_fb <- all_proteins_rna_norm_fb[,-c(1,2)]
colnames(all_proteins_rna_norm_fb)[c((ncol(all_proteins_rna_norm_fb)-4):ncol(all_proteins_rna_norm_fb))] <- c("RNA FB - 0% quartile Control", "RNA FB - 25% quartile Control", "RNA FB - 50% quartile Control", "RNA FB - 75% quartile Control", "RNA FB - 100% quartile Control")

all_proteins_rna_norm_fb <- merge(all_proteins_rna_norm_fb[,c('gene_name', "RNA FB - 0% quartile Control", "RNA FB - 25% quartile Control", 
                                                              "RNA FB - 50% quartile Control", "RNA FB - 75% quartile Control", "RNA FB - 100% quartile Control"), 
                                                           drop=F], 
                                  t(apply(all_proteins_rna_norm_fb[, subset(coldata_fb, coldata_fb$Stage == 'HD')$Sample], 1, FUN=quantile)), 
                                  by='row.names')
rownames(all_proteins_rna_norm_fb) <- all_proteins_rna_norm_fb$Row.names
all_proteins_rna_norm_fb <- all_proteins_rna_norm_fb[,-c(1)]
colnames(all_proteins_rna_norm_fb)[c((ncol(all_proteins_rna_norm_fb)-4):ncol(all_proteins_rna_norm_fb))] <- c("RNA FB - 0% quartile HD", "RNA FB - 25% quartile HD", "RNA FB - 50% quartile HD", "RNA FB - 75% quartile HD", "RNA FB - 100% quartile HD")

all_proteins <- merge(all_proteins_rna_norm_fb, all_proteins, by.y='Gene',by.x='gene_name', all.y=T)


all_proteins_rna_norm_in <- subset(rna_norm_in, rna_norm_in$gene_name %in% all_proteins$gene_name)
rownames(all_proteins_rna_norm_in) <- make.unique(all_proteins_rna_norm_in$gene_name)

all_proteins_rna_norm_in <- merge(all_proteins_rna_norm_in, 
                                  t(apply(all_proteins_rna_norm_in[, subset(coldata_in, coldata_in$Stage == 'CTRL')$Sample], 1, FUN=quantile)), 
                                  by='row.names')
rownames(all_proteins_rna_norm_in) <- all_proteins_rna_norm_in$Row.names
all_proteins_rna_norm_in <- all_proteins_rna_norm_in[,-c(1,2)]
colnames(all_proteins_rna_norm_in)[c((ncol(all_proteins_rna_norm_in)-4):ncol(all_proteins_rna_norm_in))] <- c("RNA iN - 0% quartile Control", "RNA iN - 25% quartile Control", "RNA iN - 50% quartile Control", "RNA iN - 75% quartile Control", "RNA iN - 100% quartile Control")

all_proteins_rna_norm_in <- merge(all_proteins_rna_norm_in[,c('gene_name', "RNA iN - 0% quartile Control", "RNA iN - 25% quartile Control", 
                                                              "RNA iN - 50% quartile Control", "RNA iN - 75% quartile Control", "RNA iN - 100% quartile Control"), 
                                                           drop=F], 
                                  t(apply(all_proteins_rna_norm_in[, subset(coldata_in, coldata_in$Stage == 'HD')$Sample], 1, FUN=quantile)), 
                                  by='row.names')
rownames(all_proteins_rna_norm_in) <- all_proteins_rna_norm_in$Row.names
all_proteins_rna_norm_in <- all_proteins_rna_norm_in[,-c(1)]
colnames(all_proteins_rna_norm_in)[c((ncol(all_proteins_rna_norm_in)-4):ncol(all_proteins_rna_norm_in))] <- c("RNA iN - 0% quartile HD", "RNA iN - 25% quartile HD", "RNA iN - 50% quartile HD", "RNA iN - 75% quartile HD", "RNA iN - 100% quartile HD")

all_proteins <- merge(all_proteins_rna_norm_in, all_proteins, by='gene_name', all.y=T)
colnames(all_proteins)[1] <- 'Gene'

write.xlsx(all_proteins, file='/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/allproteins_alldatasets.xlsx', row.names = F)

# Protein ----
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
# protein_test <- protein_test[which(protein_test$`Shapiro Pvalue - HD` > 0.05 & protein_test$`Shapiro Pvalue - Control` > 0.05),]
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

png('/Volumes/Seagate Backup /hd_in/07.20_hd/plots/protein_fb.in_meanplot.png', res = 200, height = 15, width = 15, units = "cm")
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

# Overlap with RNAseq ----
# UPSET 
input <- data.frame(gene_name=unique(c(rna_norm_in_test_upreg$gene_name,
                                       rna_norm_fb_test_upreg$gene_name,
                                       protein_in_test_upreg$Gene,
                                       protein_fb_test_upreg$Gene,
                                       rna_norm_in_test_dwnreg$gene_name,
                                       rna_norm_fb_test_dwnreg$gene_name,
                                       protein_in_test_dwnreg$Gene,
                                       protein_fb_test_dwnreg$Gene)))


input$RNA.Upreg.iN <- ifelse(input$gene_name %in% rna_norm_in_test_upreg$gene_name, 1, 0)
input$RNA.Upreg.FB <- ifelse(input$gene_name %in% rna_norm_fb_test_upreg$gene_name, 1, 0)
input$Protein.Upreg.iN <- ifelse(input$gene_name %in% protein_in_test_upreg$Gene, 1, 0)
input$Protein.Upreg.FB <- ifelse(input$gene_name %in% protein_fb_test_upreg$Gene, 1, 0)
input$RNA.Downreg.iN <- ifelse(input$gene_name %in% rna_norm_in_test_dwnreg$gene_name, 1, 0)
input$RNA.Downreg.FB <- ifelse(input$gene_name %in% rna_norm_fb_test_dwnreg$gene_name, 1, 0)
input$Protein.Downreg.iN <- ifelse(input$gene_name %in% protein_in_test_dwnreg$Gene, 1, 0)
input$Protein.Downreg.FB <- ifelse(input$gene_name %in% protein_fb_test_dwnreg$Gene, 1, 0)

library(UpSetR)
p <- upset(input, sets = c("RNA.Downreg.iN", # #417096
                           "RNA.Downreg.FB", # #417096
                           "Protein.Downreg.iN", # #9dc3e3
                           "Protein.Downreg.FB", # #9dc3e3
                           "RNA.Upreg.iN", # #d15685
                           "RNA.Upreg.FB", # #d15685
                           "Protein.Upreg.iN", # #f2a7c4
                           "Protein.Upreg.FB"), # #f2a7c4
           mainbar.y.label = "Gene Intersections", sets.x.label = "Genes Per Group", 
           keep.order=TRUE, 
           sets.bar.color=c('#417096', '#417096', '#9dc3e3', '#9dc3e3', '#d15685', '#d15685', '#f2a7c4', '#f2a7c4'))

p <- upset(input, sets = c("RNA.Downreg.iN", # #417096
                           # "RNA.Downreg.FB", # #417096
                           "Protein.Downreg.iN", # #9dc3e3
                           # "Protein.Downreg.FB"), # #9dc3e3
                           "RNA.Upreg.iN", # #d15685
                           # "RNA.Upreg.FB", # #d15685
                           "Protein.Upreg.iN"), #f2a7c4
                           #"Protein.Upreg.FB"), # #f2a7c4
           mainbar.y.label = "Gene Intersections", sets.x.label = "Genes Per Group", 
           keep.order=TRUE, 
           sets.bar.color=c('#417096', '#417096', '#9dc3e3', '#9dc3e3'))



pdf('/Volumes/Seagate Backup /hd_in/07.20_hd/plots/upset_rna_protein_iN.pdf', width = 7, height = 5)
print(p)
dev.off()  

# Jiovanis protein coding genes: Significant differences between control iNs and fibroblasts ----
in_fb <- c("BECN1",
           "CAMKK2",
           "PRKAA1",
           "EEF2K",
           "IRS1",
           "PPP2R5E",
           "PPP2R1B", 
           "ACTB")
previous_analysis_fb <- melt(subset(protein_fb, protein_fb$gene_name %in% in_fb), by='gene_name')
previous_analysis_fb$type <- rep("FB", nrow(previous_analysis_fb))
previous_analysis_in <- melt(subset(protein_in, protein_in$gene_name %in% in_fb), by='gene_name')
previous_analysis_in$type <- rep("IN", nrow(previous_analysis_in))

previous_analysis <- rbind(previous_analysis_fb, previous_analysis_in)
previous_analysis$variable <- as.character(previous_analysis$variable)
previous_analysis <- merge(previous_analysis, coldata[,c(1,2)], by.x='variable', by.y='Sample')
previous_analysis$stage_type <- paste(previous_analysis$type, previous_analysis$Stage, sep="_")
previous_analysis$stage_type <- factor(previous_analysis$stage_type, levels = c("IN_HD", "IN_CTRL", "FB_HD", "FB_CTRL"))
str(previous_analysis)
library(ggpubr)
comparisons <- list(c("IN_HD", "IN_CTRL"), c("FB_HD", "FB_CTRL"))
p <- ggplot(previous_analysis, aes(x=stage_type, y=value, fill=Stage)) + geom_boxplot() + 
  facet_wrap( ~ gene_name, scales="free") + theme_classic() + stat_compare_means(comparisons = comparisons, method = 't.test', exact=F)

pdf('/Volumes/Seagate Backup /hd_in/07.20_hd/plots/some_pathways_representatives.pdf', width=10, height=10)
print(p)
dev.off()

rownames(input) <- input$gene_name

# Pathway analysis STRING ----
string_db <- fread('/Volumes/Seagate Backup /hd_in/09.19_hd/9606.protein.info.v11.0.txt', data.table = F)
background <- length(unique(string_db$preferred_name))

pathway_downreg_bp <- fread('/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/downregulated/STRING/barplot/in_downreg_enrichment.bp.tsv', data.table = F)
pathway_downreg_bp$log10fdr <- -log10(pathway_downreg_bp$`false discovery rate`)
pathway_downreg_bp$expected <- ((pathway_downreg_bp$`background gene count`) / background)*length(protein_in_test_dwnreg$Gene)
pathway_downreg_bp$FC <- pathway_downreg_bp$`observed gene count` / pathway_downreg_bp$expected
pathway_downreg_bp$`term description` <- factor(pathway_downreg_bp$`term description`, levels=pathway_downreg_bp[order(pathway_downreg_bp$FC), 'term description'])

write.xlsx(pathway_downreg_bp, '/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/downregulated/STRING/barplot/in_downreg_enrichment.bp.xlsx', row.names = F)

pathway_downreg_bp_p <- ggplot(data=pathway_downreg_bp, aes(x=`term description` , y=FC)) + geom_bar(stat='identity') +
  geom_point(aes(y=log10fdr), size=3) + coord_flip() + theme_classic() + geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red') + 
  labs(y='Fold Change') + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "-log10(padj)"))

pdf('/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/downregulated/STRING/barplot/in_downreg_enrichment.bp.pdf', height = 10, width = 10)
print(pathway_downreg_bp_p)
dev.off()

pathway_downreg_reactome <- fread('/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/downregulated/STRING/barplot/in_downreg_enrichment.RCTM.tsv', data.table = F)
pathway_downreg_reactome$log10fdr <- -log10(pathway_downreg_reactome$`false discovery rate`)
pathway_downreg_reactome$expected <- ((pathway_downreg_reactome$`background gene count`) / background)*length(protein_in_test_dwnreg$Gene)
pathway_downreg_reactome$FC <- pathway_downreg_reactome$`observed gene count` / pathway_downreg_reactome$expected
pathway_downreg_reactome$`term description` <- factor(pathway_downreg_reactome$`term description`, levels=pathway_downreg_reactome[order(pathway_downreg_reactome$FC), 'term description'])

write.xlsx(pathway_downreg_reactome, '/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/downregulated/STRING/barplot/in_downreg_enrichment.RCTM.xlsx', row.names = F)

pathway_downreg_reactome_p <- ggplot(data=pathway_downreg_reactome, aes(x=`term description` , y=FC)) + geom_bar(stat='identity') +
  geom_point(aes(y=log10fdr), size=3) + coord_flip() + theme_classic() + geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red') + 
  labs(y='Fold Change') + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "-log10(padj)"))

pdf('/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/downregulated/STRING/barplot/in_downreg_enrichment.RCTM.pdf', height = 10, width = 10)
print(pathway_downreg_reactome_p)
dev.off()

pathway_downreg_bp <- fread('/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis_and_STRING/downregulated/STRING/barplot/in_downreg_enrichment.BiologicalProcess.tsv', data.table = F)
pathway_downreg_bp$log10fdr <- -log10(pathway_downreg_bp$`false discovery rate`)
pathway_downreg_bp$expected <- ((pathway_downreg_bp$`background gene count`) / background)*length(protein_in_test_dwnreg$Gene)
pathway_downreg_bp$FC <- pathway_downreg_bp$`observed gene count` / pathway_downreg_bp$expected
pathway_downreg_bp$`term description` <- factor(pathway_downreg_bp$`term description`, levels=pathway_downreg_bp[order(pathway_downreg_bp$FC), 'term description'])

write.xlsx(pathway_downreg_bp, '/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis_and_STRING/downregulated/STRING/barplot/in_downreg_enrichment.bp.xlsx', row.names = F)

pathway_downreg_bp_p <- ggplot(data=pathway_downreg_bp, aes(x=`term description` , y=FC)) + geom_bar(stat='identity') +
  geom_point(aes(y=log10fdr), size=3) + coord_flip() + theme_classic() + geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red') + 
  labs(y='Fold Change') + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "-log10(padj)"))

pdf('/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis_and_STRING/downregulated/STRING/barplot/in_downreg_enrichment.bp.pdf', height = 15, width = 10)
print(pathway_downreg_bp_p)
dev.off()


pathway_upreg_bp <- fread('/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis_and_STRING/upregulated/STRING/barplot/in_upreg_enrichment.bp.tsv', data.table = F)
pathway_upreg_bp$log10fdr <- -log10(pathway_upreg_bp$`false discovery rate`)
pathway_upreg_bp$expected <- ((pathway_upreg_bp$`background gene count`) / background)*length(protein_in_test_dwnreg$Gene)
pathway_upreg_bp$FC <- pathway_upreg_bp$`observed gene count` / pathway_upreg_bp$expected
pathway_upreg_bp$`term description` <- factor(pathway_upreg_bp$`term description`, levels=pathway_upreg_bp[order(pathway_upreg_bp$FC), 'term description'])

write.xlsx(pathway_upreg_bp, '/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis_and_STRING/upregulated/STRING/barplot/in_upreg_enrichment.bp.xlsx', row.names = F)

pathway_upreg_bp_p <- ggplot(data=pathway_upreg_bp, aes(x=`term description` , y=FC)) + geom_bar(stat='identity') +
  geom_point(aes(y=log10fdr), size=3) + coord_flip() + theme_classic() + geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red') + 
  labs(y='Fold Change') + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "-log10(padj)"))

pdf('/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis_and_STRING/upregulated/STRING/barplot/in_upreg_enrichment.bp.pdf', height = 10, width = 10)
print(pathway_upreg_bp_p)
dev.off()

pathway_upreg_reactome <- fread('/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/upregulated/STRING/barplot/in_upreg_enrichment.RCTM.tsv', data.table = F)
pathway_upreg_reactome$log10fdr <- -log10(pathway_upreg_reactome$`false discovery rate`)
pathway_upreg_reactome$expected <- ((pathway_upreg_reactome$`background gene count`) / background)*length(protein_in_test_dwnreg$Gene)
pathway_upreg_reactome$FC <- pathway_upreg_reactome$`observed gene count` / pathway_upreg_reactome$expected
pathway_upreg_reactome$`term description` <- factor(pathway_upreg_reactome$`term description`, levels=pathway_upreg_reactome[order(pathway_upreg_reactome$FC), 'term description'])

write.xlsx(pathway_upreg_reactome, '/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis/upregulated/STRING/barplot/in_upreg_enrichment.RCTM.xlsx', row.names = F)

pathway_upreg_reactome_p <- ggplot(data=pathway_upreg_reactome, aes(x=`term description` , y=FC)) + geom_bar(stat='identity') +
  geom_point(aes(y=log10fdr), size=3) + coord_flip() + theme_classic() + geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red') + 
  labs(y='Fold Change') + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "-log10(padj)"))

pdf('/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis_and_STRING/upregulated/STRING/barplot/in_upreg_enrichment.RCTM.pdf', height = 10, width = 10)
print(pathway_upreg_reactome_p)
dev.off()

pathway_upreg_bp <- fread('/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis_and_STRING/upregulated/STRING/barplot/in_upreg_enrichment.BiologicalProcess.tsv', data.table = F)
pathway_upreg_bp$log10fdr <- -log10(pathway_upreg_bp$`false discovery rate`)
pathway_upreg_bp$expected <- ((pathway_upreg_bp$`background gene count`) / background)*length(protein_in_test_dwnreg$Gene)
pathway_upreg_bp$FC <- pathway_upreg_bp$`observed gene count` / pathway_upreg_bp$expected
pathway_upreg_bp$`term description` <- factor(pathway_upreg_bp$`term description`, levels=pathway_upreg_bp[order(pathway_upreg_bp$FC), 'term description'])

write.xlsx(pathway_upreg_bp, '/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis_and_STRING/upregulated/STRING/barplot/in_upreg_enrichment.bp.xlsx', row.names = F)

pathway_upreg_bp_p <- ggplot(data=pathway_upreg_bp, aes(x=`term description` , y=FC)) + geom_bar(stat='identity') +
  geom_point(aes(y=log10fdr), size=3) + coord_flip() + theme_classic() + geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red') + 
  labs(y='Fold Change') + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "-log10(padj)"))

pdf('/Volumes/Seagate Backup /hd_in/07.20_hd/5_proteomics/GO_analysis_and_STRING/upregulated/STRING/barplot/in_upreg_enrichment.bp.pdf', height = 10, width = 10)
print(pathway_upreg_bp_p)
dev.off()

