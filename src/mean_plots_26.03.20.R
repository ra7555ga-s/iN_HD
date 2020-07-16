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
  
  ggsave(plot, file=paste(prefix, 'fdr.svg', sep = ''), width=32, height=length(results_for_plot$GO_label)/1.2, units="cm")
  
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
rna <- rna[,which(colnames(rna) %in% coldata$Sample)]

coldata_rna <- subset(coldata, coldata$Sample %in% colnames(rna))

# RNAseq FBs ----
coldata_rna_fb <- subset(coldata, coldata$Sample %in% colnames(rna) & coldata$CellType == 'FB')

# Filter out genes with less than 10 reads in one sample
rna_expressed_fb <- rna[which(apply(rna[,coldata_rna_fb$Sample], 1, more_10) > 0), coldata_rna_fb$Sample]

# DESeq DEA and normalizing with median of ratios
dds_fb <- DESeqDataSetFromMatrix(rna_expressed_fb[,coldata_rna_fb$Sample], coldata_rna_fb, design= ~Stage)
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
  
  row <- row[coldata_rna_fb$Sample,,drop=F]
  row <- merge(row, coldata_rna_fb[,2,drop=F], by='row.names')
  row[,2] <- as.numeric(as.character(row[,2]))
  
  mean_ctrl <- mean(subset(row, row$Stage == 'CTRL')[,2])
  mean_hd <- mean(subset(row, row$Stage == 'HD')[,2])

  gene_id <- colnames(row)[2]
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
    
    results <- c(gene_id, mean_ctrl, mean_hd, log2((mean_hd/mean_ctrl)), normality_ctrl$p.value, normality_hd$p.value, test$statistic, test$p.value, test$conf.int)
      
  }
  else{
    results <- c(gene_id, mean_ctrl, mean_hd, log2((mean_hd/mean_ctrl)), NA, NA, NA, NA, NA, NA)
  }

  names(results) <- c("Gene id", "Mean Control", "Mean HD", "Log2FC", 'Shapiro Pvalue - Control', 'Shapiro Pvalue - HD', 'T', 'Pvalue', 'Low conf int', 'High conf int')
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
rna_norm_fb_test$Padj <- p.adjust(rna_norm_fb_test$Pvalue, 'fdr')
# Filter out genes that are not normally distributed
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

pdf('/Volumes/Seagate Backup /hd_in/09.19_hd/plots/rna_fb_meanplot.pdf')
plot(log2(rna_norm_fb_test$`Mean Control`+0.5), 
     log2(rna_norm_fb_test$`Mean HD`+0.5), 
     col=rna_norm_fb_test$colours, 
     cex=rna_norm_fb_test$cexs, 
     pch=16, 
     xlab='log2(mean Control)', 
     ylab='log2(mean HD)')

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
            file='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/upregulated/fb_rna_upreg.tab')

rna_norm_fb_test_dwnreg <- merge(rna_norm_fb_test_signdiff[which(rna_norm_fb_test_signdiff$type == 'Downregulated'),], unique(gene_transcript[,c(2,3)]), by.x='row.names', by.y='gene_id')
write.table(rna_norm_fb_test_dwnreg$gene_name, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/downregulated/fb_rna_dwnreg.tab')

rna_norm_fb_test <- merge(rna_norm_fb_test, unique(gene_transcript[,c(2,3)]), by.x='Gene id', by.y='gene_id')
write.table(rna_norm_fb_test$gene_name, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/fb_rna_expressed.tab')


# GO analysis FB RNA ----
# Biological process
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/upregulated/slim_biological_process/fb_slim_biological_process_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/upregulated/slim_biological_process/fb_slim_biological_process_upreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 15)
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/downregulated/slim_biological_process/fb_slim_biological_process_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/downregulated/slim_biological_process/fb_slim_biological_process_dwnreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 15)
# Molecular function
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/upregulated/slim_molecular_function/fb_slim_molecular_function_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/upregulated/slim_molecular_function/fb_slim_molecular_function_upreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 15)
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/downregulated/slim_molecular_function/fb_slim_molecular_function_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/downregulated/slim_molecular_function/fb_slim_molecular_function_dwnreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth = 15)

# RNAseq iNs ----
coldata_rna_in <- subset(coldata, coldata$Sample %in% colnames(rna) & coldata$CellType == 'IN')

# DESeq DEA and normalizing with median of ratios
dds_in <- DESeqDataSetFromMatrix(rna[,coldata_rna_in$Sample], coldata_rna_in, design= ~Stage)
dds_in$Stage <- relevel(dds_in$Stage, 'CTRL')
dds_in <- DESeq(dds_in)
res_in <- results(dds_in)

rna_norm_in <- as.data.frame(counts(dds_in, normalize=T))
rna_norm_in$gene_id <- rownames(rna_norm_in)

rna_norm_in_test <- apply(rna_norm_in, 1, FUN=function(row){
  row <- as.data.frame(row)
  colnames(row) <- row['gene_id',]
  # print(row)
  row <- row[coldata_rna_in$Sample,,drop=F]
  row <- merge(row, coldata_rna_in[,2,drop=F], by='row.names')
  row[,2] <- as.numeric(as.character(row[,2]))
  
  mean_ctrl <- mean(subset(row, row$Stage == 'CTRL')[,2])
  mean_hd <- mean(subset(row, row$Stage == 'HD')[,2])
  
  gene_id <- colnames(row)[2]
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

pdf('/Volumes/Seagate Backup /hd_in/09.19_hd/plots/rna_in_meanplot.pdf')
plot(log2(rna_norm_in_test$`Mean Control` + 0.5), 
     log2(rna_norm_in_test$`Mean HD` + 0.5), 
     col=rna_norm_in_test$colours, 
     cex=rna_norm_in_test$cexs, 
     pch=16, 
     xlab='log2(mean Control)', 
     ylab='log2(mean HD)')

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
            file='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/upregulated/in_rna_upreg.tab')

rna_norm_in_test_dwnreg <- merge(rna_norm_in_test[which(rna_norm_in_test$type == 'Downregulated'),], unique(gene_transcript[,c(2,3)]), by.x='row.names', by.y='gene_id')
write.table(rna_norm_in_test_dwnreg$gene_name, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/downregulated/in_rna_dwnreg.tab')

rna_norm_in_test <- merge(rna_norm_in_test, unique(gene_transcript[,c(2,3)]), by.x='Gene id', by.y='gene_id')
write.table(rna_norm_in_test$gene_name, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/in_rna_expressed.tab')


# GO analysis iN ----
# Biological process
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/upregulated/slim_biological_process/in_slim_biological_process_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/upregulated/slim_biological_process/in_slim_biological_process_upreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 10)
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/downregulated/slim_biological_process/in_slim_biological_process_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/downregulated/slim_biological_process/in_slim_biological_process_dwnreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 25, pdf_heigth= 15)
# Molecular function
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/upregulated/slim_molecular_function/in_slim_molecular_function_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/upregulated/slim_molecular_function/in_slim_molecular_function_upreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 15)
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/downregulated/slim_molecular_function/in_slim_molecular_function_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/3_stdmapping/GO_analysis/downregulated/slim_molecular_function/in_slim_molecular_function_dwnreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth = 15)

# Protein ----
protein <- fread('/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/20200309_HDProteomics_Log2_SubtractMedianColumn_AverageSample_+20.txt', data.table=F)
protein <- protein[,c('Gene', colnames(protein)[which(colnames(protein) %in% colnames(rna))])]

protein_fb <- protein[, c(colnames(protein)[which(colnames(protein) %in% subset(coldata, coldata$CellType == 'FB')$Sample)], 'Gene')]
colnames(protein_fb)[ncol(protein_fb)] <- 'gene_name'
protein_in <- protein[, c(colnames(protein)[which(colnames(protein) %in% subset(coldata, coldata$CellType == 'IN')$Sample)], 'Gene')]
colnames(protein_in)[ncol(protein_in)] <- 'gene_name'

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
  colnames(row) <- row['gene_name',]
  
  row <- row[coldata_rna_fb$Sample,,drop=F]
  row <- merge(row, coldata_rna_fb[,2,drop=F], by='row.names')
  row[,2] <- as.numeric(as.character(row[,2]))
  
  mean_ctrl <- mean(subset(row, row$Stage == 'CTRL')[,2])
  mean_hd <- mean(subset(row, row$Stage == 'HD')[,2])
  
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
    
    results <- c(gene_name, mean_ctrl, mean_hd, log2((mean_hd/mean_ctrl)), normality_ctrl$p.value, normality_hd$p.value, test$statistic, test$p.value, test$conf.int)
    
  }
  else{
    results <- c(gene_name, mean_ctrl, mean_hd, log2((mean_hd/mean_ctrl)), NA, NA, NA, NA, NA, NA)
  }
  
  names(results) <- c("Gene name", "Mean Control", "Mean HD", "Log2FC", 'Shapiro Pvalue - Control', 'Shapiro Pvalue - HD', 'T', 'Pvalue', 'Low conf int', 'High conf int')
  
  return(results)
})

protein_fb_test <- as.data.frame(t(protein_fb_test))
protein_fb_test$`Gene name` <- as.character(protein_fb_test$`Gene name`)
protein_fb_test$`Mean Control` <- as.numeric(as.character(protein_fb_test$`Mean Control`))
protein_fb_test$`Mean HD` <- as.numeric(as.character(protein_fb_test$`Mean HD`))
protein_fb_test$`Log2FC` <- as.numeric(as.character(protein_fb_test$`Log2FC`))
protein_fb_test$`Shapiro Pvalue - HD` <- as.numeric(as.character(protein_fb_test$`Shapiro Pvalue - HD`))
protein_fb_test$`Shapiro Pvalue - Control` <- as.numeric(as.character(protein_fb_test$`Shapiro Pvalue - Control`))
protein_fb_test$T <- as.numeric(as.character(protein_fb_test$T))
# T= Mean difference 
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

protein_fb_test$`log2 Mean Control` <- log2(protein_fb_test$`Mean Control` + 0.5) 
protein_fb_test$`log2 Mean HD` <- log2(protein_fb_test$`Mean HD` + 0.5)

pdf('/Volumes/Seagate Backup /hd_in/09.19_hd/plots/protein_fb_meanplot.pdf')
plot(protein_fb_test$`Mean Control`, 
     protein_fb_test$`Mean HD`, 
     col=protein_fb_test$colours, 
     cex=protein_fb_test$cexs, 
     pch=16, 
     xlab='mean of normalized expression in Control', 
     ylab='mean of normalized expression in HD')

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
write.table(protein_fb_test_upreg$`Gene name`, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/upregulated/fb_protein_upreg.tab')

protein_fb_test_dwnreg <- protein_fb_test_signdiff[which(protein_fb_test_signdiff$type == 'Downregulated'),]
write.table(protein_fb_test_dwnreg$`Gene name`, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/downregulated/fb_protein_dwnreg.tab')

write.table(protein_fb_test$`Gene name`, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/fb_protein_expressed.tab')

# GO analysis FB protein ----
# Biological process
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/upregulated/slim_biological_process/fb_slim_biological_process_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/upregulated/slim_biological_process/fb_slim_biological_process_upreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth=15)
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/downregulated/slim_biological_process/fb_slim_biological_process_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/downregulated/slim_biological_process/fb_slim_biological_process_dwnreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth=15)
# Molecular function
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/upregulated/slim_molecular_function/fb_slim_molecular_function_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/upregulated/slim_molecular_function/fb_slim_molecular_function_upreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 15, pdf_heigth=15)
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/downregulated/slim_molecular_function/fb_slim_molecular_function_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/downregulated/slim_molecular_function/fb_slim_molecular_function_dwnreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 15, pdf_heigth=15)

# Protein iN ----
protein_in_test <- apply(protein_in, 1, FUN=function(row){
  row <- as.data.frame(row)
  colnames(row) <- row['gene_name',]
  
  row <- row[coldata_rna_in$Sample,,drop=F]
  row <- merge(row, coldata_rna_in[,2,drop=F], by='row.names')
  row[,2] <- as.numeric(as.character(row[,2]))
  
  mean_ctrl <- mean(subset(row, row$Stage == 'CTRL')[,2])
  mean_hd <- mean(subset(row, row$Stage == 'HD')[,2])
  
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
    
    results <- c(gene_name, mean_ctrl, mean_hd, log2((mean_hd/mean_ctrl)), normality_ctrl$p.value, normality_hd$p.value, test$statistic, test$p.value, test$conf.int)
    
  }
  else{
    results <- c(gene_name, mean_ctrl, mean_hd, log2((mean_hd/mean_ctrl)), NA, NA, NA, NA, NA, NA)
  }
  
  names(results) <- c("Gene name", "Mean Control", "Mean HD", "Log2FC", 'Shapiro Pvalue - Control', 'Shapiro Pvalue - HD', 'T', 'Pvalue', 'Low conf int', 'High conf int')
  
  return(results)
})
protein_in_test <- as.data.frame(t(protein_in_test))
protein_in_test$`Gene name` <- as.character(protein_in_test$`Gene name`)
protein_in_test$`Mean Control` <- as.numeric(as.character(protein_in_test$`Mean Control`))
protein_in_test$`Mean HD` <- as.numeric(as.character(protein_in_test$`Mean HD`))
protein_in_test$`Log2FC` <- as.numeric(as.character(protein_in_test$`Log2FC`))
protein_in_test$`Shapiro Pvalue - HD` <- as.numeric(as.character(protein_in_test$`Shapiro Pvalue - HD`))
protein_in_test$`Shapiro Pvalue - Control` <- as.numeric(as.character(protein_in_test$`Shapiro Pvalue - Control`))
protein_in_test$T <- as.numeric(as.character(protein_in_test$T))
protein_in_test$Pvalue <- as.numeric(as.character(protein_in_test$Pvalue))
protein_in_test$`Low conf int` <- as.numeric(as.character(protein_in_test$`Low conf int`))
protein_in_test$`High conf int` <- as.numeric(as.character(protein_in_test$`High conf int`))

table(protein_in_test$`Shapiro Pvalue - Control` > 0.05)
table(protein_in_test$`Shapiro Pvalue - HD` > 0.05)
table(protein_in_test$`Shapiro Pvalue - HD` > 0.05 | protein_in_test$`Shapiro Pvalue - Control` > 0.05)
table(protein_in_test$`Shapiro Pvalue - HD` > 0.05 & protein_in_test$`Shapiro Pvalue - Control` > 0.05)

# Filter out genes that are not normally distributed
protein_in_test <- protein_in_test[which(protein_in_test$`Shapiro Pvalue - HD` > 0.05 & protein_in_test$`Shapiro Pvalue - Control` > 0.05),]
# Get the counts of the ones that are normally distributed
protein_in <- subset(protein_in, protein_in$gene_name %in% protein_in_test$`Gene name`)
# Calculate the log2(mean per condition + 0.5)

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

# protein_in_test[which(protein_in_test$`Mean Control` == min(protein_in_test$`Mean Control`)),]
# row <- t(protein_in[which(protein_in$gene_name == 'SEMG1'),])
# tmp <- row[nrow(row),]
# row <- as.data.frame(row[coldata_rna_in$Sample,, drop=F])
# row$sample <- rownames(row)
# colnames(row)[1] <- 'expression'
# row <- merge(row, coldata_rna_in, by.x='sample', by.y='Sample')
# row$expression <- as.numeric(as.character(row$expression))
# shapiro.test(subset(row, row$Stage == 'HD')$expression)
# shapiro.test(row$expression)

pdf('/Volumes/Seagate Backup /hd_in/09.19_hd/plots/protein_in_meanplot.pdf')
plot(protein_in_test$`Mean Control`, 
     protein_in_test$`Mean HD`, 
     col=protein_in_test$colours, 
     cex=protein_in_test$cexs, 
     pch=16, 
     xlab='mean of normalized expression in Control', 
     ylab='mean of normalized expression in HD')#, ylim=c(14,30), xlim=c(14,30))

legend("bottomright", legend = c(paste("up (",as.numeric(table(protein_in_test$type)["Upregulated"]),")",sep=""),
                                 paste("down (",as.numeric(table(protein_in_test$type)["Downregulated"]),")",sep = ""),
                                 paste("not significant (",as.numeric(table(protein_in_test$type)["Not significant"]),")",sep = "")),
       pch=16,col=c("firebrick3","steelblue4","black"),cex=1)

dev.off()

# Up and down regulated for PANTHER
protein_in_test_upreg <- protein_in_test[which(protein_in_test$type == 'Upregulated'),]
write.table(protein_in_test_upreg$`Gene name`, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/upregulated/in_protein_upreg.tab')

protein_in_test_dwnreg <- protein_in_test[which(protein_in_test$type == 'Downregulated'),]
write.table(protein_in_test_dwnreg$`Gene name`, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/downregulated/in_protein_dwnreg.tab')

write.table(protein_in_test$`Gene name`, quote=F, col.names = F, row.names = F, 
            file='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/in_protein_expressed.tab')

# GO analysis iN protein ----
# Biological process
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/upregulated/slim_biological_process/in_slim_biological_process_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/upregulated/slim_biological_process/in_slim_biological_process_upreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 15, pdf_heigth = 15)
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/downregulated/slim_biological_process/in_slim_biological_process_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/downregulated/slim_biological_process/in_slim_biological_process_dwnreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 15, pdf_heigth = 15)
# Molecular function
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/upregulated/slim_molecular_function/in_slim_molecular_function_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/upregulated/slim_molecular_function/in_slim_molecular_function_upreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 15, pdf_heigth = 15)
GO_json(file_name='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/downregulated/slim_molecular_function/in_slim_molecular_function_vs2.json', prefix='/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/GO_analysis/downregulated/slim_molecular_function/in_slim_molecular_function_dwnreg_vs2', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 17, pdf_heigth = 12)


# Overlap with RNAseq ----
# Protein data has gene names as rownames
rna_norm_in_gene_names <- merge(rna_norm_in, unique(gene_transcript[,c(2,3)]), by='gene_id')
rna_norm_in_protein_subset <- rna_norm_in_gene_names[which(rna_norm_in_gene_names$gene_name %in% protein_in$gene_name),]
# rownames(rna_norm_in_protein_subset) <- make.unique(rna_norm_in_protein_subset$gene_name)
# Remove gene name and gene id column
# rna_norm_in_protein_subset <- rna_norm_in_protein_subset[,coldata_rna_in$Sample]

# Just the samples we have in the RNAseq
# rownames(protein_in) <- make.unique(protein_in$gene_name)
protein_in_rna_subset <- protein_in[, c(which(colnames(protein_in) %in% colnames(rna_norm_in_protein_subset)))]

protein_means_in <- data.frame(gene_name=protein_in_rna_subset[,'gene_name'],
                               ctrl=rowMeans(protein_in_rna_subset[,which(colnames(protein_in_rna_subset) %in% subset(coldata_rna_in, coldata_rna_in$Stage == 'CTRL')$Sample)]),
                               hd=rowMeans(protein_in_rna_subset[,which(colnames(protein_in_rna_subset) %in% subset(coldata_rna_in, coldata_rna_in$Stage == 'HD')$Sample)]))

# rna_norm_in_protein_subset <- log2(rna_norm_in_protein_subset[,coldata_rna_in$Sample]+0.5)
rna_means_in_gene_names <- data.frame(gene_name=rna_norm_in_protein_subset[, 'gene_name'],
                                      ctrl=rowMeans(rna_norm_in_protein_subset[,which(colnames(rna_norm_in_protein_subset) %in% subset(coldata_rna_in, coldata_rna_in$Stage == 'CTRL')$Sample)]),
                                      hd=rowMeans(rna_norm_in_protein_subset[,which(colnames(rna_norm_in_protein_subset) %in% subset(coldata_rna_in, coldata_rna_in$Stage == 'HD')$Sample)]))

protein_means_in_ctrl <- protein_means_in[,c("gene_name", 'ctrl')]
rna_means_in_ctrl <- rna_means_in_gene_names[,c("gene_name", 'ctrl')]

protein_vs_rna_means_in_ctrl <- merge(rna_means_in_ctrl, protein_means_in_ctrl, all.y=T, by='gene_name')
colnames(protein_vs_rna_means_in_ctrl) <- c('gene_name', 'rna', 'protein')
protein_vs_rna_means_in_ctrl$rna <- as.numeric(as.character(protein_vs_rna_means_in_ctrl$rna))
protein_vs_rna_means_in_ctrl$protein <- as.numeric(as.character(protein_vs_rna_means_in_ctrl$protein))

protein_vs_rna_means_in_hd <- merge(protein_means_in[,c('gene_name', 'hd')], rna_means_in_gene_names[,c('gene_name', 'hd')], all.x=T, by='gene_name')
colnames(protein_vs_rna_means_in_hd) <- c('gene_name', 'protein', 'rna')

protein_vs_rna_means_in_hd <- protein_vs_rna_means_in_hd[which(!is.na(protein_vs_rna_means_in_hd$rna)),]

protein_vs_rna_means_in_hd$rna <- as.numeric(as.character(protein_vs_rna_means_in_hd$rna))
protein_vs_rna_means_in_hd$protein <- as.numeric(as.character(protein_vs_rna_means_in_hd$protein))

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

protein_vs_rna_means_in_hd$log2rna <- log2(protein_vs_rna_means_in_hd$rna+0.5)
in_hd_lm <- lm(data=protein_vs_rna_means_in_hd, protein~log2rna)
p <- ggplotRegression(in_hd_lm, 'iN HD', xlab='log2(mean RNA expression)', ylab='Normalized mean protein expression')
ggsave(plot=p, '/Volumes/Seagate Backup /hd_in/09.19_hd/plots/in_hd_protein_vs_rna.pdf')

protein_vs_rna_means_in_ctrl$log2rna <- log2(protein_vs_rna_means_in_ctrl$rna+0.5)
in_ctrl_lm <- lm(data=protein_vs_rna_means_in_ctrl, protein~log2rna)
p <- ggplotRegression(in_ctrl_lm, 'iN Control', xlab='log2(mean RNA expression)', ylab='Normalized mean protein expression')
ggsave(plot=p, '/Volumes/Seagate Backup /hd_in/09.19_hd/plots/in_control_protein_vs_rna.pdf')

protein_vs_rna_means_in_ctrl$condition <- rep('Control', nrow(protein_vs_rna_means_in_ctrl))
protein_vs_rna_means_in_hd$condition <- rep('HD', nrow(protein_vs_rna_means_in_hd))
# protein_vs_rna_means_in_hd$gene_name <- rownames(protein_vs_rna_means_in_hd)
# protein_vs_rna_means_in_hd <- protein_vs_rna_means_in_hd[,colnames(protein_vs_rna_means_in_ctrl)]

protein_vs_rna_means_in <- rbind(protein_vs_rna_means_in_hd, protein_vs_rna_means_in_ctrl)
protein_vs_rna_means_in$condition <- factor(protein_vs_rna_means_in$condition, levels=c('Control', 'HD'))
in_lm <- lm(data=protein_vs_rna_means_in, protein~log2rna+condition)
ggplotRegression(in_lm, 'iN : protein~log2(rna)+condition', xlab='log2(mean RNA expression)', ylab='Normalized mean protein expression')

summary(in_lm)
summary(in_ctrl_lm)
summary(in_ko_lm)

# UPSET ----
input <- data.frame(gene_name=unique(c(rna_norm_in_test_upreg$gene_name,
                                       rna_norm_fb_test_upreg$gene_name,
                                       protein_in_test_upreg$`Gene name`,
                                       protein_fb_test_upreg$`Gene name`,
                                       rna_norm_in_test_dwnreg$gene_name,
                                       rna_norm_fb_test_dwnreg$gene_name,
                                       protein_in_test_dwnreg$`Gene name`,
                                       protein_fb_test_dwnreg$`Gene name`)))


input$RNA.Upreg.iN <- ifelse(input$gene_name %in% rna_norm_in_test_upreg$gene_name, 1, 0)
input$RNA.Upreg.FB <- ifelse(input$gene_name %in% rna_norm_fb_test_upreg$gene_name, 1, 0)
input$Protein.Upreg.iN <- ifelse(input$gene_name %in% protein_in_test_upreg$`Gene name`, 1, 0)
input$Protein.Upreg.FB <- ifelse(input$gene_name %in% protein_fb_test_upreg$`Gene name`, 1, 0)
input$RNA.Downreg.iN <- ifelse(input$gene_name %in% rna_norm_in_test_dwnreg$gene_name, 1, 0)
input$RNA.Downreg.FB <- ifelse(input$gene_name %in% rna_norm_fb_test_dwnreg$gene_name, 1, 0)
input$Protein.Downreg.iN <- ifelse(input$gene_name %in% protein_in_test_dwnreg$`Gene name`, 1, 0)
input$Protein.Downreg.FB <- ifelse(input$gene_name %in% protein_fb_test_dwnreg$`Gene name`, 1, 0)

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

p

pdf('/Volumes/Seagate Backup /hd_in/09.19_hd/plots/upset_rna_protein.pdf', width = 7, height = 5)
print(p)
dev.off()  

# Pathways representatives - Protein ----
protein <- fread('/Volumes/Seagate Backup /hd_in/09.19_hd/5_proteomics/20200309_HDProteomics_Log2_SubtractMedianColumn_AverageSample_+20.txt', data.table=F)
protein <- protein[,c('Gene', colnames(protein)[which(colnames(protein) %in% colnames(rna))])]

protein_fb <- protein[, c(colnames(protein)[which(colnames(protein) %in% subset(coldata, coldata$CellType == 'FB')$Sample)], 'Gene')]
colnames(protein_fb)[ncol(protein_fb)] <- 'gene_name'
protein_in <- protein[, c(colnames(protein)[which(colnames(protein) %in% subset(coldata, coldata$CellType == 'IN')$Sample)], 'Gene')]
colnames(protein_in)[ncol(protein_in)] <- 'gene_name'

protein_in[is.na(protein_in)] <- 0
protein_fb[is.na(protein_fb)] <- 0

 in_fb <- c("BECN1",
           "CAMKK2",
           "PRKAA1",
           "EEF2K",
           "IRS1",
           "PPP2R5E",
           "PPP2R1B")#,
           # "ACTBL2")
previous_analysis_fb <- melt(subset(protein_fb, protein_fb$gene_name %in% in_fb), by='gene_name')
previous_analysis_fb$type <- rep("FB", nrow(previous_analysis_fb))
previous_analysis_in <- melt(subset(protein_in, protein_in$gene_name %in% in_fb), by='gene_name')
previous_analysis_in$type <- rep("IN", nrow(previous_analysis_in))

previous_analysis <- rbind(previous_analysis_fb, previous_analysis_in)
previous_analysis$variable <- as.character(previous_analysis$variable)
previous_analysis <- merge(previous_analysis, coldata_rna[,c(1,2)], by.x='variable', by.y='Sample')
previous_analysis$stage_type <- paste(previous_analysis$type, previous_analysis$Stage, sep="_")
previous_analysis$stage_type <- factor(previous_analysis$stage_type, levels = c("IN_HD", "IN_CTRL", "FB_HD", "FB_CTRL"))
str(previous_analysis)
library(ggpubr)
comparisons <- list(c("IN_HD", "IN_CTRL"), c("FB_HD", "FB_CTRL"))
p <- ggplot(previous_analysis, aes(x=stage_type, y=value, fill=Stage)) + geom_boxplot() + 
  facet_wrap( ~ gene_name, scales="free") + theme_classic() + stat_compare_means(comparisons = comparisons)

pdf('/Volumes/Seagate Backup /hd_in/09.19_hd/plots/some_pathways_representatives_protein.pdf', width=10, height=10)
# pdf('/Volumes/Seagate Backup /hd_in/09.19_hd/plots/actbl2_protein.pdf', width=10, height=10)
print(p)
dev.off()


# Pathways representatives - RNA ----
rna_norm_fb <- as.data.frame(counts(dds_fb, normalize=T))
rna_norm_fb$gene_id <- rownames(rna_norm_fb)
rna_norm_in <- as.data.frame(counts(dds_in, normalize=T))
rna_norm_in$gene_id <- rownames(rna_norm_in)
in_fb <- 'ACTBL2'
rna_norm_fb <- merge(rna_norm_fb, unique(gene_transcript[,c('gene_id', 'gene_name')]), by='gene_id')
rna_norm_in <- merge(rna_norm_in, unique(gene_transcript[,c('gene_id', 'gene_name')]), by='gene_id')
previous_analysis_fb <- melt(subset(rna_norm_fb, rna_norm_fb$gene_name %in% in_fb), by='gene_name')
previous_analysis_fb$type <- rep("FB", nrow(previous_analysis_fb))
previous_analysis_in <- melt(subset(rna_norm_in, rna_norm_in$gene_name %in% in_fb), by='gene_name')
previous_analysis_in$type <- rep("IN", nrow(previous_analysis_in))

previous_analysis <- rbind(previous_analysis_fb, previous_analysis_in)
previous_analysis$variable <- as.character(previous_analysis$variable)
previous_analysis <- merge(previous_analysis, coldata_rna[,c(1,2)], by.x='variable', by.y='Sample')
previous_analysis$stage_type <- paste(previous_analysis$type, previous_analysis$Stage, sep="_")
previous_analysis$stage_type <- factor(previous_analysis$stage_type, levels = c("IN_HD", "IN_CTRL", "FB_HD", "FB_CTRL"))
str(previous_analysis)
library(ggpubr)
comparisons <- list(c("IN_HD", "IN_CTRL"), c("FB_HD", "FB_CTRL"))
p <- ggplot(previous_analysis, aes(x=stage_type, y=value, fill=Stage)) + geom_boxplot() +
  facet_wrap( ~ gene_name, scales="free") + theme_classic() + stat_compare_means(comparisons = comparisons)

# compute lower and upper whiskers
ylim1 = boxplot.stats(previous_analysis$y)$stats[c(1, 5)]

# scale y limits based on ylim1
p1 = p0 + coord_cartesian(ylim = ylim1*1.05)

pdf('/Volumes/Seagate Backup /hd_in/09.19_hd/plots/actbl2_rna.pdf', width=10, height=10)
print(p)
dev.off()


# Venn diagrams ----
library(venneuler)
vd <- venneuler(c(A=0.3, B=0.3, C=1.1, "A&B"=0.1, "A&C"=0.2, "B&C"=0.1 ,"A&B&C"=0.1))
plot(vd)
