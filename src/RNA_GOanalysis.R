setwd('/Volumes/My Passport/hd_in/24.02.20/')
load('in_RNA_hd_vs_ctrl.RData')

# Settings ----
library(ggplot2)
library(RColorBrewer)
library(data.table)
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
    transcript_gene <- fread('/Volumes/My Passport/annotation/human/gencode/gencode.v30.annotation.transcr.gene.tab', data.table = FALSE, header=FALSE)
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

GO_json(file_name='/Volumes/My Passport/hd_in/24.02.20/3_stdmapping/GO_analysis/upregulated/slim_biological_process/in_rna_hd.ctrl_upreg_biological_process.json', prefix='/Volumes/My Passport/hd_in/24.02.20/3_stdmapping/GO_analysis/upregulated/in_rna_hd.ctrl_upreg_bp_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 10)
GO_json(file_name='/Volumes/My Passport/hd_in/24.02.20/3_stdmapping/GO_analysis/downregulated/slim_biological_process/in_rna_hd.ctrl_dwnreg_biological_process.json', prefix='/Volumes/My Passport/hd_in/24.02.20/3_stdmapping/GO_analysis/downregulated/in_rna_hd.ctrl_dwnreg_bp_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 15, pdf_heigth= 15)

GO_json(file_name='/Volumes/My Passport/hd_in/24.02.20/3_stdmapping/GO_analysis/upregulated/slim_biological_process/fb_rna_hd.ctrl_upreg_biological_process.json', prefix='/Volumes/My Passport/hd_in/24.02.20/3_stdmapping/GO_analysis/upregulated/fb_rna_hd.ctrl_upreg_bp_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 10)
GO_json(file_name='/Volumes/My Passport/hd_in/24.02.20/3_stdmapping/GO_analysis/downregulated/slim_biological_process/fb_rna_hd.ctrl_dwnreg_biological_process.json', prefix='/Volumes/My Passport/hd_in/24.02.20/3_stdmapping/GO_analysis/downregulated/fb_rna_hd.ctrl_dwnreg_bp_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 15, pdf_heigth= 15)


# GO_json(file_name='/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/upregulated/fb_protein_hd.ctrl_upreg_biological_process.json', prefix='/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/upregulated/fb_protein_hd.ctrl_upreg_bp_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 15)
# GO_json(file_name='/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/upregulated/in_protein_hd.ctrl_upreg_biological_process.json', prefix='/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/upregulated/in_protein_hd.ctrl_upreg_bp_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth = 15)
# 
# GO_json(file_name='/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/downregulated/fb_protein_hd.ctrl_dwnreg_biological_process.json', prefix='/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/downregulated/fb_protein_hd.ctrl_dwnreg_bp_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 15)
# GO_json(file_name='/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/downregulated/in_protein_hd.ctrl_dwnreg_biological_process.json', prefix='/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/downregulated/in_protein_hd.ctrl_dwnreg_bp_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth = 15)
# 
# GO_json(file_name='/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/upregulated/protein_fb.in_upreg_biological_process.json', prefix='/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/upregulated/protein_fb.in_upreg_bp_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth= 15)
# GO_json(file_name='/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/downregulated/protein_fb.in_dwnreg_biological_process.json', prefix='/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/downregulated/protein_fb.in_dwnreg_bp_', morethan = 2, botlog2fc = -0.5, toplog2fc = 0.5, pdf_width = 20, pdf_heigth = 15)


save.image('RNA_GOanalysis.RData')
