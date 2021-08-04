library(clusterProfiler)
library(msigdbr)
set.seed(10)
# GSEA ----
m_tfs <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
  dplyr::select(gs_name, human_gene_symbol)#, human_gene_symbol)

## feature 1: numeric vector
geneList <- rna_norm_in_test$Log2FC

## feature 2: named vector
names(geneList) <- as.character(rna_norm_in_test$gene_name)

## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)
em2 <- GSEA(geneList, TERM2GENE = m_tfs, seed = T)
BiocManager::install("enrichplot", force=T)
library(enrichplot)
# write.table(em2@result, sep='\t', file = "~/Desktop/gsea.tsv", quote = F)
gseaplot2(em2, geneSetID = c("ZNF23_TARGET_GENES", "ZNF16_TARGET_GENES"), 
          title = "Enrichment of ZNF23 and ZNF16 target genes")

genes_genesets_upreg <- unlist(str_split(paste(em2@result[which(em2@result$enrichmentScore > 0),"core_enrichment"], collapse = "/"), "/"))

gseaplot2(em2, geneSetID = c("ATGTCAC_MIR489", "MIR508_3P", "MIR433_3P"), 
          title = "Enrichment of ZNF23 and ZNF16 target genes")
# ZNF16 --> GATA-1 PPAR-gamma1 PPAR-gamma2
# ZNF23 --> AREB6 COUP COUP-TF COUP-TF1 HFH-1 HNF-4alpha1 HNF-4alpha2 LUN-1 PPAR-gamma1 PPAR-gamma2
tmp <- protein
rownames(tmp) <- tmp$gene_unique

# Genes belonging to the gene sets found to be enriched in the downregulated genes ----
genes_genesets_dwnreg <- unlist(str_split(paste(em2@result[which(em2@result$enrichmentScore < 0),"core_enrichment"], collapse = "/"), "/"))
genes_genesets_dwnreg_df <- merge(data.frame(gene_name = genes_genesets_dwnreg), unique(gene_transcript[,c(3,4)]), by='gene_name')
tmp1 <- genes_genesets_dwnreg_df[which(genes_genesets_dwnreg_df$gene_type == "protein_coding"),]
tmp2 <- genes_genesets_dwnreg_df[which(!genes_genesets_dwnreg_df[which(genes_genesets_dwnreg_df$gene_type != "protein_coding"),"gene_name"] %in% tmp1$gene_name),]
tmp2 <- tmp2[which(!duplicated(tmp2$gene_name)),]
# How many of those were protein coding genes?
table(rbind(tmp1, tmp2)$gene_type)

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
# Without filtering ----
length(genes_genesets_dwnreg[which(genes_genesets_dwnreg %in% protein_in_test$`Gene name`)])

genes_genesets_dwnreg_test <- protein_in_test[which(protein_in_test$`Gene name` %in% genes_genesets_dwnreg[which(genes_genesets_dwnreg %in% protein_in_test$`Gene name`)]),]
ggplot(genes_genesets_dwnreg_test, aes(x=1, y=Log2FC)) + 
  geom_boxplot() + geom_hline(yintercept = 0, linetype = "dashed", colour="red") +
  theme_classic() + theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()) + ggtitle("Log2FC of protein targets of enriched TF among\ndownregulated genes in RNAseq")

# What is the mean of these log2FC?
# mean(genes_genesets_dwnreg_test[which(!is.na(genes_genesets_dwnreg_test$Log2FC)), "Log2FC"])
# Now with the upregulated ----
genes_genesets_upreg_df <- tmp[genes_genesets_upreg[which(genes_genesets_upreg %in% tmp$Gene)],]
genes_genesets_upreg_df <- merge(data.frame(gene_name = genes_genesets_upreg), unique(gene_transcript[,c(3,4)]), by='gene_name')
genes_genesets_upreg_df[which(!duplicated(genes_genesets_upreg_df$gene_name)),]

tmp1 <- unique(genes_genesets_upreg_df[which(genes_genesets_upreg_df$gene_type == "protein_coding"),])
tmp2 <- unique(genes_genesets_upreg_df[which(!genes_genesets_upreg_df[which(genes_genesets_upreg_df$gene_type != "protein_coding"),"gene_name"] %in% tmp1$gene_name),])
tmp2 <- tmp2[which(!duplicated(tmp2$gene_name)),]
table(rbind(tmp1, tmp2)$gene_type)
 
# How many were captured in our data?
length(genes_genesets_upreg[which(genes_genesets_upreg %in% protein_in_test$`Gene name`)])

genes_genesets_upreg_test <- protein_in_test[which(protein_in_test$`Gene name` %in% genes_genesets_upreg[which(genes_genesets_upreg %in% protein_in_test$`Gene name`)]),]
ggplot(genes_genesets_upreg_test, aes(x=1, y=Log2FC)) + 
  geom_boxplot() + geom_hline(yintercept = 0, linetype = "dashed", colour="red") +
  theme_classic() + theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()) + ggtitle("Log2FC of protein targets of enriched TF among\nupregulated genes in RNAseq")


# c("OCT1", "POU2F1", "TAF9", "TAF9B", "ZNF581") %in% rownames(tmp)

taf9_proteomics <- as.data.frame(t(tmp[c("TAF9"),coldata_in$Sample]))
taf9_proteomics <- merge(taf9_proteomics, coldata_in, by='row.names')
library(ggpubr)
ggplot(taf9_proteomics, aes(x=Stage, y=TAF9, fill=Stage)) + geom_boxplot() + theme_classic() +
  labs(x="Condition", fill = "Condition", y="Protein levels") + ggtitle("TAF9 (iN proteomics)")


# Heterogeneity of samples? ----
top10_upreg_protein_in <- head(protein_in_test_upreg[order(protein_in_test_upreg$Log2FC, decreasing = T),], 10)$Gene
top10_dwnreg_protein_in <- head(protein_in_test_dwnreg[order(protein_in_test_dwnreg$Log2FC, decreasing = T),], 10)$Gene

# Top 10 upregulated proteins in RNAseq (iNs)
coldata_in <- coldata_in[order(coldata_in$Stage),]
top10_upreg_protein_in_df <- rna_counts_norm[top10_upreg_protein_in,]
top10_upreg_protein_in_df <- reshape2::melt(top10_upreg_protein_in_df[,-1], by="gene_name")
top10_upreg_protein_in_df <- merge(top10_upreg_protein_in_df, coldata_in, by.x='variable', by.y="Sample")

ggplot(top10_upreg_protein_in_df, aes(x=Stage, y=value, fill=Stage)) + geom_boxplot() + 
  facet_wrap(.~gene_name, scales = "free_y") + theme_classic() + 
  ggtitle("Top 10 upregulated proteins in RNAseq (iNs)") + 
  labs(y="Normalized expression (median of ratios)", fill="Condition", x="Condition")

pheatmap(log2(rna_counts_norm[top10_upreg_protein_in, coldata_in$Sample]+0.5), 
         cluster_cols = F, gaps_col = 7, main = "Top 10 upregulated proteins in RNAseq (iNs)")

# Top 10 downregulated proteins in RNAseq (iNs)
top10_dwnreg_protein_in_df <- rna_counts_norm[top10_dwnreg_protein_in,]
top10_dwnreg_protein_in_df <- reshape2::melt(top10_dwnreg_protein_in_df[,-1], by="gene_name")
top10_dwnreg_protein_in_df <- merge(top10_dwnreg_protein_in_df, coldata_in, by.x='variable', by.y="Sample")

ggplot(top10_dwnreg_protein_in_df, aes(x=Stage, y=value, fill=Stage)) + geom_boxplot() + 
  facet_wrap(.~gene_name, scales = "free_y") + theme_classic() + 
  ggtitle("Top 10 downregulated proteins in RNAseq (iNs)") + 
  labs(y="Normalized expression (median of ratios)", fill="Condition", x="Condition")

pheatmap(log2(rna_counts_norm[top10_dwnreg_protein_in, coldata_in$Sample]+0.5), 
         cluster_cols = F, gaps_col = 7, main = "Top 10 downregulated proteins in RNAseq (iNs)")


# Top 10 downregulated genes in proteomics (iNs)
top10_upreg_rna_in <- head(rna_norm_in_test_upreg[order(rna_norm_in_test_upreg$Log2FC, decreasing = T),], 10)$Gene
top10_dwnreg_rna_in <- head(rna_norm_in_test_dwnreg[order(rna_norm_in_test_dwnreg$Log2FC, decreasing = T),], 10)$Gene

table(top10_upreg_rna_in %in% rownames(tmp))
top10_dwnreg_rna_in <- top10_dwnreg_rna_in[which(top10_dwnreg_rna_in %in% rownames(tmp))]

top10_dwnreg_rna_in_df <- tmp[top10_dwnreg_rna_in,]
top10_dwnreg_rna_in_df <- reshape2::melt(top10_dwnreg_rna_in_df[,-1], by="gene_unique")
top10_dwnreg_rna_in_df <- merge(top10_dwnreg_rna_in_df, coldata_in, by.x='variable', by.y="Sample")

ggplot(top10_dwnreg_rna_in_df, aes(x=Stage, y=value, fill=Stage)) + geom_boxplot() + 
  facet_wrap(.~gene_unique, scales = "free_y") + theme_classic() + 
  ggtitle("Top 10 downregulated genes in proteomics (iNs)") + 
  labs(y="Protein levels", fill="Condition", x="Condition")

pheatmap(log2(tmp[top10_dwnreg_rna_in, coldata_in$Sample]+0.5), 
         cluster_rows = F, cluster_cols = F, gaps_col = 7, main = "Top 10 downregulated genes in proteomics (iNs)",
         labels_row = top10_dwnreg_rna_in)

# Cortical vs striatal ----
striatal <- toupper(unique(c("Drd1", "Ido1", "Ppp1r1b", "Adora2a", "Ido1", "Adora2a", "Gpr83", "Penk", "Drd1", "Dclk3", "Mhrt", "Drd1", "Lrpprc", "Tac1", "Wnt2", "A230065H16Rik", "Zic1", "2810459M11Rik", "Gm14964")))
striatal <- striatal[which(striatal %in% tmp$Gene)]

cortical <- toupper(unique(c("Myl4", "Cpne4", "Rprm", "Rell1", "Cplx3", "Nxph3", "Hs3st4", "Sulf1", "Prss12", "Igfbp6", "C130074G19Rik", "Nr4a2", "Col24a1", "Oprk1", "Hs3st2", "Vipr1", "Pde1a", "Dkkl1", "Galntl6", "Krt12", "Tcap", "A830009L08Rik", "Gm12371", "Lamp5", "Tshz2", "Ddit4l", "Wfs1", "RP24-134N2.1", "Vwc2l", "Cbln1", "Ntsr1", "Rxfp1", "Fermt1", "RP23-231J2.1", "Dbpht2", "Cdkl4", "Scube1", "Tekt5", "Trim54", "Slc30a3", "Lmo3", "Abi3bp", "4930426D05Rik", "Ndst4", "Klhl14", "Rgs14", "Oxtr", "Bmp3", "Fezf2", "RP24-134N2.1", "Kcng1", "Pthlh", "Cox6a2", "Rbp4", "Cox6a2", "Cort", "Sst", "Ccna1", "Crhbp", "Tacr1", "Chodl", "Lhx6", "Pde11a", "Npy", "Ptchd2", "Cplx3", "Teddm3", "Krt73", "Stk32b", "5330429C05Rik", "Yjefn3", "Hdhd3", "Npas1", "Col19a1", "Slc17a8", "Vip", "Gm17750", "Cxcl14", "Vip", "Tiam1", "Tac2", "Epha7")))
cortical <- cortical[which(cortical %in% tmp$Gene)]

markers <- data.frame(markers = c(cortical,striatal),
           type = c(rep("cortical", length(cortical)), rep("striatal", length(striatal))))
rownames(markers) <- markers$markers


pheatmap(log2(tmp[markers$markers,coldata_in$Sample]+0.5), 
         cluster_rows = F, 
         cluster_cols = F,
         annotation_row = markers[,"type", drop=F],
         annotation_col = coldata_in[,"Stage", drop=F],
         gaps_col = 7, fontsize = 8, show_colnames = F)

