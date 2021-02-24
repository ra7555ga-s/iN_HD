setwd('/Volumes/My Passport/hd_in/24.02.20/')
load('in_protein_hd_vs_ctrl.RData')

# Overlap with RNAseq ----
to_file <- c("Gene", "50% quartile Control", "50% quartile HD", "Log2FC", "Shapiro Pvalue - Control",
             "Shapiro Pvalue - HD", "T", "Pvalue", "Low conf int", "High conf int", "type")
colnames(rna_norm_in_test_upreg) <- c("gene_id","Gene id","50% quartile Control","50% quartile HD",
                                      "Log2FC","Shapiro Pvalue - Control",
                                      "Shapiro Pvalue - HD","T","Pvalue","Low conf int","High conf int","type",
                                      "colours","cexs","Gene")
colnames(rna_norm_in_test_dwnreg) <- c("gene_id","Gene id","50% quartile Control","50% quartile HD",
                                       "Log2FC","Shapiro Pvalue - Control",
                                       "Shapiro Pvalue - HD","T","Pvalue","Low conf int","High conf int","type",
                                       "colours","cexs","Gene")
colnames(rna_norm_fb_test_upreg) <- c("gene_id","Gene id","50% quartile Control","50% quartile HD",
                                      "Log2FC","Shapiro Pvalue - Control",
                                      "Shapiro Pvalue - HD","T","Pvalue","Low conf int","High conf int","type",
                                      "colours","cexs","Gene")
colnames(rna_norm_fb_test_dwnreg) <- c("gene_id","Gene id","50% quartile Control","50% quartile HD",
                                       "Log2FC","Shapiro Pvalue - Control",
                                       "Shapiro Pvalue - HD","T","Pvalue","Low conf int","High conf int","type",
                                       "colours","cexs","Gene")

write.xlsx(x = rna_norm_in_test_upreg[,to_file], row.names = F,
           file = paste(getwd(), "/tables/protein_rna_signdiff_genes.xlsx", sep=''),
           sheetName = "iN (RNA) - Upreg in HD")
write.xlsx(x = rna_norm_in_test_dwnreg[,to_file],  row.names = F,
           file = paste(getwd(), "/tables/protein_rna_signdiff_genes.xlsx", sep=''),
           sheetName = "iN (RNA) - Downreg in HD", append=TRUE)
write.xlsx(x = rna_norm_fb_test_upreg[,to_file],  row.names = F,
           file = paste(getwd(), "/tables/protein_rna_signdiff_genes.xlsx", sep=''),
           sheetName = "FB (RNA) - Upreg in HD", append=TRUE)
write.xlsx(x = rna_norm_fb_test_dwnreg[,to_file],  row.names = F,
           file = paste(getwd(), "/tables/protein_rna_signdiff_genes.xlsx", sep=''),
           sheetName = "FB (RNA) - Downreg in HD", append=TRUE)
write.xlsx(x = protein_in_test_upreg[,to_file], row.names = F,
           file = paste(getwd(), "/tables/protein_rna_signdiff_genes.xlsx", sep=''),
           sheetName = "iN (Protein) - Upreg in HD", append=TRUE)
write.xlsx(x = protein_in_test_dwnreg[,to_file], row.names = F,
           file = paste(getwd(), "/tables/protein_rna_signdiff_genes.xlsx", sep=''),
           sheetName = "iN (Protein) - Downreg in HD", append=TRUE)
write.xlsx(x = protein_fb_test_upreg[,to_file], row.names = F,
           file = paste(getwd(), "/tables/protein_rna_signdiff_genes.xlsx", sep=''),
           sheetName = "FB (Protein) - Upreg in HD", append=TRUE)
write.xlsx(x = protein_fb_test_dwnreg[,to_file], row.names = F,
           file = paste(getwd(), "/tables/protein_rna_signdiff_genes.xlsx", sep = ''),
           sheetName = "FB (Protein) - Downreg in HD", append=TRUE)
# UPSET 
input <- data.frame(gene_name=unique(c(rna_norm_in_test_upreg$Gene,
                                       rna_norm_fb_test_upreg$Gene,
                                       protein_in_test_upreg$Gene,
                                       protein_fb_test_upreg$Gene,
                                       rna_norm_in_test_dwnreg$Gene,
                                       rna_norm_fb_test_dwnreg$Gene,
                                       protein_in_test_dwnreg$Gene,
                                       protein_fb_test_dwnreg$Gene)))

input$RNA.Upreg.iN <- ifelse(input$gene_name %in% rna_norm_in_test_upreg$Gene, 1, 0)
input$RNA.Upreg.FB <- ifelse(input$gene_name %in% rna_norm_fb_test_upreg$Gene, 1, 0)
input$Protein.Upreg.iN <- ifelse(input$gene_name %in% protein_in_test_upreg$Gene, 1, 0)
input$Protein.Upreg.FB <- ifelse(input$gene_name %in% protein_fb_test_upreg$Gene, 1, 0)
input$RNA.Downreg.iN <- ifelse(input$gene_name %in% rna_norm_in_test_dwnreg$Gene, 1, 0)
input$RNA.Downreg.FB <- ifelse(input$gene_name %in% rna_norm_fb_test_dwnreg$Gene, 1, 0)
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

