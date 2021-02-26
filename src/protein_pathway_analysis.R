setwd('/Volumes/My Passport/hd_in/24.02.20/')
load('RNA_protein_upset_plot.RData')

# Pathway analysis STRING ----
string_db <- fread('/Volumes/My Passport/hd_in/09.19_hd/9606.protein.info.v11.0.txt', data.table = F)
background <- length(unique(string_db$preferred_name))

in_pathway_downreg_bp <- fread('/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/downregulated/in_protein_hd.ctrl_dwnreg_enrichment.BiologicalProcess.tsv', data.table = F)
in_pathway_downreg_bp$log10fdr <- -log10(in_pathway_downreg_bp$`false discovery rate`)
# Expected is the number of genes you would expect in your list for this category, based on the reference list.
# We calculate the ratio of the genes that belong to that term (using all genes annotated in the STRING database as reference)
# We multiply by the number of genes we have on the list
# We get how many genes we would expect if we would randomly sample the database
in_pathway_downreg_bp$expected <- ((in_pathway_downreg_bp$`background gene count`) / background)*length(protein_in_test_dwnreg$Gene)
in_pathway_downreg_bp$FC <- in_pathway_downreg_bp$`observed gene count` / in_pathway_downreg_bp$expected
in_pathway_downreg_bp$`term description` <- factor(in_pathway_downreg_bp$`term description`, levels=in_pathway_downreg_bp[order(in_pathway_downreg_bp$FC), 'term description'])

write.xlsx(in_pathway_downreg_bp, '/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/downregulated/in_protein_hd.ctrl_dwnreg_enrichment.BiologicalProcess.xlsx', row.names = F)

in_pathway_downreg_bp_p <- ggplot(data=in_pathway_downreg_bp, aes(x=`term description` , y=FC)) + geom_bar(stat='identity') +
  geom_point(aes(y=log10fdr), size=3) + coord_flip() + theme_classic() + geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red') + 
  labs(y='Fold Change') + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "-log10(padj)"))


svg(paste(getwd(), "/plots/in_pathway_downreg_bp.svg", sep=''))
print(in_pathway_downreg_bp_p)
dev.off()

pdf(paste(getwd(), "/plots/in_pathway_downreg_bp.pdf", sep=''))
print(in_pathway_downreg_bp_p)
dev.off()

fb_pathway_downreg_bp <- fread('/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/downregulated/fb_protein_hd.ctrl_dwnreg_enrichment.BiologicalProcess.tsv', data.table = F)
fb_pathway_downreg_bp$log10fdr <- -log10(fb_pathway_downreg_bp$`false discovery rate`)
fb_pathway_downreg_bp$expected <- ((fb_pathway_downreg_bp$`background gene count`) / background)*length(protein_fb_test_dwnreg$Gene)
fb_pathway_downreg_bp$FC <- fb_pathway_downreg_bp$`observed gene count` / fb_pathway_downreg_bp$expected
fb_pathway_downreg_bp$`term description` <- factor(fb_pathway_downreg_bp$`term description`, levels=fb_pathway_downreg_bp[order(fb_pathway_downreg_bp$FC), 'term description'])

write.xlsx(fb_pathway_downreg_bp, '/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/downregulated/fb_protein_hd.ctrl_dwnreg_enrichment.BiologicalProcess.xlsx', row.names = F)

fb_pathway_downreg_bp_p <- ggplot(data=fb_pathway_downreg_bp, aes(x=`term description` , y=FC)) + geom_bar(stat='identity') +
  geom_point(aes(y=log10fdr), size=3) + coord_flip() + theme_classic() + geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red') + 
  labs(y='Fold Change') + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "-log10(padj)"))

in_pathway_upreg_bp <- fread('/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/upregulated/in_protein_hd.ctrl_upreg_enrichment.BiologicalProcess.tsv', data.table = F)
in_pathway_upreg_bp$log10fdr <- -log10(in_pathway_upreg_bp$`false discovery rate`)
in_pathway_upreg_bp$expected <- ((in_pathway_upreg_bp$`background gene count`) / background)*length(protein_in_test_upreg$Gene)
in_pathway_upreg_bp$FC <- in_pathway_upreg_bp$`observed gene count` / in_pathway_upreg_bp$expected
in_pathway_upreg_bp$`term description` <- factor(in_pathway_upreg_bp$`term description`, levels=in_pathway_upreg_bp[order(in_pathway_upreg_bp$FC), 'term description'])

write.xlsx(in_pathway_upreg_bp, '/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/upregulated/in_protein_hd.ctrl_upreg_enrichment.BiologicalProcess.xlsx', row.names = F)

in_pathway_upreg_bp_p <- ggplot(data=in_pathway_upreg_bp, aes(x=`term description` , y=FC)) + geom_bar(stat='identity') +
  geom_point(aes(y=log10fdr), size=3) + coord_flip() + theme_classic() + geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red') + 
  labs(y='Fold Change') + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "-log10(padj)"))

fb_pathway_upreg_bp <- fread('/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/upregulated/fb_protein_hd.ctrl_upreg_enrichment.BiologicalProcess.tsv', data.table = F)
fb_pathway_upreg_bp$log10fdr <- -log10(fb_pathway_upreg_bp$`false discovery rate`)
fb_pathway_upreg_bp$expected <- ((fb_pathway_upreg_bp$`background gene count`) / background)*length(protein_fb_test_upreg$Gene)
fb_pathway_upreg_bp$FC <- fb_pathway_upreg_bp$`observed gene count` / fb_pathway_upreg_bp$expected
fb_pathway_upreg_bp$`term description` <- factor(fb_pathway_upreg_bp$`term description`, levels=unique(fb_pathway_upreg_bp[order(fb_pathway_upreg_bp$FC), 'term description']))

write.xlsx(fb_pathway_upreg_bp, '/Volumes/My Passport/hd_in/24.02.20/5_proteomics/GO_analysis/upregulated/fb_protein_hd.ctrl_upreg_enrichment.BiologicalProcess.xlsx', row.names = F)

fb_pathway_upreg_bp_p <- ggplot(data=fb_pathway_upreg_bp, aes(x=`term description` , y=FC)) + geom_bar(stat='identity') +
  geom_point(aes(y=log10fdr), size=3) + coord_flip() + theme_classic() + geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red') + 
  labs(y='Fold Change') + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "-log10(padj)"))


save.image('protein_pathway_analysis.RData')

