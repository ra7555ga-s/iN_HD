Here we present the methods to test for differential expression analysis in RNAseq and proteomics that was reported in "Distinct sub-cellular autophagy impairments occur independently of protein aggregation in aged induced neurons from patients with Huntington’s disease" (Pircs, et al., not yet published).

We have organized the papers results into eight sections, each of them with their corresponding R script. Please note that these scripts are dependent to one another, so they save their environment and the next one loads the previous one, they are described (in order) below: 

1. Figure 1d presents a mean plot with the differentially expressed genes (pvalue < 0.05 in an unpaired t-tests). The test and the scatter plot can be generated with `RNA_in_vs_fb.R` which will generate `RNA_in_vs_fb.RData`.
2. Differential gene expression in fibroblasts' RNAseq between Huntington's patients and controls can be found in `fb_RNA_hd_vs_ctrl.R`. This script loads `RNA_in_vs_fb.RData` and generates `fb_RNA_hd_vs_ctrl.RData`. 
3. Differential gene expression in iNs' RNAseq between Huntington's patients and controls can be found in `in_RNA_hd_vs_ctrl.R`. This will generate the scatter plot shown in Figure 2b. This script loads `fb_RNA_hd_vs_ctrl.RData` and generates `in_RNA_hd_vs_ctrl.RData`.
4. For the RNA GO analysis, we performed an overrepresentation test using the Slim Biological Process Database of PANTHER. The visualization for the GO analysis results from points 1, 2, and 3 were done using `RNA_GOanalysis.R`. This script loads `in_RNA_hd_vs_ctrl.RData` and generates `RNA_GOanalysis.RData`.
5. Figure 1h presents a mean plot with the differentially expressed proteins (pvalue < 0.05 in an unpaired t-tests). The test and the scatter plot can be generated with `protein_in_vs_fb.R`. This script loads `RNA_GOanalysis.RData` and generates `protein_in_vs_fb.RData`. 
6. Differential protein expression in fibroblasts' RNAseq between Huntington's patients and controls can be found in `fb_protein_hd_vs_ctrl.R`. This script loads `protein_in_vs_fb.RData` and generates `fb_protein_hd_vs_ctrl.RData`. 
7. Differential protein expression in iNs' RNAseq between Huntington's patients and controls can be found in `in_protein_hd_vs_ctrl.R`. This will generate the scatter plot shown in Figure 2c. This script loads `fb_protein_hd_vs_ctrl.RData` and generates `in_protein_hd_vs_ctrl.RData`.
8. Upset plot to visualize the intersections between the following groups:
    * Fibroblasts proteomics upregulated in HD
    * iNs proteomics upregulated in HD
    * Fibroblasts RNAseq upregulated in HD
    * iNs RNAseq upregulated in HD
    * Fibroblasts proteomics downregulated in HD
    * iNs proteomics downregulated in HD
    * Fibroblasts RNAseq downregulated in HD
    * iNs RNAseq downregulated in HD
Can be generated using `RNA_protein_upset_plot.R`. This script loads `in_protein_hd_vs_ctrl.RData`.