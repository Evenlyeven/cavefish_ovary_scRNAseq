# Reproductive adaptation of Astyanax mexicanus under nutrient limitation
Xia, F., Santacruz, A., Wu, D., Bertho, S., Fritz, E., Morales-Sosa, P., McKinney, S., Nowotarski, S. H., & Rohner, N.
<br>

## :desktop_computer:  R script for scRNAseq analysis
Refer to [scRNAseq_analysis.rmd](https://github.com/Evenlyeven/cavefish_ovary_scRNAseq/blob/main/scRNAseq_analysis.rmd) or the rendered version [scRNAseq_analysis.html](https://github.com/Evenlyeven/cavefish_ovary_scRNAseq/blob/main/scRNAseq_analysis.html) for the analysis code associated with the scRNA-seq data. 

## :test_tube:  Raw data for scRNAseq

The corresponding raw data, including processed seurat object, has been deposited in GEO under accession number [GSE282933](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE282933).

## :sparkles:  Shiny App

The interactive Shiny app for exploring the data is available [here](https://simrcompbio.shinyapps.io/astmex_ovary_scrnaseq/). The script used to build the Shiny app can be accessed [here](https://github.com/Evenlyeven/cavefish_ovary_scRNAseq/blob/main/app.R).

## :dna: RBHs between cavefish and zebrafish

[rbh_am_dr.csv](https://github.com/Evenlyeven/cavefish_ovary_scRNAseq/blob/main/rbh_am_dr.csv) contains the reciprocal best hits (RBHs) between cavefish and zebrafish. RBHs were identified from tblastx results between GRCz11 (Ens110 annotation) and Astyanax_mexicanus-2.0 (Ens110 annotation) by filtering for E-value < 0.001, bit score > 100, and alignment length > 40, then selecting gene pairs with the highest bit scores as best hits and retaining mutual best hits between species.
