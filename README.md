# PoritesWPS
global gene expression patterns in WPS affected tissues

Files in this repository 
-----------

1. Analyses_final.R: R script for plotting physiological data and generating linear mixed models 
	- Input file: Patchy_Porites_Feb2019_for_R.csv
	
2. DESeq_PoritesPatch.R: R script for conducting differential gene expression analysis and KOG term enrichments
	- Input file: AllCountsHost.txt
	- Input file: plob_iso2kogClassNR.tab
	- Input file: amil_iso2kogClassNR.tab
	- Input file: MetaAnalysisFiles.RData
	- Input file: plob_iso2gene.tab
	
	- Output file: hostVSDandPVALS_no_g4_deseq1_4jun_plusPC1.csv
	- Output file: GOpatchHost.csv
	- Output file: GObinaryHostPatch.csv
	- Output file: VSDs_GObinaryHostPatch.csv

3. For gene ontology enrichment scripts and example input files, please visit https://github.com/ckenkel/GO_MWU
