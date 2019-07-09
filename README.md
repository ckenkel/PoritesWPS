# PoritesWPS
global gene expression patterns in WPS affected tissues

Files in this repository 
-----------

1. Put all this into the same directory:
	- scripts: GO_MWU.R, gomwu_a.pl, gomwu_b.pl, gomwu.functions.R
	- GO hierarchy file 
		(version 1.0, http://www.geneontology.org/GO.downloads.ontology.shtml)
	- table of GO annotations for your sequences: two-column (gene id - GO terms), 
		tab-delimited, one line per gene, multiple GO terms separated by semicolon. 
		If you have multiple lines per gene, use nrify_GOtable.pl to merge them.
	- table of measure of interest for your sequences: two columns of comma-separated 
		values: gene id, continuous measure of significance such as log(fold-change) or
		-log(p-value). To perform standard GO enrichment analysis based on Fisher's 
		exact test, use binary measure (1 or 0, i.e., either sgnificant or not). To analyze modules derived from WGCNA, specify 0 for genes not included in the module and the kME value (number between 0 and 1, module membership score) for genes included in the module.
	
	It is important to have the latter two tables representing the whole 
	genome (or transcriptome) - at least the portion that was measured -
	rather than some select group of genes since the test relies on comparing
	the behavior of individual GO categories to the whole.

2. For completing gene ontology enrichments, please visit 
