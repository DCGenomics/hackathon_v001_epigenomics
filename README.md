# Creating Testable and Integrative Epigenetic models of Transcription

### Prerequites:
"There is an abundance of epigenetic data encompassing a wide variety of markers over many different cell types.
This surfeit of information empowers labs to infer how these different markers contribute to gene expression and 
chromatin state. However, there is currently no existing methodology to elucidate relationships between epigenetic 
modifiers.  We sought to rectify this by developing an analytic pipeline that efficiently  allows labs to form 
models of relationships between RNA-seq data, ChIP-seq data, and DNA methylation data from their own cell lines. 
These models can then be used to predict changes in gene expression with respect to changes in these epigenetic 
signals. This pipeline is also flexible in the sense that labs need not provide all three types of sequence reads. 
Publicly available datasets can be utilized to generate a model, which a lab can then use to predict the state of 
the chromatin based on their own epigenetic data. The pipeline uses a combination of python, R, and command line 
based tools." 


Prerequites:
R version >3.0 (Packages: data.table, plyr, reshape2, glmnet, preprocessCore, GenomicRanges) 
Python version >2.7 but not 3.0 

### Inputs:
2+ sets (patient/replicates) of bigWig files. Inputs are given in an input text file of format (including a header, **tab separated**):

[patient/cell name] [assay type] [absolute file location]  [aggregation function]

### Example:
Patient	| DataType | DataFile	| Aggregation
--------|----------|------------|------------
Patient1 | RNA | /data/patient1.RNA.bigWig | mean
Patient2 | RNA | /data/patient2.RNA.bigWig | mean
Patient3 | RNA | /data/patient3.RNA.bigWig | mean
Patient1 | H3K4me1	| /data/patient1.H3K4me1.bigWig | sum
Patient2 | H3K4me1	| /data/patient2.H3K4me1.bigWig | sum
Patient3 | H3K4me1	| /data/patient3.H3K4me1.bigWig | sum
Patient1 | H3K4me3	| /data/patient1.H3K4me3.bigWig | sum
Patient2 | H3K4me3	| /data/patient2.H3K4me3.bigWig | sum
Patient3 | H3K4me3	| /data/patient3.H3K4me3.bigWig | sum
Patient1 | RRBS	| /data/patient1.RRBS.bigWig | mean
Patient2 | RRBS	| /data/patient2.RRBS.bigWig | mean
Patient3 | RRBS	| /data/patient3.RRBS.bigWig | mean

Aggregation values: mean, mean0, sum, max, min


Each patient MUST have a row with DataType of "RNA". There should be at least one more row for each patient (with common DataType across all the patients).


**Current call:**
```
python pipeline_top.py bigwigFeaturesFile
```
	hardcoded inputs:
		refgene_fname = "refGene.txt"
	    enhancerRangeFile = "sandelin_enh_expanded.bed"

# Model Output format:
	# Index	[ integer ]
	# Gene_name	[ string ]
	# Feature_chr [ chr# or chr ##]
	# Feature_start [ integer ]
	# Feature_end  [ integer ] 
	# Feature_pos_relative_to_gene  [ up,down,body ]
	# Feature_score [ float ]
