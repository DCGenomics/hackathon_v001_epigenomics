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
2+ sets (patient/replicates) of bigWig files. Inputs are given in an input text of format (including a header):
[patient/cell name] [assay type] [absolute file location]  [aggregation function]

### Example:
Patient	DataType	DataFile	Aggregation
Patient1	RNA	/epigenomes/teamdata/h1_data/GSM484nnn/GSM484408/suppl/GSM484408_UCSF-UBC.H1.mRNA-Seq.H1EScd1_batch2_vial2.wig.gz.bw	mean
Patient2	RNA	/epigenomes/teamdata/h1_data/GSM915nnn/GSM915328/suppl/GSM915328_UCSD.H1.mRNA-Seq.polyA-RNA-seq_h1_r1a.wig.gz.bw	mean
Patient3	RNA	/epigenomes/teamdata/h1_data/GSM484nnn/GSM484408/suppl/GSM484408_UCSF-UBC.H1.mRNA-Seq.H1EScd1_batch2_vial2.wig.gz.bw	mean
Patient1	H3K4me1	/epigenomes/teamdata/h1_data/GSM409nnn/GSM409307/suppl/GSM409307_UCSD.H1.H3K4me1.LL228.wig.gz.bw	sum
Patient2	H3K4me1	/epigenomes/teamdata/h1_data/GSM433nnn/GSM433177/suppl/GSM433177_BI.H1.H3K4me1.Solexa-10529.wig.gz.bw	sum
Patient3	H3K4me1	/epigenomes/teamdata/h1_data/GSM434nnn/GSM434762/suppl/GSM434762_UCSF-UBC.H1.H3K4me1.H1EScd1-me1K4-D.wig.gz.bw	sum
Patient1	H3K4me3	/epigenomes/teamdata/h1_data/GSM409nnn/GSM409308/suppl/GSM409308_UCSD.H1.H3K4me3.LL227.wig.gz.bw	sum
Patient2	H3K4me3	/epigenomes/teamdata/h1_data/GSM410nnn/GSM410808/suppl/GSM410808_UCSF-UBC.H1.H3K4me3.H1EScd1-me3K4-A.wig.gz.bw	sum
Patient3	H3K4me3	/epigenomes/teamdata/h1_data/GSM432nnn/GSM432392/suppl/GSM432392_UCSF-UBC.H1.H3K4me3.H1EScd2-me3K4-C.wig.gz.bw	sum
Patient1	RRBS	/epigenomes/teamdata/h1_data/GSM621nnn/GSM621357/suppl/GSM621357_BI.H1.RRBS.RRBS360-361.wig.gz.bw	mean
Patient2	RRBS	/epigenomes/teamdata/h1_data/GSM621nnn/GSM621705/suppl/GSM621705_BI.H1.RRBS.RRBS776-779.wig.gz.bw	mean
Patient3	RRBS	/epigenomes/teamdata/h1_data/GSM621nnn/GSM621763/suppl/GSM621763_BI.H1.RRBS.RRBS534-535.wig.gz.bw	mean

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
