Prerequites:
R version >3.0 or 
Python version >2.7 

Current call:
python pipeline_top.py bigwigFeaturesFile
	hardcoded inputs:
			refgene_fname = "refGene.txt"
		  enhancerRangeFile = "sandelin_enh_expanded.bed"

# File input list format (tab-separated/tsv):
[patient/cell name] [assay type] [absolute file location]  [aggregation function]



# Model Output format:
	# Index	[ integer ]
	# Gene_name	[ string ]
	# Feature_chr [ chr# or chr ##]
	# Feature_start [ integer ]
	# Feature_end  [ integer ] 
	# Feature_pos_relative_to_gene  [ up,down,body ]
	# Feature_score [ float ]




