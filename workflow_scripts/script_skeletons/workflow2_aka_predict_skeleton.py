#! /bin/bash
# Expecting to get called as:
#	predict user_input_list.txt  model_output gene_model.bed.txt




# Read in model_output, create model object(s)

for each patient in user_input:
	#read in user data
	for each gene in gene_list:
		#predict transcription
		line=(patient, gene, prediction)
		output.write(line)





