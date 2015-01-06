library(data.table)
library(plyr)
library(reshape2)
library(glmnet)
library(preprocessCore)

rangeFile = "/epigenomes/teamdata/regions_typed.tab"
rangeEnhancerFile = "/epigenomes/teamdata/sandelin_enh_expanded.bed"
inputListFile = "/home/cemmeydan/inputTest3.txt"
outputFile = "/epigenomes/teamdata/H1_dummyData.txt"
numCores = 10


read.tsv = function(file, sep="\t", header=T)
{
	return(read.table(file, sep=sep, header=header, stringsAsFactors=F,  na.strings = c("NA","#N/A","N/A")))
}

write.tsv = function(object, file, sep="\t", header=T)
{
	write.table(object, file, sep=sep, col.names=header, quote=F, row.names=F)
}


############------ Featurization Start ------############

## Summarizes a list of signal from a list of files
## Input: 
GetRegionSignal = function(ranges, inputList)
{
	rangesUniq = unique(data.frame(ranges[,c("chr","start","end","id")]))
	write.tsv(rangesUniq, rangeFileTmpOut, header=F)

	summaryData = data.frame()
	for(i in 1:nrow(inputList))
	{
		inputRow = inputList[i,]
		inFile = inputRow$DataFile
		tempOut = paste0( basename(inFile), ".summary.bed")
		command = paste0("/home/cemmeydan/bigWigAverageOverBed ", inFile, " ", rangeFileTmpOut, " ", tempOut, " -minMax")
		system(command)
		curSummary = read.table(tempOut, sep="\t", header=F, stringsAsFactors=F)
		colnames(curSummary) = c("id","size","covered","sum","mean0","mean","min","max")
		curSummary$patient = inputRow$Patient
			
		if( ! inputRow$Aggregation %in% c("min","max","mean","mean0","sum")) inputRow$Aggregation = "mean"
									   
		curSummary = curSummary[, c("id","patient",inputRow$Aggregation)]
		colnames(curSummary)[3] = inputRow$DataType
		curSummary = melt(curSummary, c("id","patient"))
		
		summaryData = rbind(summaryData, curSummary)
		file.remove(tempOut, showWarnings=F)
	}
	file.remove(rangeFileTmpOut, showWarnings=F)
	summaryData2 = dcast(summaryData, id+patient ~ variable, fun.aggregate=mean) ## If everything is correct we shouldn't need fun.aggregate now
	return(summaryData2)
}

rangeFileTmpOut = paste0(basename(rangeFile), ".uniq.bed")

inputList = read.tsv(inputListFile)
ranges = read.tsv(rangeFile, header=F)
colnames(ranges) = c("id","gene","region","chr","start","end")
ranges$id = paste(ranges$chr, ranges$start, ranges$end, sep=".")

###
#ranges = ranges[1:600,]

rangesEnhancers = read.tsv(rangeEnhancerFile, header=T, sep=" ")
rangesEnhancers = rangesEnhancers[,1:3]
rownames(rangesEnhancers) = NULL
colnames(rangesEnhancers) = c("chr","start","end")
rangesEnhancers$id = paste(rangesEnhancers$chr, rangesEnhancers$start, rangesEnhancers$end, sep=".")

summaryData = GetRegionSignal(ranges, inputList)
#summaryDataEnh = GetRegionSignal(rangesEnhancers, inputList)

geneData = merge(ranges, summaryData, by="id")
geneData2 = geneData[,c(1:2,6:ncol(geneData))]

write.tsv(geneData2, outputFile)

############------ Modeling Start ------############


modelRNA = function(i, geneDataList, rna)
{
    geneData = geneDataList[[i]]
    geneData = dcast(data = geneData, formula = gene + patient ~ id + variable)
    covari <- as.matrix(geneData[,-c(1:3)])
    # first local cpg model
    step1 <- cv.glmnet(covari, rna, standardize=TRUE)
    s1.c <- predict(step1, type="coefficients", s="lambda.1se")
    return(s1.c)
}



# possible normalization strategy for at least histones
aqn <- function(dF)
{
    cn <- colnames(dF)
    rn <- rownames(dF)
    varstab <- cbind(apply(dF, 2, asinh))
    dF.scale <- normalize.quantile(varstab)
    colnames(dF.scale) <- cn
    rownames(dF.scale) <- rn
    dF.scale
}

better.scale <- function(mat)
{
    nmat <- apply(mat, 2, function(xx){
        if (all(xx == 0)){ return(xx) }
        else{ return(scale(xx)) } })
    return(nmat)
}

geneData2 = read.tsv(outputFile)
geneDataList = split(geneData2, f = geneData2$gene )

#res = t(simplify2array(parallel::mclapply(1:length(geneDataList), modelRNA, geneDataList, mc.cores=numCores)))
res = t(parallel::mclapply(1:length(geneDataList), modelRNA, geneDataList,  mc.cores=4))





