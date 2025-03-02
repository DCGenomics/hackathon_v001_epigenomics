args = commandArgs(trailingOnly = TRUE)
if(length(args) < 3)
{
	cat("Usage: \nRscript featurizeModel.r geneRangeFile enhancerRangeFile inputBigWigFile [signalOutput] [coefficientOutput] [numCores] [maxDistanceFromEnhancer]\n\n")
	cat("Defaults:\nsignalOutput = ./outputSignal.txt\ncoefficientOutput = ./outputCoefficient.txt & ./outputCoefficient.txt.filtered \nnumCores = 8 \nmaxDistanceFromEnhancer = 100000\n")
	stop("Error: Need at least 3 inputs\n\n")
}
rangeGeneFile = args[1]
rangeEnhancerFile = args[2]
inputListFile = args[3]
outputSummaryFile = "outputSignal.txt"
outputCoefsFile = "outputCoefficient.txt"
numCores = 8
enhancerProximity = 1e+05

if(length(args) >= 4) outputSummaryFile = args[4]
if(length(args) >= 5) outputCoefsFile = args[5]
if(length(args) >= 6) numCores = as.integer(args[6])
if(length(args) >= 7) enhancerProximity = as.integer(args[7])

suppressMessages(library(data.table)    )
suppressMessages(library(plyr)          )
suppressMessages(library(reshape2)      )
suppressMessages(library(glmnet)        )
suppressMessages(library(preprocessCore))
suppressMessages(library(GenomicRanges) )



#rangeGeneFile = "/epigenomes/teamdata/regions_typed.tab"
#rangeEnhancerFile = "/epigenomes/teamdata/sandelin_enh_expanded.bed"
#inputListFile = "/epigenomes/teamdata/test_data/bigwig_files.txt"
#outputSummaryFile = "allBigWigFiles.txt"
#outputCoefsFile = "allBigWigFiles_coefs.txt"
#numCores = 10
#enhancerProximity = 1e+05





read.tsv = function(file, sep="\t", header=T)
{
	return(read.table(file, sep=sep, header=header, stringsAsFactors=F,  na.strings = c("NA","#N/A","N/A")))
}

write.tsv = function(object, file, sep="\t", header=T)
{
	write.table(object, file, sep=sep, col.names=header, quote=F, row.names=F)
}

############------ Featurization Start ------############

## GetRegionSignal: Summarizes a list of signal for a list of given intervals
## Input: 
## ranges = data.frame, giving the regions the signals will be summarized over (gene bodies, introns, enhancers...). MUST contain the fields chr/start/end and a unique field named id
## inputList = data.frame, list of *bigWig* files that contains the signal (histone peaks, TFBS, DNA methylation, RNAseq). 
##            	Row format is of (Patient, DataType, DataFile, Aggregation) where each field is as described below. Each unique patient MUST contain a row with DataType of "RNA".
##					Patient: A patient id
##					DataType: A text describing the signal. 
##					DataFile: Path to the bigWig file
##					Aggregation: how to summarize the data over the peak. can take values (mean, min, max, sum, mean0). mean will give the mean value in only the covered bases whereas mean0 will average over zeroes as well
## Output: data.frame containing tuples of (gene, region, patient, DataType, Summarized value)
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



print("Reading inputs...")

rangeFileTmpOut = paste0(basename(rangeGeneFile), ".uniq.bed")

## process gene ranges
inputList = read.tsv(inputListFile, header=T)
colnames(inputList) = c("Patient", "DataType", "DataFile", "Aggregation")
rangesGene = read.tsv(rangeGeneFile, header=F)
colnames(rangesGene) = c("id","gene","region","chr","start","end")
rangesGene$id = paste(rangesGene$chr, rangesGene$start, rangesGene$end, sep=".")

## process enhancer ranges
rangesEnhancers = read.tsv(rangeEnhancerFile, header=T, sep=" ")
rangesEnhancers = rangesEnhancers[,1:3]
rownames(rangesEnhancers) = NULL
colnames(rangesEnhancers) = c("chr","start","end")
rangesEnhancers$id = paste(rangesEnhancers$chr, rangesEnhancers$start, rangesEnhancers$end, sep=".")
rangesEnhancers$region = paste0("enhancer.", rangesEnhancers$id)

## find enhancers within {enhancerProximity} distance to each gene
rangesGenebody = unique(rangesGene[rangesGene$region == "body",c("gene","chr","start","end")])
grGeneFlanking = GRanges(seqnames=rangesGenebody$chr, ranges = IRanges(start=rangesGenebody$start-enhancerProximity, end=rangesGenebody$end+enhancerProximity), gene=rangesGenebody$gene)
grEnhancer = GRanges(seqnames=rangesEnhancers$chr, ranges = IRanges(start=rangesEnhancers$start, end=rangesEnhancers$end), id=rangesEnhancers$id, region=rangesEnhancers$region)
overlapGeneEnh = findOverlaps(grEnhancer, grGeneFlanking)
overlapGeneEnh2 = data.frame(id=grEnhancer$id[overlapGeneEnh@queryHits], gene=grGeneFlanking$gene[ overlapGeneEnh@subjectHits], region=grEnhancer$region[overlapGeneEnh@queryHits], chr=grEnhancer@seqnames[overlapGeneEnh@queryHits], start=grEnhancer@ranges@start[overlapGeneEnh@queryHits], end=grEnhancer@ranges@start[overlapGeneEnh@queryHits]+grEnhancer@ranges@width[overlapGeneEnh@queryHits]-1)

rangesGenesEnh = rbind(rangesGene, overlapGeneEnh2)
write.tsv(rangesGenesEnh, "ranges_GenesAndProximalEnhancers.txt")

## Calculate the signal for the genic and enhancer regions

print("Calculating signal for the inputs...")

summaryData = GetRegionSignal(rangesGenesEnh, inputList)
geneData = merge(rangesGenesEnh, summaryData, by="id")
geneData2 = geneData[,c(1:3,7:ncol(geneData))]
geneData2 = melt(geneData2, c("id", "gene", "region", "patient"))
geneData2 = geneData2[ ! (geneData2$variable == "RNA" & geneData2$region != "body"), ]

write.tsv(geneData2, outputSummaryFile)


############------ Modeling Start ------############
modelRNA.bak <- function(i, geneNames, geneDataAll)
{
	geneName = geneNames[i]
	geneData = geneDataAll[geneDataAll$gene == geneName, ]
    #geneData = geneDataList[[i]]


    rnaData = geneData[ geneData$variable=="RNA", c("gene","patient","value")]
    colnames(rnaData)[3] = "RNA"
    geneData = geneData[geneData$variable != "RNA", ]
    geneData = dcast(geneData, gene+patient~region+variable, fun.aggregate=mean)
    geneData = merge(rnaData, geneData, by=c("gene","patient"))	
    metaData.id <- grep("gene|patient", colnames(geneData))
    rna.id <- which(colnames(geneData) == "RNA")
    enh.id <- grep("enhancer", colnames(geneData))
    covari <- better.scale(as.matrix(geneData[,-c(metaData.id, rna.id, enh.id)]))
    enh <- better.scale(as.matrix(geneData[,enh.id]))
    rna <- as.matrix(geneData[,rna.id])
    if (all(rna == 0))
	{
		if (length(enh.id) != 0){ covari <- cbind(covari, enh) }
		return(data.frame(gene=unique(geneData$gene)[1],
						   variable=c("(Intercept)", colnames(covari)),
						   coefficient=rep.int(0, 1+ncol(covari)), row.names=NULL))
	}
    rna <- log2((rna+0.5)/sum(rna+1)*1e6)
    # first local cpg model
    tryCatch({
        step1 <- cv.glmnet(covari, rna, standardize=TRUE)
        s1.c <- predict(step1, type="coefficients", s="lambda.1se")
        rnames.s.c <- rownames(s1.c)
        c.s.c <- s1.c[,1]
        s1.fit <- predict(step1, newx=covari, s="lambda.1se")
        residual <- rna - s1.fit
        if (length(enh.id) != 0){
            step2 <- cv.glmnet(enh, residual, standardize=TRUE, intercept=FALSE)
            s2.c <- predict(step2, type="coefficients", s="lambda.1se")
            rnames.s.c <- c(rnames.s.c, rownames(s2.c)[-1])
            c.s.c <- c(c.s.c, s2.c[-1,1])
        }
        coefs = data.frame(gene=unique(geneData$gene), variable=rnames.s.c, coefficient=c.s.c, row.names=NULL)
        return(coefs)
    }, error=function(e){
                if (length(enh.id) != 0){ covari <- cbind(covari, enh) }
	    	    return(data.frame(gene=unique(geneData$gene),
                                      variable=c("(Intercept)", colnames(covari)),
                                      coefficient=rep.int(0, ncol(covari)+1), row.names=NULL))
    })## end tryCatch
}

modelRNA <- function(i, geneNames, geneDataAll)
{
    geneName = geneNames[i]
    geneData = geneDataAll[geneDataAll$gene == geneName, ]
    #geneData = geneDataList[[i]]

    rnaData = geneData[ geneData$variable=="RNA", c("gene","patient","value")]
    colnames(rnaData)[3] = "RNA"
    geneData = geneData[geneData$variable != "RNA", ]
    geneData = dcast(geneData, gene+patient~region+variable, fun.aggregate=mean)
    geneData = merge(rnaData, geneData, by=c("gene","patient"))	
    metaData.id <- grep("gene|patient", colnames(geneData))
    rna.id <- which(colnames(geneData) == "RNA")
    covari <- better.scale(as.matrix(geneData[,-c(metaData.id, rna.id)]))
    rna <- as.matrix(geneData[,rna.id])
    if (all(rna == 0))
	{
		return(data.frame(gene=unique(geneData$gene)[1],
						   variable=c("(Intercept)", colnames(covari)),
						   coefficient=rep.int(0, 1+ncol(covari)), row.names=NULL))
	}
    rna <- log2((rna+0.5)/sum(rna+1)*1e6)
    # first local cpg model
    tryCatch({
        step1 <- cv.glmnet(covari, rna, standardize=TRUE)
        s1.c <- predict(step1, type="coefficients", s="lambda.1se")
        rnames.s.c <- rownames(s1.c)
        c.s.c <- s1.c[,1]
        coefs = data.frame(gene=unique(geneData$gene), variable=rnames.s.c, coefficient=c.s.c, row.names=NULL)
        return(coefs)
    }, error=function(e){
	    	    return(data.frame(gene=unique(geneData$gene),
                                      variable=c("(Intercept)", colnames(covari)),
                                      coefficient=rep.int(0, ncol(covari)+1), row.names=NULL))
    })## end tryCatch
}


# modelRNA.fast <- function(i, geneDataList){
#     geneData = geneDataList[[i]]
#     rnaData = geneData[ geneData$variable=="RNA", c("gene","patient","value")]
#     colnames(rnaData)[3] = "RNA"
#     geneData = geneData[geneData$variable != "RNA", ]
#     geneData = dcast(geneData, gene+patient~region+variable, fun.aggregate=mean)
#     geneData = merge(rnaData, geneData, by=c("gene","patient"))
# 	apply(geneData, 2, sd)
# 	
#     metaData.id <- grep("gene|patient", colnames(geneData))
#     rna.id <- which(colnames(geneData) == "RNA")
#     covari <- as.matrix(geneData[,-c(metaData.id, rna.id)])
#     rna <- as.matrix(geneData[,rna.id])
#     if (all(rna == 0))
# 	{
# 		return(data.frame(gene=unique(geneData$gene), variable=c("(Intercept)", colnames(covari)), coefficient=rep.int(0, ncol(covari)+1), row.names=NULL))
# 	}
#     rna <- log2((rna+0.5)/sum(rna+1)*1e6)
#     # first local cpg model
#     step1 <- cv.glmnet(covari, rna, standardize=TRUE)
#     s1.c <- predict(step1, type="coefficients", s="lambda.1se")
#     coefs = data.frame(gene=unique(geneData$gene), variable=rownames(s1.c), coefficient=s1.c[,1], row.names=NULL)
# 	return(coefs)
# }


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
        if (length(unique(xx) == 1)){ return(xx) }
        else{ return(scale(xx)) } })
    return(nmat)
}

#geneData2 = read.tsv(outputSummaryFile)
geneNames = unique(geneData2$gene)
#geneDataList = split(geneData2, f = geneData2$gene )

print("Building the model...")

results = parallel::mclapply(1:length(geneNames), modelRNA, geneNames, geneData2, mc.cores=numCores)

results2 = rbindlist(results)
write.tsv(results2, outputCoefsFile)
results3 = results2[results2$variable != "(Intercept)" & results2$coefficient != 0, ]
write.tsv(results3, paste0(outputCoefsFile, ".filtered"))

print("Done!")

#res = rbindlist(parallel::mclapply(1:100, modelRNA, geneDataList, mc.cores=numCores))

results = parallel::mclapply(1:length(geneNames), function(xx){
    tryCatch({
      modelRNA(xx, geneNames, geneData2)
    }, error=function(e){
        print(e)
        print(xx)
    })
}, mc.cores=numCores)

