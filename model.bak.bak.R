library(data.table)
library(plyr)
library(reshape2)
library(glmnet)
library(preprocessCore)

read.tsv = function(file, header=T)
{
	return(read.table(file, sep="\t", header=header, stringsAsFactors=F,  na.strings = c("NA","#N/A","N/A")))
}

write.tsv = function(object, file, header=T)
{
	write.table(object, file, sep="\t", col.names=header, quote=F, row.names=F)
}


geneData2 = read.tsv("/home/cemmeydan/H1_dummyData.txt")

geneDataList = split(geneData2, f = geneData2$gene )

res = t(parallel::mclapply(1:100, modelRNA, geneDataList,  mc.cores=4))

modelRNA = function(i, geneDataList, rna){
    geneData = geneDataList[[i]]
    geneData = dcast(data = geneData, formula = gene + patient ~ id + variable)
    covari <- as.matrix(geneData[,-c(1:3)])
    # first local cpg model
    step1 <- cv.glmnet(covari, rna, standardize=TRUE)
    s1.c <- predict(step1, type="coefficients", s="lambda.1se")
    return(s1.c)
}



# possible normalization strategy for at least histones
aqn <- function(dF){
    cn <- colnames(dF)
    rn <- rownames(dF)
    varstab <- cbind(apply(dF, 2, asinh))
    dF.scale <- normalize.quantile(varstab)
    colnames(dF.scale) <- cn
    rownames(dF.scale) <- rn
    dF.scale
}

better.scale <- function(mat){
    nmat <- apply(mat, 2, function(xx){
        if (all(xx == 0)){ return(xx) }
        else{ return(scale(xx)) } })
    return(nmat)
}
