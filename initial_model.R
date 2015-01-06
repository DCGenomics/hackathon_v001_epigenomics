library(lars)
library(GenomicRanges)
library(preprocessCore)
library(bsseq)

# competing models :
# simple multi-term linear model
# lasso model

# remember: bw summary * length of region

aqn <- function(dF){
    cn <- colnames(dF)
    rn <- rownames(dF)
    varstab <- cbind(apply(dF, 2, asinh))
    dF.scale <- normalize.quantile(varstab)
    colnames(dF.scale) <- cn
    rownames(dF.scale) <- rn
    dF.scale
}

getSum <- function(gr, mat){
    return(mat * width(gr))
}

lars.resid <- function(lar, y, x){
    fi <- predict.lars(lar, newx=x, type="fit")
    ## Mallow's Cp for selecting on model size not RSS
    return(y - fi$fit[,which.min(summary(lar)$Cp)])
}

term.split <- function(term){
    return(strsplit(term, ":|-"))
}

## from sketch in cem's notebook i took a picture of
## format: geneID region:cpg [region:histX]*N
proc.line <- function(l){
    l <- strsplit(l,split=" +")[[1]]
    gene <- l[1]
    cpg <- l[2]
    hists <- l[3:length(l)]
    cpg <- term.split(cpg)
    hists <- lapply(hists, term.split)
    return(list(gene=gene, cpg=cpg, hists=hists))
}

## input is a single person's data.frame with lines in format above
## regions are presumably constant across individuals due to
## our region building strategy
proc.dF <- function(dF){
    lins <- apply(dF, 1, proc.line)
    ## make regions
    cpg.regions <- do.call(rbind, lapply(lins, function(xx){
        return(xx$cpg[[1]]) 
    }))
    cpg.vals <- cpg.regions[,4]
    cpg.regions <- GRanges(seqnames=Rle(cpg.regions[,1]),
                           start=as.numeric(cpg.regions[,2]),
                           end=as.numeric(cpg.regions[,3]))
    histones <- lapply(lins, function(xx){
        return(xx$hists)
    })
    


## data object generations


# first local cpg model
step1 <- lm(y ~ cpg)
step1.resid <- residuals(step1)
# then lasso
step2 <- lars(step1.resid, hist1,type="lasso")
# l.resid <- lars.resid(step2, y, step1.resid)
# some debate about using mallow's cp vs cv - cp performs exceptionally well
# according to efron if model is true...
# some boilerplate for lasso with best model
# h/t http://stats.stackexchange.com/a/121535/10530

idx <- which.min(summary(step2)$Cp)
## vs
cv <- lars::cv.lars(x, y, plot.it=FALSE,mode="lasso")
idx <- which.max(cv$cv - cv$cv.error <= min(cv$cv))

## feature extraction
pthresh <- 0.05 # should be modified for multiple testing
which(coef(step2)[idx,] > 0)
which(summary(step1)$coefficients[,4] < p-thresh)
# GOF
summary(step1)$r.squared
step2$R2[idx]
