library(lars)
library(GenomicRanges)
library(preprocessCore)

# competing models :
# simple multi-term linear model
# lasso model

# remember: bw summary * length of region

hist.norm <- function(dF){
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
