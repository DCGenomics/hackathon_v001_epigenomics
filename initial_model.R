library(lars)
library(GenomicRanges)

# competing models :
# simple multi-term linear model
# lasso model

# remember: bw summary * length of region

getSum <- function(gr, mat){
    return(mat * width(gr))
}

lars.resid <- function(lar, y, x){
    fi <- predict.lars(lar, newx=x, type="fit")
    return(y - fi$fit[,which.min(summary(lar)$Cp)])
}


# general lm - lm(y ~ x, data=as.data.frame(cbind(y, x))) with dim(x) = 442, 10
base.mod <- lm(y ~ cpg + hist1, data=dat)
resid <- residuals(base.mod)

# other mod
# first local cpg model
step1 <- lm(y ~ cpg)
step1.resid <- residuals(step1)
l.lars <- lars(step1.resid, hist1,type="lasso")
l.resid <- lars.resid(l.lars, y, step1.resid)

# some boilerplate for lasso with best model
# h/t http://stats.stackexchange.com/a/121535/10530
cv <- lars::cv.lars(x, y, plot.it=FALSE,mode="step")
idx <- which.max(cv$cv - cv$cv.error <= min(cv$cv))
coef(lars::lars(x,y))[idx,]

# LRT - lasso(Y - Est.Y)/lm(Y - Est.Y)
sigma.base <- (t(resid) %*% resid)/(nrow(resid) - 




lars.ss/resid


## feature extraction
base.sigs <- summary(base.mod)$coefficients[,4]
adv.sigs <- 
