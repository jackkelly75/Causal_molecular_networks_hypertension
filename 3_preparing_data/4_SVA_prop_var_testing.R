#module load apps/gcc/R/4.1.0
#module load tools/env/proxy
library(sva)
library(limma)
#import the data
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
##import the adjusted data
load("data_adj_limma.Rdata") #data_adj
load("covariate.Rdata") #mod1, covariate_data

###Identify Proportion of Variance for each SV
#svd() function
#dat = data
dat = data_adj
mod = mod1
mod0 = NULL
n.sv= 100
B = 5

#functions required
edge.lfdr <- function(p, trunc=TRUE, monotone=TRUE, transf=c("probit", "logit"), adj=1.5, eps=10^-8, lambda=0.8, ...) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0 <- min(pi0, 1)
    
    n <- length(p)
    transf <- match.arg(transf)
    
    if(transf=="probit") {
        p <- pmax(p, eps)
        p <- pmin(p, 1-eps)
        x <- qnorm(p)
        myd <- density(x, adjust=adj)
        mys <- smooth.spline(x=myd$x, y=myd$y)
        y <- predict(mys, x)$y
        lfdr <- pi0*dnorm(x)/y
    }
    
    if(transf=="logit") {
        x <- log((p+eps)/(1-p+eps))
        myd <- density(x, adjust=adj)
        mys <- smooth.spline(x=myd$x, y=myd$y)
        y <- predict(mys, x)$y
        dx <- exp(x) / (1+exp(x))^2
        lfdr <- pi0 * dx/y
    }
    
    if(trunc) {
        lfdr[lfdr > 1] <- 1
    }
    if(monotone) {  
        lfdr <- lfdr[order(p)]
        lfdr <- mono(lfdr)
        lfdr <- lfdr[rank(p)]
    }
    
    return(lfdr)
}
mono <- function(lfdr){
    .Call("monotone", as.numeric(lfdr), PACKAGE="sva")
}

#SVA() function
n <- ncol(dat)
m <- nrow(dat)
mod0 <- mod[, 1]
Id <- diag(n)
resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
uu <- eigen(t(resid) %*% resid)
vv <- uu$vectors
ndf <- n - dim(mod)[2]
pprob <- rep(1, m)
one <- rep(1, n)
Id <- diag(n)
df1 <- dim(mod)[2] + n.sv
df0 <- dim(mod0)[2] + n.sv
rm(resid)
cat(paste("Iteration (out of", B, "):"))
#irwsva function within SVA
#This function is the implementation of the iteratively re-weighted least squares approach for estimating surrogate variables
#As a buy product, this function produces estimates of the probability of being an empirical control.
for (i in 1:B) {
    mod.b <- cbind(mod, uu$vectors[, 1:n.sv])
    mod0.b <- cbind(mod0, uu$vectors[, 1:n.sv])
    ptmp <- f.pvalue(dat, mod.b, mod0.b) #parametric F-test P-values adjusted for surrogate variables
    pprob.b <- (1 - edge.lfdr(ptmp))
    mod.gam <- cbind(mod0, uu$vectors[, 1:n.sv])
    mod0.gam <- cbind(mod0)
    ptmp <- f.pvalue(dat, mod.gam, mod0.gam)
    pprob.gam <- (1 - edge.lfdr(ptmp))
    pprob <- pprob.gam * (1 - pprob.b)
    dats <- dat * pprob
    dats <- dats - rowMeans(dats)
    uu <- eigen(t(dats) %*% dats)
    cat(paste(i, " "))
}

#Singular Value Decomposition of adjusted data
temp <- svd(dats)
sv = temp$v[, 1:n.sv, drop = FALSE] #extract the surrogate variables for number wanted
var <- temp$d #get singular values of dats
#create the standard output of sva() function
retval <- list(sv = sv, pprob.gam = pprob.gam, pprob.b = pprob.b, n.sv = n.sv)

#prop of variance
var_explained = var^2 / sum(var^2)
#cumulative proportion
cumulative <- cumsum(var_explained)


#plotting the results
library(ggplot2)
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")

#create scree plot
var_explained <- var_explained[1:100]
png(paste0("scree_plot_SVs.png"), width = 10, height =10 , units = "in", res = 800) 
qplot(c(1:100), var_explained) + 
  geom_line() + 
  xlab("SVs") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 0.4) #limit to 0.4 as variance explained does not go above this for this example
dev.off()

#limit to top 20 for plotting purposes 
var_explained <- var_explained[1:20]
#create scree plot
png(paste0("scree_plot_SVs_cropped.png"), width = 10, height =10 , units = "in", res = 800)
qplot(c(1:20), var_explained) + 
  geom_line() + 
  xlab("SVs") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 0.4) #limit to 0.4 as variance explained does not go above this for this example
dev.off()





####rough way of finding elbow of the scree plots
#adapted from quick.elbow() from bigpca (old r package)
varpc <- var_explained
low=.05 #what is considered low proportion of variance in this data
max.pc=.8 #if sum of prop of variance hits this then consider it the elbow

low_prop <- (varpc<low) #TRUE/FALSE for if prop of variance is lower than 0.05
high_prop <- varpc>low/8  #those that are false are very low prop of var that will defintely not be included
low.ones <- which(low_prop & high_prop) #gets which are lower than low but higher than low/8
keeping <- length(which(!low_prop)) #keeps those that are above low cutoff for prop of variance

set <- varpc[low.ones] #set of the SVs that are lower than low but higher than low/8
pc.drops <- abs(diff(set))/(set[1:(length(set)-1)]) #difference between each of low prop/set without last
max_drop <- which(pc.drops==max(pc.drops,na.rm=T)) #max difference in prop of var (finding biggest drop between SVs in this subset)
elbow <- max_drop+keeping #gives the number we are keeping + the number before max drop in those that are low

#test if cumulative prop is over 0.8 for any
if(tail(cumsum(varpc[1:elbow]),1)>max.pc) {
    elbow <- which(cumsum(varpc)>max.pc)[1]-1
}else{
    print("cumulative not used")
}

elbow
#6




eigenvalues <- var^2
eigenvalues <- eigenvalues[1:300]
qplot(c(1:300), eigenvalues) + 
  geom_line() + 
  xlab("SVs") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1000)
