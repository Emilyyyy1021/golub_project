
## Compare results from different methods of controlling for overall false discovery rate

#  This R script is used for a simulation to estimate the correct
##   p-value and t-statistic to use when adjusting the significance
##   level of a test to control experiment-wise error when conducting
##   many t-tests in a single study.

#  The overall design is as follows:
##   1. Select 100 genes for analysis (you might want to modify the script to keep all genes ## after some filtering for genes with low sd across samples)
##   2. Use the 100 genes to estimate the correlation among the genes
##      by calculating a covariance matrix.  The correlation among
##      the genes will induce a correlation among simulated t-tests
##   3. Draw repeated samples from simulated sets of 100 genes assumed
##      to have mean zero expression and correlation/covariance matching
##      the earlier estimates.
##   4. Calculate and store p-values for the repeated samples.
##   5. Search the empirical distribution of p-values to find the mininum p-value
##      that it would lead to rejection of the hypopthesis of no difference among
##      any of the genes about 5% of the time, i.e., has the correct overall
##      experiment-wise error.

#  Since the genes in the simulation are assumed to have mean expression zero, the 
##    simulated expression levels come from a setting where no genes are associated
##    with leukemia type.

#  First, clear the working directory

rm(list = ls())  # empty the workspace to start fresh


# MASS package needed to simulate a vector of Normal RVs;
##    Install  in R Studio using the Install button in
##    the packages file pane

require("MASS")

load("golub_exprs_pheno.Rdata")
Golub <- golub.exprs.pheno

#create gene.matrix, trimmed version of Golub dataset
gene.matrix = as.matrix(Golub[,-(1:6)])

#create logical variable for cancer type
leuk.type = (Golub$cancer == "aml")

#create a vector of integers from 1 to the total number of genes
gene.columns = 1:ncol(gene.matrix)

#set the seed for a pseudo-random sample
set.seed(2024)

#sample 100 numbers from gene.columns, without replacement
gene.index.set = sample(gene.columns, size = 100, replace = FALSE)

#contains expression info from the rows corresponding to gene.index.set
gene.matrix.sample = gene.matrix[,gene.index.set]

#  In the code below set

##    gene.set = 1:2 for 2 tests using 2 genes
##    gene.set = 1:10 for 10 tests using 10 genes
##    gene.set = 1:25 for 25 tests using 25 genes
##    gene.set = 1:100 for 100 tests using 100 genes


# gene.set = 1:2
# gene.set = 1:10
gene.set = 1:25
# gene.set = 1:100

gene.expression.dim = length(gene.set)
sample.size = nrow(Golub)
num.replicates = 1000

# select appropriate number of genes and set theoretical means to 0

golub.expression.sample = gene.matrix.sample[,1:gene.expression.dim]
gene.expression.mean = rep(0,len = gene.expression.dim)

# use the sample selected to estimate the covariance (similar to correlation) 
# among these genes.  The correlation is important because it effects the extent
# to which t-tests for two genes tend to reject the null hypothesis together

cov.matrix = cov(golub.expression.sample)

## initialize p-value vectors

min.p.value = rep(1, len = num.replicates)
gene.p.value = rep(1, len= gene.expression.dim)

for (ii in 1:num.replicates){
  gene.expression.mat = mvrnorm(sample.size, mu = gene.expression.mean, Sigma = cov.matrix)
  for (jj in 1:gene.expression.dim){
    gene.p.value[jj] = t.test(gene.expression.mat[,jj] ~ leuk.type)$p.val
  }
  min.p.value[ii] = min(gene.p.value)
}

summary(min.p.value)

## explore effect of lowering individual p-value

a = (min.p.value < 0.050)
sum(a)/num.replicates

b = (min.p.value < 0.005)
sum(b)/num.replicates

c = (min.p.value < 0.001)
sum(c)/num.replicates

##  search for correct individual p-value to
##     yield approximately 0.05 experiment-wise
##     error

cut.off.p = quantile(min.p.value, 0.05)
cut.off.p


#### Compare the list of genes you find using the simulation approach above to q value approach for controlling the FDR. 

