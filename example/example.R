# ===================================================
# Multivariate GWAS Example Script using profiling_gwas.R
# ===================================================

# Load required libraries
library(doParallel) 
library(foreach) 
library(RSpectra) 

# Load the multivariate GWAS functions
source("profiling_gwas.R")

# ---------------------------------------------------
# Step 1: Read and prepare phenotype data
# ---------------------------------------------------
# traits are in columns 2 to 5 of the CSV file
Phenotypes = read.csv("RIL-Phenotypes.csv")
P = scale(Phenotypes[,2:5,drop = FALSE]) # standardize the traits
head(P)

# ---------------------------------------------------
# Step 2: Read and prepare genotype data
# ---------------------------------------------------
# Assumes the first column is SNP ID; the rest are genotype values (-1/0/1)
Genotypes <- read.csv("RIL-Genotypes.csv")
Z <- as.matrix(Genotypes[,2:ncol(Genotypes)]) # convert to numeric matrix
head(Z)

# ---------------------------------------------------
# Step 3: Construct the kinship matrix (K)
# ---------------------------------------------------
KK = t(Z)%*%Z
K = KK/mean(diag(KK)) 

# ---------------------------------------------------
# Option A: Run the multivariate GWAS all in one step
# ---------------------------------------------------
# This convenience function wraps everything: eigen decomposition, variance estimation, genome scanning
result <- gwas_all(K, P, Z, name="gwas",num_cores=1)

# ---------------------------------------------------
# Option B: Run in detailed steps (more control)
# ---------------------------------------------------
# Step 1: Eigen decomposition of the kinship matrix
r <- eigen_trans(K)

# Step 2: Estimate variance component lambda for the multivariate model
lambda <- vc_estimation(P, r, name="gwas") 

# Step 3: Perform the semi-exact multivariate GWAS test
result <- semi_gwas(P, r, Z, lambda, name="gwas",num_cores=3)


# ---------------------------------------------------
# Optional: Manhattan Plot Visualization
# ---------------------------------------------------

# Load SNP position information
map <- read.csv("map.csv")

# Build the GWAS result data frame in qqman format
library(qqman)
k= nrow(result) # number of SNPs

gwasinfo <- data.frame(    
  SNP = map$Bin,            # SNP names 
  CHR =  map$Chr,           # chromosome number
  BP = map$Stop,            # SNP base-pair position
  P =result$log10p          # -log10-transformed p-values
)

# Define Bonferroni threshold (suggestive line)
suggestiveline <- -log10(0.05/k)

# Plot Manhattan plot using qqman package
manhattan(gwasinfo,
          main="Multivariate GWAS",
          logp = FALSE,
          genomewideline=F,
          suggestiveline=suggestiveline) # ,genomewideline  = F,suggestiveline = 4.51

