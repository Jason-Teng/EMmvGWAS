# Multivariate GWAS in R

This repository contains an R implementation of a multivariate genome-wide association study (GWAS) method designed to identify genetic markers associated with multiple correlated phenotypes simultaneously.

## Features
- Supports multi-trait and single-trait analysis of continuous traits
- Incorporates a kinship matrix to account for genetic relatedness
- Efficient implementation optimized for moderate to high phenotype settings

## Inputs
**Phenotype, P**
- `P` should be a numeric matrix of dimensions `n x m`, where:
  - `n` = number of individuals (rows)
  - `m` = number of traits (columns)
- Column names should correspond to trait names.

  
**Genotype, Z** (SNPs as rows)
- `P` should be a numeric matrix of dimensions `p x n`, where:
  - `p` = number of SNPs (rows)
  - `n` = number of individuals (columns)

**Kinship matrix, K** 
- `K` should be a square matrix of dimensions `n x n`
  - `n` = number of individuals
  - Row/column names should match the individual order in P and Z

## Installation

Clone the repository and make sure the profiling_gwas.R file is in the same folder where you're working.

```bash
git clone https://github.com/Jason-Teng/mvGWAS.git
```

In R:

```r
# Ensure you're in the same folder or provide full path
source("profiling_gwas.R")
```

## Requirements

This method is implemented in base R and depends on the following packages:

```r
install.packages(c("doParallel", "foreach", "RSpectra"))
```

Tested on R version ≥ 4.0.

## Example Usage

```r
library(doParallel)
library(foreach)
library(RSpectra)
result <- gwas_all(K, P, Z)
```

## Output Format

The output from `semi_gwas()` and `gwas_all()` is a data frame with one row per SNP and the following columns:

| Column     | Description                                                  |
|------------|--------------------------------------------------------------|
| `snp`      | SNP index (1-based order)                                    |
| `trait1`, `trait2`, ..., `trait_m` | Estimated intercept/effect for each trait              |
| `wald`     | Wald test statistic for the multivariate test                |
| `p`        | Raw p-value from the Wald test                               |
| `log10p`   | Negative log10 of the p-value (−log₁₀(p))                    |

This object is returned as an R data frame and is also written to a CSV file named <name>_association.csv 

### Example Output for two traits
```r
head(res)
```

| snp |   yd     |    tp     |   wald   |    p     | log10p  |
|-----|----------|-----------|----------|----------|---------|
|  1  | -0.0283  | 0.0532    | 0.1473   | 0.9290   | 0.0320  |
|  2  | 0.0232   | -0.0006   | 0.0190   | 0.9905   | 0.0041  |
|  3  | -0.0281  | -0.0352   | 0.0729   | 0.9642   | 0.0158  |
|  4  | -0.0699  | -0.0352   | 0.2083   | 0.9019   | 0.0452  |
|  5  | -0.1854  | -0.0820   | 1.3718   | 0.5036   | 0.2979  |
|  6  | -0.1668  | -0.1287   | 1.5141   | 0.4690   | 0.3288  |


## Example

Example scripts and data are available in the `example/` directory.

## Citation

If you use this method in your research, please cite:

> Teng, C.-S., Xu, S., et al. (2025). A Multivariate GWAS Method for High-Dimensional Phenotypes. *In preparation.*

## License

This project is licensed under the MIT License. See the `LICENSE` file for full terms.

## Contact

For questions, suggestions, or collaboration inquiries, please contact:

- Chin-Sheng Teng (Jason): chinsheng.teng@email.ucr.edu  
- Shizhong Xu: shxu@ucr.edu

---

Developed at the University of California, Riverside.

