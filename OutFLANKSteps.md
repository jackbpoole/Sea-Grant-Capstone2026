---
title: "OutFLANK Analysis"
author: "Larissa Firmansyah"
date: "2026-05-29"
output: html_document
---

This R Markdown file runs OutFLANK on the filtered and clumped Olympia oyster SNP dataset. The final outputs are:

- `outflank_results_(geno_label).csv`
- `outflank_outliers_(geno_label).csv`

Set your parameters at the top.

```{r setup}
library(vcfR)
library(OutFLANK)

# -----------------------------
# Parameters
# -----------------------------

# Choose which filtered/clumped VCF to run
geno_label <- "geno10"
# geno_label <- "geno075"

# Input VCF
vcf_file <- "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/Preprocessing/n218_filtered_snps_geno_10_clumped.recode.vcf"
# vcf_file <- "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/Preprocessing/n218_filtered_snps_geno_075_clumped.recode.vcf"

# Metadata file
metadata_file <- "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/outflank/sample_table_218_clean (2).csv"

# Output directory
out_dir <- "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/outflank"

# OutFLANK settings
hmin_value <- 0.1
q_threshold <- 0.05
```

## Read metadata

```{r metadata}
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

metadata$Estuary <- trimws(metadata$Estuary)

cat("Number of samples in metadata:", nrow(metadata), "\n")
print(table(metadata$Estuary))
```

## Read VCF and extract genotype calls

```{r read-vcf}
vcf <- read.vcfR(vcf_file, verbose = FALSE)

fixdf <- getFIX(vcf)
locusNames <- paste0(fixdf[, "CHROM"], ":", fixdf[, "POS"])

gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

cat("VCF dimensions, loci x individuals:\n")
print(dim(gt))
```

## Convert genotypes to dosage format

OutFLANK needs diploid genotype data coded as 0, 1, or 2, where the value represents the number of alternate alleles. Missing genotypes are converted to `NA` first.

```{r genotype-dosage}
gt_to_dosage <- function(x) {
  if (is.na(x) || x %in% c(".", "./.", ".|.")) return(NA_real_)
  x <- gsub("\\|", "/", x)
  a <- strsplit(x, "/", fixed = TRUE)[[1]]
  if (length(a) != 2 || any(a == ".")) return(NA_real_)
  vals <- suppressWarnings(as.numeric(a))
  if (any(is.na(vals))) return(NA_real_)
  s <- sum(vals)
  if (!s %in% 0:2) return(NA_real_)
  s
}

dos <- apply(gt, c(1, 2), gt_to_dosage)

# OutFLANK expects rows = individuals and columns = SNPs
SNPmat <- t(dos)

cat("SNP matrix dimensions, individuals x loci:\n")
print(dim(SNPmat))
```

## Match samples to estuary populations

```{r match-populations}
sampleNames <- colnames(gt)

# Sample names are formatted like 1_1, 10_10, etc.
sampleNums <- as.integer(sub("_.*", "", sampleNames))

popNames <- metadata$Estuary[
  match(sampleNums, metadata$Sample.number)
]

cat("Number of unmatched samples:", sum(is.na(popNames)), "\n")
print(table(popNames))

head(data.frame(sampleNames, sampleNums, popNames), 20)
```

## Prepare data for OutFLANK

OutFLANK uses `9` as the missing data code.

```{r prepare-outflank}
stopifnot(length(popNames) == nrow(SNPmat))
stopifnot(length(locusNames) == ncol(SNPmat))
stopifnot(sum(is.na(popNames)) == 0)

SNPmat[is.na(SNPmat)] <- 9
storage.mode(SNPmat) <- "numeric"
```

## Calculate FST values

```{r calculate-fst}
FSTdf <- MakeDiploidFSTMat(
  SNPmat = SNPmat,
  locusNames = locusNames,
  popNames = popNames
)

cat("FST table dimensions:\n")
print(dim(FSTdf))

summary(FSTdf$FST)
```

## Run OutFLANK

OutFLANK identifies loci with elevated FST relative to the genome-wide neutral distribution. SNPs with q-values below 0.05 were considered significant outliers.

```{r run-outflank}
out <- OutFLANK(
  FSTdf,
  NumberOfSamples = length(unique(popNames)),
  LeftTrimFraction = 0.05,
  RightTrimFraction = 0.05,
  Hmin = hmin_value,
  qthreshold = q_threshold
)

cat("Outlier summary:\n")
print(table(out$results$OutlierFlag, useNA = "ifany"))

cat("q-value summary:\n")
print(summary(out$results$qvalues))
```

## Save results

```{r save-results}
outliers <- subset(out$results, OutlierFlag == TRUE)
outliers <- outliers[order(outliers$qvalues), ]

all_results_file <- file.path(
  out_dir,
  paste0("outflank_results_n218_", geno_label, "_clumped.csv")
)

outliers_file <- file.path(
  out_dir,
  paste0("outflank_outliers_n218_", geno_label, "_clumped.csv")
)

write.csv(out$results, all_results_file, row.names = FALSE)
write.csv(outliers, outliers_file, row.names = FALSE)

cat("Saved all results to:", all_results_file, "\n")
cat("Saved outliers to:", outliers_file, "\n")
cat("Number of outliers:", nrow(outliers), "\n")

head(outliers)
```
