# OutFLANK Analysis

## Purpose

This script runs OutFLANK to identify SNPs with unusually high genetic differentiation among Olympia oyster estuary populations.

## Required Packages

```r
library(vcfR)
library(OutFLANK)
```

* `vcfR` is used to read VCF files and extract genotype information.
* `OutFLANK` is used to calculate FST values and identify outlier loci.

## Input Files

```r
vcf_file <- "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/Preprocessing/n218_filtered_snps_geno_075_clumped.recode.vcf"

metadata_file <- "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/outflank/sample_table_218_clean (2).csv"
```

The metadata file contains:

* `Sample.number`
* `Estuary`

## Read Metadata

```r
metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

metadata$Estuary <- trimws(metadata$Estuary)

head(metadata)
table(metadata$Estuary)
```

This imports the sample metadata and verifies the number of oysters assigned to each estuary.

## Read VCF File

```r
vcf <- read.vcfR(vcf_file, verbose = FALSE)

fixdf <- getFIX(vcf)

locusNames <- paste0(
  fixdf[, "CHROM"],
  ":",
  fixdf[, "POS"]
)
```

This loads the VCF and creates a unique identifier for each SNP using scaffold and genomic position.

## Extract Genotypes

```r
gt <- extract.gt(
  vcf,
  element = "GT",
  as.numeric = FALSE
)

dim(gt)
```

This extracts genotype calls for all individuals and loci.

## Match Individuals to Estuaries

```r
sampleNames <- colnames(gt)

sampleNums <- as.integer(
  sub("_.*", "", sampleNames)
)

popNames <- metadata$Estuary[
  match(
    sampleNums,
    metadata$Sample.number
  )
]

sum(is.na(popNames))
table(popNames)
```

This links each oyster in the VCF to its estuary of origin.

`sum(is.na(popNames))` should return 0.

## Convert Genotypes to Dosage Format

```r
gt_to_dosage <- function(x) {

  if (is.na(x) || x %in% c(".", "./.", ".|."))
    return(NA_real_)

  x <- gsub("\\|", "/", x)

  a <- strsplit(x, "/", fixed = TRUE)[[1]]

  if (length(a) != 2 || any(a == "."))
    return(NA_real_)

  vals <- suppressWarnings(as.numeric(a))

  if (any(is.na(vals)))
    return(NA_real_)

  s <- sum(vals)

  if (!s %in% 0:2)
    return(NA_real_)

  return(s)
}
```

This function converts genotypes into:

* 0 = homozygous reference
* 1 = heterozygous
* 2 = homozygous alternate

## Build SNP Matrix

```r
dos <- apply(gt, c(1, 2), gt_to_dosage)

SNPmat <- t(dos)

SNPmat[is.na(SNPmat)] <- 9

storage.mode(SNPmat) <- "numeric"
```

OutFLANK requires individuals as rows and SNPs as columns.

Missing values are encoded as 9.

## Calculate FST Values

```r
FSTdf <- MakeDiploidFSTMat(
  SNPmat = SNPmat,
  locusNames = locusNames,
  popNames = popNames
)
```

This calculates locus-specific FST values across estuary populations.

## Run OutFLANK

```r
out <- OutFLANK(
  FSTdf,
  NumberOfSamples = length(unique(popNames)),
  LeftTrimFraction = 0.05,
  RightTrimFraction = 0.05,
  Hmin = 0.1,
  qthreshold = 0.05
)
```

OutFLANK estimates a neutral FST distribution and identifies loci with significantly elevated differentiation.

## Extract Outliers

```r
outliers <- subset(
  out$results,
  OutlierFlag == TRUE
)

outliers <- outliers[
  order(outliers$qvalues),
]
```

Only loci with q-values below 0.05 are retained.

## Save Results

```r
write.csv(
  out$results,
  "outflank_results_n218_geno075_clumped.csv",
  row.names = FALSE
)

write.csv(
  outliers,
  "outflank_outliers_n218_geno075_clumped.csv",
  row.names = FALSE
)
```

The full results file contains all SNPs tested by OutFLANK.

The outlier file contains only significant loci identified as candidate outliers.

