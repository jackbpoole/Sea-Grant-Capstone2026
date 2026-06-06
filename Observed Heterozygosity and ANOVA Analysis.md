# Genetic Diversity ANOVA

## Purpose

This script calculates observed heterozygosity for each Olympia oyster and uses ANOVA to test whether genetic diversity differs among estuary populations.

Observed heterozygosity was used as the genetic diversity metric.

## Required Package

```r
library(vcfR)
```

`vcfR` is used to read the filtered VCF file and extract genotype information.

## Input Files

```r
vcf_file <- "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/Preprocessing/n218_filtered_snps_geno_075_clumped.recode.vcf"

metadata_file <- "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/outflank/sample_table_218_clean (2).csv"
```

The VCF file contains the filtered and LD-clumped SNP dataset.

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

This loads the metadata and checks how many oysters are assigned to each estuary.

## Read VCF File

```r
vcf <- read.vcfR(
  vcf_file,
  verbose = FALSE
)
```

This imports the filtered SNP dataset into R.

## Extract Genotypes

```r
gt <- extract.gt(
  vcf,
  element = "GT",
  as.numeric = FALSE
)

dim(gt)
```

This extracts genotype calls from the VCF.

Rows represent SNPs and columns represent individual oysters.

## Match Samples to Estuaries

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

This matches each oyster in the VCF to its estuary population using the sample number.

`sum(is.na(popNames))` should return 0, meaning all samples were successfully matched.

## Clean Genotype Calls

```r
gt_clean <- gt

gt_clean[] <- gsub(
  "\\|",
  "/",
  gt_clean
)
```

This converts phased genotypes such as `0|1` into unphased genotypes such as `0/1`.

This makes heterozygous genotypes easier to count.

## Identify Heterozygous Genotypes

```r
is_het <- gt_clean %in% c(
  "0/1",
  "1/0"
)

dim(is_het) <- dim(gt_clean)
```

This identifies heterozygous SNPs for each individual.

A heterozygous genotype means the individual has two different alleles at that SNP.

## Identify Successfully Called Genotypes

```r
is_called <- !(
  is.na(gt_clean) |
    gt_clean %in% c(
      ".",
      "./.",
      ".|."
    )
)

dim(is_called) <- dim(gt_clean)
```

This identifies genotypes that were successfully called.

Missing genotypes are excluded from the heterozygosity calculation.

## Calculate Observed Heterozygosity

```r
Ho <- colSums(
  is_het,
  na.rm = TRUE
) /
  colSums(
    is_called,
    na.rm = TRUE
  )
```

Observed heterozygosity was calculated for each individual as:

```text
Ho = number of heterozygous SNPs / number of successfully called SNPs
```

Higher Ho values indicate that a larger proportion of SNPs are heterozygous in that individual.

## Create Diversity Dataframe

```r
diversity_df <- data.frame(
  Sample = sampleNames,
  Sample.number = sampleNums,
  Estuary = popNames,
  Ho = Ho
)

head(diversity_df)
```

This creates a dataframe containing each oyster, its estuary population, and its observed heterozygosity value.

## Calculate Mean Ho by Estuary

```r
aggregate(
  Ho ~ Estuary,
  diversity_df,
  mean
)
```

This calculates the average observed heterozygosity for each estuary population.

## Run One-Way ANOVA

```r
anova_Ho <- aov(
  Ho ~ Estuary,
  data = diversity_df
)

summary(anova_Ho)
```

This tests whether mean observed heterozygosity differs among estuaries.

The response variable is `Ho`.

The grouping variable is `Estuary`.

## Run Tukey Post-Hoc Test

```r
TukeyHSD(anova_Ho)
```

Tukey's post-hoc test compares all pairs of estuaries to determine whether any specific populations differ significantly from each other.

## Save Results

```r
write.csv(
  diversity_df,
  "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/outflank/heterozygosity_by_individual_n218_geno075.csv",
  row.names = FALSE
)
```

This saves the individual heterozygosity values for visualization and downstream analysis.

## Results Summary

For the `geno075` dataset, no significant differences in observed heterozygosity were detected among estuaries.

```text
ANOVA: F6,209 = 1.011, p = 0.419
```

Mean observed heterozygosity values were similar across estuaries, ranging from approximately 0.179 to 0.188.

These results suggest that overall genome-wide genetic diversity is relatively consistent across Olympia oyster populations, even though specific loci showed population differentiation in OutFLANK and pcadapt analyses.
