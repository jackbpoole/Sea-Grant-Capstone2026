# -----------------------------
# 1. Load required packages
# -----------------------------

# vcfR is used to read and extract genotype information from VCF files
library(vcfR)

# OutFLANK is used to calculate FST and identify outlier SNPs
library(OutFLANK)


# -----------------------------
# 2. Define input file paths
# -----------------------------

# Path to the filtered and LD-clumped VCF file
vcf_file <- "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/Preprocessing/n218_filtered_snps_geno_075_clumped.recode.vcf"

# Path to metadata file containing sample numbers and estuary names
metadata_file <- "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/outflank/sample_table_218_clean (2).csv"


# -----------------------------
# 3. Read metadata
# -----------------------------

# Read metadata into R
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Remove extra spaces from estuary names
metadata$Estuary <- trimws(metadata$Estuary)

# Check metadata formatting
head(metadata)

# Check sample counts by estuary
table(metadata$Estuary)


# -----------------------------
# 4. Read VCF file
# -----------------------------

# Load the VCF file
vcf <- read.vcfR(vcf_file, verbose = FALSE)

# Extract fixed VCF information, including chromosome/scaffold and position
fixdf <- getFIX(vcf)

# Create unique locus names using scaffold and position
locusNames <- paste0(fixdf[, "CHROM"], ":", fixdf[, "POS"])


# -----------------------------
# 5. Extract genotype calls
# -----------------------------

# Extract genotype field from the VCF
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# Check SNP and sample dimensions
dim(gt)


# -----------------------------
# 6. Match samples to estuary populations
# -----------------------------

# Get sample names from the VCF
sampleNames <- colnames(gt)

# Extract sample number from names like "1_1", "10_10", etc.
sampleNums <- as.integer(sub("_.*", "", sampleNames))

# Match sample numbers to estuary names in metadata
popNames <- metadata$Estuary[
  match(sampleNums, metadata$Sample.number)
]

# Confirm every sample matched to an estuary
sum(is.na(popNames))

# View population counts used in OutFLANK
table(popNames)


# -----------------------------
# 7. Convert genotypes to dosage format
# -----------------------------

# Function to convert VCF genotype calls into 0/1/2 dosage format
# 0 = homozygous reference
# 1 = heterozygous
# 2 = homozygous alternate
# NA = missing or unsupported genotype
gt_to_dosage <- function(x) {
  if (is.na(x) || x %in% c(".", "./.", ".|.")) return(NA_real_)
  x <- gsub("\\|", "/", x)
  a <- strsplit(x, "/", fixed = TRUE)[[1]]
  if (length(a) != 2 || any(a == ".")) return(NA_real_)
  vals <- suppressWarnings(as.numeric(a))
  if (any(is.na(vals))) return(NA_real_)
  s <- sum(vals)
  if (!s %in% 0:2) return(NA_real_)
  return(s)
}

# Apply dosage conversion to all SNPs and individuals
dos <- apply(gt, c(1, 2), gt_to_dosage)


# -----------------------------
# 8. Format SNP matrix for OutFLANK
# -----------------------------

# Transpose dosage matrix so rows = individuals and columns = SNPs
SNPmat <- t(dos)

# Replace missing genotypes with 9, which OutFLANK uses for missing data
SNPmat[is.na(SNPmat)] <- 9

# Make sure matrix is numeric
storage.mode(SNPmat) <- "numeric"

# Confirm dimensions match
stopifnot(length(popNames) == nrow(SNPmat))
stopifnot(length(locusNames) == ncol(SNPmat))


# -----------------------------
# 9. Calculate FST values
# -----------------------------

# Create the FST dataframe required by OutFLANK
FSTdf <- MakeDiploidFSTMat(
  SNPmat = SNPmat,
  locusNames = locusNames,
  popNames = popNames
)


# -----------------------------
# 10. Run OutFLANK
# -----------------------------

# Run OutFLANK using seven estuary populations
# Hmin = 0.1 removes very low heterozygosity SNPs
# qthreshold = 0.05 flags significant outliers after multiple-testing correction
out <- OutFLANK(
  FSTdf,
  NumberOfSamples = length(unique(popNames)),
  LeftTrimFraction = 0.05,
  RightTrimFraction = 0.05,
  Hmin = 0.1,
  qthreshold = 0.05
)


# -----------------------------
# 11. View and summarize results
# -----------------------------

# View first few rows of results
head(out$results)

# Count outlier and non-outlier SNPs
table(out$results$OutlierFlag, useNA = "ifany")

# Summarize q-values
summary(out$results$qvalues)

# Extract only significant outliers
outliers <- subset(out$results, OutlierFlag == TRUE)

# Sort outliers by q-value
outliers <- outliers[order(outliers$qvalues), ]

# View top outliers
head(outliers, 20)


# -----------------------------
# 12. Save output files
# -----------------------------

# Save all OutFLANK results
write.csv(
  out$results,
  "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/outflank/outflank_results_n218_geno075_clumped.csv",
  row.names = FALSE
)

# Save only significant outliers
write.csv(
  outliers,
  "/mnt/jupiter/johnsonlab/Capstone_proj/Capstone_2026/outflank/outflank_outliers_n218_geno075_clumped.csv",
  row.names = FALSE
)
