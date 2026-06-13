# Sea-Grant-Capstone2026
Repository for the Sea Grant California Johnson Lab Bioinformatics Capstone 2026 Team 

Contact information:
- Jack Poole: 
jbpoole@calpoly.edu, jackbpoole@gmail.com, +1 (503)-858-0105

Repository Structure:

- SNP_Preprocessing Folder — Folder for required files to get SNP data ready for outlier detection
    - Preprocessing Steps.Rmd — Markdown document for running all necessary filtering steps
    - filter_hwe_by_pop.pl — dDocent's script for analyzing Hardy-Weinberg equilibrium of SNPs within each pop
    - n218_popmap.txt — population map for the n218 vcf of genotyped oysters. Needed for multiple sections of Preprocessing Steps.Rmd
    - pop_missing_filter.sh — dDocent's script for removing SNPs with low call rates within each population
- Observed Heterozygosity and ANOVA Analysis.md — script for calculating heterozygosity within each population
- pcadapt Steps.Rmd — PCA-based genome scan for candidate loci associated with population structure.
- OutFLANK Steps.Rmd — FST-based genome scan used to identify loci showing elevated genetic differentiation among estuary populations.

Overview--
  This project was completed as part of the California Sea Grant Bioinformatics Capstone. We analyzed single nucleotide polymorphism (SNP) data from Olympia oysters (Ostrea lurida) collected across         multiple California estuaries to investigate genetic diversity, population structure, and potential adaptive variation.
  
Background--
  The Olympia oyster is the only oyster species native to the U.S. West Coast. Historical overharvesting, habitat degradation, invasive species, and changing environmental conditions have contributed to    population declines. Restoration programs rely heavily on hatchery-raised oysters, making it important to understand how their genetic diversity compares to wild populations.
  
Objectives--
  Characterize genetic variation among California estuary populations
  Identify loci potentially under selection
  Detect population differentiation using SNP data
  
Methods--
  Data Processing--
    Quality filtering of SNP datasets
    Linkage disequilibrium (LD) pruning/clumping
    Population assignment using estuary metadata

Dataset--
  Genomic data were collected from 218 Olympia oysters sampled across multiple California estuaries. Following quality control filtering and linkage disequilibrium pruning, the resulting SNP dataset was    used for population genomic analyses.

pcadapt--
  Used to identify candidate loci under selection through principal component analysis-based genome scans.
  
OutFLANK--
  Used to identify SNPs with unusually high FST values compared to the neutral genome-wide distribution.

Results--
  Multiple complementary approaches were used to identify candidate loci associated with population differentiation among California estuaries. Genome scan methods detected a subset of SNPs showing         elevated differentiation relative to neutral expectations, suggesting potential signatures of local adaptation. Candidate loci identified across methods were prioritized for downstream analyses and       marker development. The resulting set of SNP markers provides a foundation for developing a targeted GTseq panel that can be used to assess genetic diversity, population structure, and connectivity       among Olympia oyster populations. These tools will support future conservation, restoration, and hatchery management efforts.
