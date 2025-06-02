#############
# File paths
#############
beauce_file <- "clusterbeauce.hwe"
uqc_file <- "clusteruqc.hwe"
variants_file <- "all_variant_enriched_in_beauce_compared_UrbanQc_tresh_0.1.txt"  # Output from identify_variants_with_higher_freq.R
#############
# Read input data
#############
beauce <- read.table(beauce_file, header = TRUE)
uqc <- read.table(uqc_file, header = TRUE)
variants <- read.table(variants_file, header = TRUE)
# Extract SNP identifiers and reformat to match genotype data
SNPs <- variants$SNP
variants$SNP2 <- sub("^chr(\\d+)_(\\d+)_.*$", "\\1:\\2", variants$SNP)

#############
# Process genotype data - Beauce cohort
#############
# Extract and split genotype counts
hw_beauce <- strsplit(beauce$GENO, "/")
hw_beauce <- do.call(rbind, lapply(hw_beauce, as.numeric))
hw_beauce <- data.frame(SNP = beauce$SNP, homo1 = hw_beauce[,1], hetero = hw_beauce[,2], homo2 = hw_beauce[,3])
# Keep only enriched variants
hw_beauce <- subset(hw_beauce, SNP %in% variants$SNP2)
# Calculate carrier rate: CR = 1 / (heterozygotes / n)
hw_beauce$CR <- 1 / (hw_beauce$hetero / 317)
# Keep relevant columns
hw_beauce <- data.frame(SNP = hw_beauce$SNP, CR_beauce = hw_beauce$CR)

#############
# Process genotype data - UrbanQc cohort
#############
# Extract and split genotype counts
hw_uqc <- strsplit(uqc$GENO, "/")
hw_uqc <- do.call(rbind, lapply(hw_uqc, as.numeric))
hw_uqc <- data.frame(SNP = uqc$SNP, homo1 = hw_uqc[,1], hetero = hw_uqc[,2], homo2 = hw_uqc[,3])
# Keep only enriched variants
hw_uqc <- subset(hw_uqc, SNP %in% variants$SNP2)
# Calculate carrier rate: CR = 1 / (heterozygotes / n)
hw_uqc$CR <- 1 / (hw_uqc$hetero / 893)
# Keep relevant columns
hw_uqc <- data.frame(SNP = hw_uqc$SNP, CR_UrbanQc = hw_uqc$CR)

#############
# Merge carrier rate data with enriched variants table
#############
# Remove existing CR columns if present
variants$CR_beauce <- NULL
variants$CR_UQC <- NULL
# Merge with UrbanQc data
variants <- merge(variants, hw_uqc, by.x = "SNP2", by.y = "SNP")
# Merge with Beauce data
variants <- merge(variants, hw_beauce, by.x = "SNP2", by.y = "SNP")
# Clean up temporary column
variants$SNP2 <- NULL

write.table(variants, "final_variant_results_with_carrier_rates_beauceVSUrbanQc.txt", sep = ";", col.names = TRUE, row.names = FALSE)


