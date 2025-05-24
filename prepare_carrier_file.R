##### Prepare file with carriers of variants with relative frequency difference > 0.1 and at least 5 carriers
rm(list = ls());
library(stringr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(dplyr)
library(data.table)

############### Select inds from Beauce with variant ############################
variants_file <- "/lustre03/project/6033529/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/results_clinvar/variant_enriched_in_beauce_compared_UrbanQc_Imput_CaG_tresh_0.1_V2.tped"
indiv_file <- "/lustre03/project/6033529/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/results_clinvar/variant_enriched_in_beauce_compared_UrbanQc_Imput_CaG_tresh_0.1_V2.tfam"
# Read the .tped file and only keep enriched variants
tped_data <- fread(variants_file, header = FALSE)
# Identify the minor alleles
#VAT_variants_with_information_NTNFE <- read.table("/lustre03/project/6033529/schizo/mylgag/scripts/beauce_project/balsac/results_clinvar/carrier_rates_variant_enriched_in_beauce_compared_UrbanQc_Imput_CaG_tresh_0.1.txt", header = TRUE)
VAT_variants_with_information_NTNFE <- read.table("/lustre03/project/6033529/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/results_clinvar/all_variant_enriched_in_beauce_compared_UrbanQc_Imput_CaG_tresh_0.1_V2.txt", header = TRUE)
minor_alleles <- data.frame(SNP = paste0(VAT_variants_with_information_NTNFE$CHROM, ":", VAT_variants_with_information_NTNFE$POS), ALT = VAT_variants_with_information_NTNFE$A1)
minor_alleles <- unique(minor_alleles)
# Keep only variants in the list
tped_data <- subset(tped_data, V2 %in% minor_alleles$SNP)

# Extract the genotype columns (5th to last columns)
genotype_data <- tped_data[, -c(1:4)]
genotype_data <- cbind(tped_data$V2, genotype_data)
# Create new column names
indiv_data <- read.table(indiv_file)
new_colnames <- unlist(lapply(seq_along(indiv_data$V2), function(i) {
  c(paste0(indiv_data$V2[i], "_1"), paste0(indiv_data$V2[i], "_2"))
}))
names(genotype_data) <- c("SNP", new_colnames)

## Identify all alleles with ALT
df <- genotype_data %>%
  left_join(minor_alleles, by = "SNP")
matched_individuals <- df %>%
  rowwise() %>%
  mutate(matching_individuals = list(names(select(cur_data(), -SNP, -ALT))[select(cur_data(), -SNP, -ALT) == ALT])) %>%
  pull(matching_individuals)
names(matched_individuals) <- df$SNP
## Remove suffixes to keep only the IIDs carrying
matched_individuals_corrected <- lapply(matched_individuals, function(individuals) {
  sub("(_[0-9]+)$", "", individuals)
})
## Removed doubled IDs (homozygotes)
matched_individuals_corrected <- lapply(matched_individuals_corrected, unique)

## Keep only SNPs with 5 or more carriers
filtered_list <- Filter(function(x) length(x) >= 5, matched_individuals_corrected)
## Prepare the list of individuals
carriers <- data.frame(
  Pair = unlist(filtered_list)
)
variant <- row.names(carriers)
carriers <- cbind(variant, carriers)
names(carriers) <- c("variant", "IID")

good_variants <- minor_alleles$SNP
# Function to replace values based on the match with the original vector (numbers introduced after variant ID)
replace_values <- function(column, original_values) {
  sapply(column, function(x) {
    match <- original_values[str_detect(x, paste0("^", original_values))]
    if(length(match) > 0) {
      return(match[1])
    } else {
      return(x)
    }
  })
}
# Apply the function to the column in the data frame
carriers_corrected <- carriers
carriers_corrected$variant <- replace_values(carriers$variant, good_variants)
## Modify ID_allele to be ID_ID
# remove the allele identifier in ID
carriers_corrected$IID <- sub("_.*", "", carriers_corrected$IID)
# Double the ID to match IBD file
carriers_corrected <- carriers_corrected %>%
  mutate(across(IID, ~ ifelse(!is.na(.), paste0(., "_", .), NA)))
write.table(carriers_corrected, "/lustre03/project/6033529/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/carriers_enriched_variants_beauceVSurbanQc_V2.txt", sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

carrier_file <- "/lustre03/project/6033529/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/carriers_enriched_variants_beauceVSurbanQc_in_UrbanQc_V2.txt"
carriers_corrected <- read.table(carrier_file, header = TRUE)
path_out <- "/home/mylgag/projects/rrg-girardsi/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/results_clinvar/"
file_names <- readLines("/home/mylgag/projects/rrg-girardsi/schizo/mylgag/Beauce_founder_effect/enriched_variants/scripts/ibd_files_list.txt")

# Loop through each IBD file and chromosome number
chr_numbers <- 1:22
for (i in seq(1, length(file_names))) {
  file <- file_names[i]
  chr <- chr_numbers[i]  # Get the chromosome number for this file
  print(paste("Reading file:", file))

  # Read the file and name columns
  file_data <- read.table(file, header = TRUE)
  names(file_data) <- c('ind1', 'chr1', 'ind2', 'chr2', 'chr', 'start', 'end', 'LOD', 'length')

  # Filter rows based on the IBD length (length >= 2)
  file_data <- file_data[file_data$length >= 2, ]
  # Loop through each variant
  for (variant in unique(carriers_corrected$variant)) {
    # Get the IDs for the current variant
    ids <- carriers_corrected$IID[carriers_corrected$variant == variant]

    # Filter IBD sharing where both ID1 and ID2 are in the list of IDs for the current variant
    filtered <- file_data %>%
      filter(ind1 %in% ids & ind2 %in% ids) %>%
      mutate(variant = variant)

    # Define the output file path using the specified pattern
    output_file <- paste0(path_out, "ibd_sharing_carrier_BeauceVSUrbanQc_", variant, "_chr", chr, "_in_UrbanQc.txt")

    # Write the filtered data to the output file
    write.table(filtered, output_file, row.names = FALSE, quote = FALSE)

    print(paste("Written file:", output_file))
  }
}

## Save the list of variants to scroll through
unique_variants <- data.frame(variant = unique(carriers_corrected$variant))
write.table(unique_variants, paste0(path_out, "list_of_variants_BeauceVSUrbanQc_in_UrbanQc_V2.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

