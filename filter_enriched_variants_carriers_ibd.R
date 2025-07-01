##### Identify and extract individuals from Beauce carrying variants with:
# - Relative frequency difference > 0.1
# - At least 5 carriers

rm(list = ls())
# Load required libraries
library(stringr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(dplyr)
library(data.table)

# Paths to .tped and .tfam files
variants_file <- "variant_enriched_in_beauce_compared_UrbanQc_tresh_0.1.tped"
indiv_file <- "variant_enriched_in_beauce_compared_UrbanQc_tresh_0.1.tfam"

# Read enriched variants
tped_data <- fread(variants_file, header = FALSE)

# Read variant annotations to identify ALT alleles
variants_with_information <- read.table(
  "all_variant_enriched_in_beauce_compared_UrbanQc_tresh_0.1.txt", 
  header = TRUE
)  # Output from identify_variants_with_higher_freq

# Build dataframe of ALT alleles
minor_alleles <- data.frame(
  SNP = paste0(variants_with_information$CHROM, ":", variants_with_information$POS),
  ALT = variants_with_information$A1
) %>% unique()

# Keep only variants of interest
tped_data <- subset(tped_data, V2 %in% minor_alleles$SNP)

##### Prepare genotype data
# Extract genotype columns
genotype_data <- tped_data[, -c(1:4)]
genotype_data <- cbind(tped_data$V2, genotype_data)

# Assign individual-level column names
indiv_data <- read.table(indiv_file)
new_colnames <- unlist(lapply(seq_along(indiv_data$V2), function(i) {
  c(paste0(indiv_data$V2[i], "_1"), paste0(indiv_data$V2[i], "_2"))
}))
names(genotype_data) <- c("SNP", new_colnames)

##### Match individuals carrying ALT alleles

df <- genotype_data %>%
  left_join(minor_alleles, by = "SNP")

# Identify carriers of each ALT
matched_individuals <- df %>%
  rowwise() %>%
  mutate(matching_individuals = list(names(select(cur_data(), -SNP, -ALT))[select(cur_data(), -SNP, -ALT) == ALT])) %>%
  pull(matching_individuals)

names(matched_individuals) <- df$SNP

# Clean individual IDs: remove allele suffixes and deduplicate homozygotes
matched_individuals_corrected <- lapply(matched_individuals, function(individuals) {
  sub("(_[0-9]+)$", "", individuals)
}) %>% lapply(unique)

# Keep only variants with at least 5 carriers
filtered_list <- Filter(function(x) length(x) >= 5, matched_individuals_corrected)

# Prepare carriers dataframe
carriers <- data.frame(Pair = unlist(filtered_list))
variant <- row.names(carriers)
carriers <- cbind(variant, carriers)
names(carriers) <- c("variant", "IID")

##### Clean and reformat variant and individual IDs

good_variants <- minor_alleles$SNP

# Function to clean variant names
replace_values <- function(column, original_values) {
  sapply(column, function(x) {
    match <- original_values[str_detect(x, paste0("^", original_values))]
    if (length(match) > 0) return(match[1]) else return(x)
  })
}

# Apply cleaning
carriers_corrected <- carriers
carriers_corrected$variant <- replace_values(carriers$variant, good_variants)
carriers_corrected$IID <- sub("_.*", "", carriers_corrected$IID)
carriers_corrected <- carriers_corrected %>%
  mutate(across(IID, ~ ifelse(!is.na(.), paste0(., "_", .), NA)))

# Save carriers file
write.table(
  carriers_corrected, 
  "carriers_enriched_variants_beauceVSurbanQc.txt", 
  sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE
)

##### Load final carrier list and IBD file list
carrier_file <- "carriers_enriched_variants_beauceVSurbanQc.txt"
carriers_corrected <- read.table(carrier_file, header = TRUE)

file_names <- readLines("ibd_files_list.txt")

##### Filter IBD segments shared by carriers for each chromosome
chr_numbers <- 1:22

for (i in seq_along(file_names)) {
  file <- file_names[i]
  chr <- chr_numbers[i]
  print(paste("Reading file:", file))

  # Load and filter IBD segments by length
  file_data <- read.table(file, header = TRUE)
  names(file_data) <- c('ind1', 'chr1', 'ind2', 'chr2', 'chr', 'start', 'end', 'LOD', 'length')
  file_data <- file_data[file_data$length >= 2, ]

  # Process each variant
  for (variant in unique(carriers_corrected$variant)) {
    ids <- carriers_corrected$IID[carriers_corrected$variant == variant]

    filtered <- file_data %>%
      filter(ind1 %in% ids & ind2 %in% ids) %>%
      mutate(variant = variant)

    output_file <- paste0(path_out, "ibd_sharing_carrier_BeauceVSUrbanQc_", variant, "_chr", chr, ".txt")
    write.table(filtered, output_file, row.names = FALSE, quote = FALSE)
    print(paste("Written file:", output_file))
  }
}
                        
##### Save final list of variants
unique_variants <- data.frame(variant = unique(carriers_corrected$variant))
write.table(
  unique_variants,
  paste0(path_out, "list_of_variants_BeauceVSUrbanQc.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)
