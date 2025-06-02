################################
# Classify variants with a higher frequency in Beauce
# Familial : < 5 unrelated carriers
# Founder : >= 5 unrelated carriers + IBD sharing in >= 50% of pairs
# Multiple introductions : Neither familial nor founder
################################
### Load libraries
library('ggplot2')
library("data.table")
library("dplyr")
library("stringr")

# Load sharing proportions data (nind = # pairs sharing the variant;  # Output from IBDsharing_variantPosition.py
# ncarriers = # of carriers; proportion = nind / (carriers * (carriers - 1) / 2))
load(paste0(path1, "/IDB_sharing_proportion_enriched_variants_BeauceVSUrbanQc.RData"))
# Load variant carriers table
carriers <- read.table("carriers_enriched_variants_beauceVSurbanQc.txt", header = TRUE)

### Update variant format and extract chromosomal position
updated_list <- lapply(graph_list, function(df) {
  df$variant <- gsub(":", "_", df$variant)
  var <- do.call(rbind, strsplit(df$variant, "_", fixed = TRUE))
  var <- data.frame(var, stringsAsFactors = FALSE)
  colnames(var) <- c("chr_var", "pos_var")
  df <- cbind(var, df)
  return(df)
})

# Keep only rows with highest sharing proportion and largest distance
updated_list <- lapply(updated_list, function(df) {
  df <- subset(df, as.character(df$chr_var) == as.character(chr))
  df$distance <- abs(as.numeric(df$pos_var) - as.numeric(df$pos))
  df <- subset(df, proportion == max(proportion))
  df <- subset(df, distance == max(distance))
  return(df)
})
variant_sharing <- do.call(rbind, updated_list)


## Generate all possible carrier pairs per variant
pairs <- carriers %>%
  group_by(variant) %>%
  summarise(pairs = list(as.data.frame(t(combn(IID, 2))))) %>%
  unnest(pairs) %>%
  rename(ind1 = V1, ind2 = V2)
pairs$variant <- gsub(":", "_", pairs$variant)

### Load related groups data (based on IBD ≥ 0.125 or kinship ≥ 0.0625)
related_file <- "related_groups_refinedIBD0.125_kinship0.0625.csv"  # Output from identify_related_individuals.py
related <- readLines(related_file)
related_list <- lapply(strsplit(related, ","), function(group) {
  sapply(trimws(group), function(id) paste0(id, "_", id))
})

### Determine if each pair is related
related_pairs <- function(ind1, ind2, groupes) {
  any(sapply(groupes, function(groupe) {
    ind1 %in% groupe && ind2 %in% groupe
  }))
}
pairs$related <- mapply(related_pairs, pairs$ind1, pairs$ind2, MoreArgs = list(groupes = related_list))

### Separate related and unrelated individuals
# Unrelated
unrelated_df <- pairs %>% filter(!related)
unrelated_individuals <- unrelated_df %>%
  select(variant, ind1, ind2) %>%
  pivot_longer(cols = c(ind1, ind2), names_to = "role", values_to = "individual")
unrelated_counts <- unrelated_individuals %>%
  distinct(variant, individual) %>%
  dplyr::group_by(variant) %>%
  dplyr::summarise(unrelated_individuals_count = n(), .groups = "drop")
# Related
related_df <- pairs %>% filter(related)
related_individuals <- related_df %>%
  select(variant, ind1, ind2) %>%
  pivot_longer(cols = c(ind1, ind2), names_to = "role", values_to = "individual")
related_counts <- related_individuals %>%
  distinct(variant, individual) %>%
  dplyr::group_by(variant) %>%
  dplyr::summarise(related_individuals_count = n(), .groups = "drop")
# Merge related/unrelated counts
relationship_count <- merge(unrelated_counts, related_counts, by = "variant", all.x = TRUE, all.y = TRUE)
relationship_count <- relationship_count %>%
  replace_na(list(unrelated_individuals_count = 0, related_individuals_count = 0))

### Classify variants by status
# Construct base variant table
all_variant_df <- data.frame(
  variant = variant_sharing$variant,
  A1 = variant_sharing$A1,
  A2 = variant_sharing$A2,
  proportion = variant_sharing$proportion,
  ncarriers = variant_sharing$ncarriers
)
# Identify familial variants (less than 5 unrelated carriers)
familial_variants <- subset(relationship_count, unrelated_individuals_count < 5)
familial_variants <- merge(familial_variants, all_variant_df, by = "variant")
# Identify founder variants (proportion ≥ 0.5)
non_fam_variants <- subset(all_variant_df, !(variant %in% familial_variants$variant))
non_fam_founder <- subset(non_fam_variants, proportion >= 0.5)
non_fam_founder <- merge(non_fam_founder, relationship_count, by = "variant")
non_fam_founder <- unique(non_fam_founder)
# Identify multiple introductions variants (proportion < 0.5)
multi_intro_variants <- subset(non_fam_variants, proportion < 0.5)
multi_intro_variants <- merge(multi_intro_variants, relationship_count, by = "variant")
multi_intro_variants <- unique(multi_intro_variants)

### Label and combine all variants
familial_variants$category <- "Familial"
non_fam_founder$category <- "Founder"
multi_intro_variants$category <- "Multiple introductions"
variants <- rbind(non_fam_founder, multi_intro_variants, familial_variants)

# Final variant summary table
variants <- data.frame(
  variant = variants$variant,
  ncarriers = variants$ncarriers,
  proportion = variants$proportion,
  category = variants$category
)

### Add carrier rate and frequency info
cr_file <- "final_variant_results_with_carrier_rates_beauceVSUrbanQc.txt"
cr <- read.table(cr_file, header = TRUE, sep = ";")
cr$variant <- paste0(cr$CHROM, "_",  cr$POS)

var_info <- merge(variants, cr, by = "variant")
var_info <- subset(var_info, ALT == A1)

# Final data frame for export
final_df <- data.frame(
  SNP = var_info$SNP,
  patho = var_info$CLNDN,
  clinvarID = var_info$ID,
  MAF_Beauce = var_info$MAF_beauce,
  MAF_UQC = var_info$MAF_UrbanQc,
  CR_beauce = var_info$CR_beauce,
  CR_UrbanQc = var_info$CR_UrbanQc,
  statut = var_info$category,
  sharing_proportion = var_info$proportion,
  rfd = var_info$enrich_stand_beauce_UrbanQc
)

outfile <- "final_table_all_variants_with_info.txt"
write.table(final_df, outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ";")
