#### Identify clinvar variants with a higher frequency in Beauce compared to UrbanQc
# Load ClinVar reference variants
variants_info <- read.table("ClinVar_variants_to_analyse_well_ref_and_known_variants_added.txt", header = TRUE)
variants_info$variant <- paste0(variants_info$CHROM, "_", variants_info$POS)

# Load Beauce and UrbanQc allele frequency data (output from PLINK --freq)
beauce <- read.table(subpop1, header = TRUE)
UrbanQc <- read.table(UrbanQc, header = TRUE)

# -----------------------------------------------------------------------------
# Format SNP column and extract relevant fields
# -----------------------------------------------------------------------------
# Extract POS from the SNP column
beauce$POS <- str_split_fixed(beauce$SNP, ":", 4)[, 2]
UrbanQc$POS <- str_split_fixed(UrbanQc$SNP, ":", 4)[, 2]

# Reconstruct SNP ID in the desired format
beauce$SNP <- paste0("chr", beauce$CHR, "_", beauce$POS, "_", beauce$A2, "_", beauce$A1)
UrbanQc$SNP <- paste0("chr", UrbanQc$CHR, "_", UrbanQc$POS, "_", UrbanQc$A2, "_", UrbanQc$A1)

# Keep only columns of interest: SNP, MAF, NCHROBS
beauce <- beauce[, -c(1, 3, 4, 7)]
UrbanQc <- UrbanQc[, -c(1, 3, 4, 7)]

# -----------------------------------------------------------------------------
# Merge datasets and compute relative frequency difference
# -----------------------------------------------------------------------------
comp <- merge(beauce, UrbanQc, by = "SNP")
names(comp) <- c("SNP", "MAF_beauce", "NCHROBS_beauce", "MAF_UrbanQc", "NCHROBS_UrbanQc")

# Compute standardized enrichment: relative difference in MAF
m <- nrow(comp)
for (i in 1:m){
  comp$enrich_stand_beauce_UrbanQc[i] = (as.numeric(comp$MAF_beauce[i]) - as.numeric(comp$MAF_UrbanQc[i])) / as.numeric(comp$MAF_beauce[i])
}

# Filter variants enriched in Beauce (≥ 0.1) and remove incomplete rows
comp_clean <- comp[comp$enrich_stand_beauce_UrbanQc >= 0.1, ]
comp_clean <- comp_clean[complete.cases(comp_clean), ]

# -----------------------------------------------------------------------------
# Extract variant info and prepare for merge with ClinVar
# -----------------------------------------------------------------------------
# Decompose SNP into CHROM, POS, A1, A2, and recompute 'variant'
comp_clean$CHROM <- sub("^chr([0-9]+)_.*", "\\1", comp_clean$SNP)
comp_clean$POS <- sub("^chr[0-9]+_([0-9]+)_.*", "\\1", comp_clean$SNP)
comp_clean$A2 <- sub("^.*_([A-Z]+)_.*", "\\1", comp_clean$SNP)
comp_clean$A1 <- sub("^.*_([A-Z]+)$", "\\1", comp_clean$SNP)
comp_clean$variant <- paste0(comp_clean$CHROM, "_", comp_clean$POS)

# -----------------------------------------------------------------------------
# Merge with ClinVar data, accounting for possible allele switches
# -----------------------------------------------------------------------------
switch_handler <- function(comp_clean, variants_info) {
  comp_clean %>%
    left_join(variants_info, by = "variant") %>%
    mutate(
      matched = case_when(
        A1 == REF & A2 == ALT ~ TRUE,  # Switched
        A1 == ALT & A2 == REF ~ TRUE   # Exact match
      )
    ) %>%
    filter(matched == TRUE)  # Keep only valid matches
}

new_beauce_df <- switch_handler(comp_clean, variants_info)

# -----------------------------------------------------------------------------
# Final dataset with relevant annotations
# -----------------------------------------------------------------------------
final_beauce_df <- data.frame(
  SNP = new_beauce_df$SNP,
  CHROM = new_beauce_df$'CHROM.x',
  POS = new_beauce_df$'POS.x',
  ID = new_beauce_df$ID,
  REF = new_beauce_df$REF,
  ALT = new_beauce_df$ALT,
  CLNDN = new_beauce_df$CLNDN,
  MAF_beauce = new_beauce_df$MAF_beauce,
  MAF_UrbanQc = new_beauce_df$MAF_UrbanQc,
  NCHROBs_beauce = new_beauce_df$NCHROBS_beauce,
  NCHROBs_UrbanQc = new_beauce_df$NCHROBS_UrbanQc,
  enrich_stand_beauce_UrbanQc = new_beauce_df$enrich_stand_beauce_UrbanQc,
  A1 = new_beauce_df$A1,
  A2 = new_beauce_df$A2
)
# Remove duplicate entries
final_beauce_df <- unique(final_beauce_df)

# -----------------------------------------------------------------------------
# Plotting: Allele frequencies comparison between UrbanQc and Beauce
# -----------------------------------------------------------------------------
p <- ggplot(final_beauce_df, aes(x = MAF_UrbanQc, y = MAF_beauce)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    title = "Comparison of Allele Frequencies: Beauce cluster vs UrbanQc cluster",
    x = "Allele Frequency in UrbanQc cluster",
    y = "Allele Frequency in Beauce cluster"
  ) +
  theme(
    plot.title = element_text(size = 24),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
# Save plot to file
ggsave("differential_frequency_BeauceVSUrbanQc.png", 
       p, width = 10, height = 10, units = "in")

# -----------------------------------------------------------------------------
# Output result tables
# -----------------------------------------------------------------------------
# Save final annotated variant list
write.table(final_beauce_df,
            "all_variant_enriched_in_beauce_compared_UrbanQc_tresh_0.1.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Prepare simple variant ID list (CHR:POS) for extraction
var <- paste0(final_beauce_df$CHR, ":", final_beauce_df$POS)
var <- unique(data.frame(var = var))

# Save variant ID list
write.table(var,
            "to_extract_variant_enriched_in_beauce_compared_UrbanQc_Imput_CaG_tresh_0.1.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

