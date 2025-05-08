subpop <- 
wholepop <- 

## Open ClinVar file
VAT_variants_info <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/ClinVar_reference/ClinVar_variants_to_analyse_well_ref_and_known_variants_added.txt", header = TRUE)
VAT_variants_info$variant <- paste0(VAT_variants_info$CHROM, "_", VAT_variants_info$POS)

beauce <- read.table(subpop1, header = T)
wholepop <- read.table(wholepop, header = T)
## Put right format SNP
# Extract the POS column from SNP
beauce$POS <-str_split_fixed(beauce$SNP, ":",4)[,2]
wholepop$POS <-str_split_fixed(wholepop$SNP, ":",4)[,2]
# Reformat the SNP column
beauce$SNP <- paste0("chr",beauce$CHR,"_",beauce$POS,"_",beauce$A2,"_",beauce$A1)
wholepop$SNP <- paste0("chr",wholepop$CHR,"_",wholepop$POS,"_",wholepop$A2,"_",wholepop$A1)
# Keep only SNP, MAF, and NCHROBS
beauce <- beauce[,-c(1,3,4,7)]
wholepop <- wholepop[,-c(1,3,4,7)]

## Merge by SNPs
comp <- merge(beauce, wholepop, by = "SNP")
names(comp) <- c("SNP", "MAF_beauce", "NCHROBS_beauce", "MAF_UrbanQc", "NCHROBS_UrbanQc")

# Calculate differential frequency
m <- nrow(comp)
for (i in 1:m){
  comp$enrich_stand_beauce_UrbanQc[i]= (as.numeric(comp$MAF_beauce[i]) - as.numeric(comp$MAF_UrbanQc[i]))/ as.numeric(comp$MAF_beauce[i])
}

# Apply filter (keep rows where the enrichment is at least 0.1 and then removes any rows with missing values)
comp_clean <- comp[comp$enrich_stand_beauce_UrbanQc >= 0.1, ]
comp_clean <- comp_clean[complete.cases(comp_clean), ]

## Merge with clinvar infos
comp_clean$CHROM <- sub("^chr([0-9]+)_.*", "\\1", comp_clean$SNP)
comp_clean$POS <- sub("^chr[0-9]+_([0-9]+)_.*", "\\1", comp_clean$SNP)
comp_clean$A2 <- sub("^.*_([A-Z]+)_.*", "\\1", comp_clean$SNP)
comp_clean$A1 <- sub("^.*_([A-Z]+)$", "\\1", comp_clean$SNP)
comp_clean$variant <- paste0(comp_clean$CHROM, "_", comp_clean$POS)

## Keep switches while merging
switch_handler <- function(comp_clean, VAT_variants_info) {
  comp_clean %>%
    left_join(VAT_variants_info, by = "variant") %>%
    mutate(
      # Check for exact match, and switches
      matched = case_when(
        A1 == REF & A2 == ALT ~ TRUE,  # Switched
        A1 == ALT & A2 == REF ~ TRUE,  # Exact match
      )
    ) %>%
    filter(matched == TRUE)  # Keep only matching rows
}
new_beauce_df <- switch_handler(comp_clean, VAT_variants_info)

final_beauce_df <- data.frame(SNP = new_beauce_df$SNP, CHROM = new_beauce_df$'CHROM.x', POS = new_beauce_df$'POS.x',
                            ID = new_beauce_df$ID, REF = new_beauce_df$REF, ALT = new_beauce_df$ALT,
                            CLNDN = new_beauce_df$CLNDN, MAF_beauce = new_beauce_df$MAF_beauce, MAF_UrbanQc = new_beauce_df$MAF_UrbanQc,
                            NCHROBs_beauce = new_beauce_df$NCHROBS_beauce, NCHROBs_UrbanQc = new_beauce_df$NCHROBS_UrbanQc,
                            enrich_stand_beauce_UrbanQc = new_beauce_df$enrich_stand_beauce_UrbanQc, A1 = new_beauce_df$A1, A2 = new_beauce_df$A2)

final_beauce_df <- unique(final_beauce_df)

final_beauce_df <- unique(final_beauce_df)
### graph freq vs
p <- ggplot(final_beauce_df, aes(x=MAF_UrbanQc, y=MAF_beauce)) +
  geom_point(size = 3) +
  geom_abline(slope=1, intercept=0, color="red") +
  labs(title=paste("Comparison of Allele Frequencies: Beauce cluster vs UrbanQc cluster"),
       x=paste("Allele Frequency in UrbanQc cluster"),
       y="Allele Frequency in Beauce cluster") +
  theme(
    plot.title = element_text(size = 24),  # Increase the size of the plot title
    axis.title.x = element_text(size = 20),  # Increase the size of the x-axis title
    axis.title.y = element_text(size = 20)  # Increase the size of the y-axis title
  )
ggsave("/lustre03/project/6033529/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/results_clinvar/differential_frequency_BeauceVSUrbanQcV2.png", p, width = 10, height = 10, units = "in")
write.table(final_beauce_df, "/lustre03/project/6033529/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/results_clinvar/all_variant_enriched_in_beauce_compared_UrbanQc_Imput_CaG_tresh_0.1_V2.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

var <- paste0(final_beauce_df$CHR,":",final_beauce_df$POS)
var <- data.frame(var = var)
var <- unique(var)
write.table(var, "/lustre03/project/6033529/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/results_clinvar/to_extract_variant_enriched_in_beauce_compared_UrbanQc_Imput_CaG_tresh_0.1_V2.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)




