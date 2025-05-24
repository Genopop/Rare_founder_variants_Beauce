#### PC-AiR Analysis on LD-pruned Genotype Data
# Load required libraries
library(GWASTools)
library(SNPRelate)
library(GENESIS)

# -----------------------------------------------------------------------------
# Step 1: Convert PLINK files to GDS format
# -----------------------------------------------------------------------------
snpgdsBED2GDS(
  bed.fn = paste0(dir, prefix, ".bed"),
  fam.fn = paste0(dir, prefix, ".fam"),
  bim.fn = paste0(dir, prefix, ".bim"),
  out.gdsfn = paste0(dir, prefix, "GenomicDataStructures.gds"),
  verbose = TRUE
)

# -----------------------------------------------------------------------------
# Step 2: Load GDS file as a GenotypeData object
# -----------------------------------------------------------------------------
geno <- GdsGenotypeReader(filename = paste0(dir, prefix, "GenomicDataStructures.gds"))
genoData <- GenotypeData(geno)
iids <- getScanID(genoData)

# -----------------------------------------------------------------------------
# Step 3: Read all individual IDs from PLINK .fam file
# -----------------------------------------------------------------------------
all_inds <- read.table(paste0(dir, prefix, ".fam"))

# -----------------------------------------------------------------------------
# Step 4: Filter out related individuals
# -----------------------------------------------------------------------------
## refined IBD (cutoff 0.125) + pedigree (cutoff 0.0625)
to_exclude <- "/home/mylgag/projects/rrg-girardsi/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/IBD_related_groups0.125_and_pedigree_groups0.0625_inds_to_exclude.txt"
to_exclude <- read.table(to_exclude, header = FALSE)

# Keep only unrelated individuals
unrel <- subset(all_inds, !(V2 %in% to_exclude$V1))
unrel_ids <- unrel$V2

# -----------------------------------------------------------------------------
# Step 5: Perform PC-AiR analysis
# -----------------------------------------------------------------------------
unrel_pcair <- pcair(
  genoData,
  sample.include = iids,
  unrel.set = unrel_ids
)

# Summarize PC-AiR object
print("Summarise pcair object:")
summary(unrel_pcair)

# Save PC-AiR object
save(
  unrel_pcair,
  file = paste0(dir, prefix, "_PC-AiR_IBD0.125_pedrigree0.0625_related.RData")
)

