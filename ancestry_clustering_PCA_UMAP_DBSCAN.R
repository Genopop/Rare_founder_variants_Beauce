#!/bash/env R

# Load required libraries
library(dplyr)
library(ggplot2)
library(umap)
library(dbscan)
library(gridExtra)
library(grid)
library(changepoint)

# =========================
# 1. Load PCA Object
# =========================
load("schizo.cag.hg38.commonsnps.mind0.05.hwe10-6_LDpruned_maf0.05_canadianIDs_PC-AiR_IBD0.125_pedrigree0.0625_related.RData")
mypcair <- unrel_pcair
pcs <- data.frame(mypcair$vectors)
colnames(pcs) <- paste0("PC", 1:ncol(pcs))
pcs <- pcs %>% tibble::rownames_to_column(var = "IID")
ids <- data.frame(IID = pcs$IID)

# =========================
# 2. Scree Plot with Change Points
# =========================
eigens <- mypcair$values
prop_var <- eigens / sum(eigens)
cgp <- cpt.meanvar(prop_var)@cpts
cgp <- cgp[cgp != max(cgp)]
eigens_df <- data.frame(V1 = eigens)
eigens_df <- (eigens_df / sum(eigens_df)) * 100

gg <- ggplot(eigens_df, aes(x = seq_along(V1), y = V1)) + 
  geom_point(size = 1.5) + 
  geom_line(linewidth = 1) +
  geom_vline(xintercept = cgp, linetype = "dashed", color = "red") +
  labs(x = "Principal component", y = "Explained variance (%)") +
  theme(text = element_text(family = "Arial"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

ggsave("figS2_PCs_eigenvals.svg", gg, width = 5, height = 5, dpi = 300, units = "in")

# =========================
# 3. Region and Phenotype Information
# =========================
# Load metadata
reg <- read.table("Regions_familles.csv", header = TRUE, sep = ";", fileEncoding = "UTF-8")
reg <- data.frame(FID = reg$Famille, region = reg$Region)

region_carta <- read.table("CaG_regions.txt", header = TRUE, fileEncoding = "UTF-8")
region_carta30k <- read.table("region_cartagene.30k.txt", header = FALSE, fileEncoding = "UTF-8")
region_carta30k <- data.frame(IID = region_carta30k$V2, region = region_carta30k$V3)

continents <- read.table("cartagene_some_continent.txt", header = TRUE, fileEncoding = "UTF-8")
country <- read.table("cartagene_birth_country.30k.txt", header = TRUE)

# Combine country and continent data
idscountry <- merge(ids, country, by.x = "IID", by.y = "file_111", all.x = TRUE) %>%
  mutate(COUNTRY_BIRTH = coalesce(COUNTRY_BIRTH, COUNTRY_BIRTH_OTHER),
         COUNTRY_BIRTH = if_else(COUNTRY_BIRTH == 99, COUNTRY_BIRTH_OTHER, COUNTRY_BIRTH)) %>%
  merge(continents, by.x = "COUNTRY_BIRTH", by.y = "country_num", all.x = TRUE) %>%
  select(-COUNTRY_BIRTH, -COUNTRY_BIRTH_OTHER, -continent_num)

# Load phenotype information
phen <- read.table("pheno14052021.csv", sep = ";", header = TRUE, fileEncoding = "UTF-8")
pheno_df <- phen %>% transmute(FID = fam, IID = id, BP = pheno_BPbr, SZ = pheno_SZbr)

# Classify individuals
pheno_df <- pheno_df %>%
  mutate(Category = case_when(
    BP == 2 & SZ == 1 ~ "Bipolar",
    BP == 1 & SZ == 2 ~ "Schizophrenic",
    BP == 1 & SZ == 1 ~ "Unaffected",
    BP == 0 & SZ == 0 ~ "Unknown",
    TRUE ~ NA_character_)) %>%
  drop_na(Category)

# Merge metadata
idscountry <- merge(idscountry, pheno_df[, c("FID", "IID", "Category")], by = "IID", all.x = TRUE)
idsreg <- idscountry %>%
  left_join(region_carta30k, by = "IID") %>%
  left_join(reg, by = "FID") %>%
  mutate(regions = coalesce(region, region.y)) %>%
  mutate(cohorte = case_when(
    IID %in% pheno_df$IID ~ "PDC",
    IID %in% region_carta30k$IID ~ "CaG",
    TRUE ~ NA_character_
  )) %>%
  mutate(origin = case_when(
    country == "NO" ~ NA_character_,
    country == "CANADA" | is.na(country) ~ regions,
    TRUE ~ country
  )) %>%
  mutate(origin = recode(origin,
    "HAITI" = "Haiti", "MOROCCO" = "Morocco", "NB_Iles" = "Qc (Iles-de-la-Madeleine)",
    "Beauce" = "Qc (Beauce)", "Trois-Rivières" = "Qc (Trois-Rivieres)",
    "Saguenay" = "Qc (SLSJ)", "Sherbrooke" = "Qc (Sherbrooke)", 
    "Montréal" = "Qc (Montreal)", "Gatineau" = "Qc (Gatineau)", 
    "Québec" = "Qc (Quebec City)"
  )) %>%
  mutate(origin = if_else(is.na(origin), continent, origin)) %>%
  mutate(origin = if_else(is.na(origin), "Unknown", origin)) %>%
  distinct()

# =========================
# 4. PCA Plotting
# =========================
pca_df <- merge(idsreg, pcs, by = "IID")
colors <- c("Qc (Montreal)" = "orange", "Qc (Beauce)" = "red", "Qc (SLSJ)" = "blue", 
            "Qc (Iles-de-la-Madeleine)" = "#00FF00", "Qc (Quebec City)" = "#FFD700",
            "Western Europe" = "#00FFFF", "Haiti" = "#FF00FF", "Morocco" = "#40E0D0",
            "Qc (Sherbrooke)" = "#8A2BE2", "Qc (Gatineau)" = "#006400",
            "South America" = "#FF7F50", "Asia" = "#87CEEB", "Eastern Europe" = "#be8ac2", 
            "Qc (Trois-Rivieres)" = "#800000", "Middle East" = "#FFE5B4", 
            "North America" = "#98FF98", "South Africa" = "#808000", 
            "North Africa" = "#DC143C", "Northern Europe" = "#4682B4", "Unknown" = "grey")

theme_custom <- theme(text = element_text(family = "Arial"),
                      axis.text = element_text(size = 8),
                      axis.title = element_text(size = 10),
                      legend.text = element_text(size = 8),
                      legend.title = element_text(size = 10),
                      legend.position = "none")

create_plot <- function(data, x, y, eigens, origin_colors) {
  ggplot(data, aes_string(x = x, y = y, color = "origin")) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = origin_colors, na.value = "grey") +
    labs(x = paste0(x, " (", round(eigens[as.numeric(gsub("PC", "", x))], 2), "%)"),
         y = paste0(y, " (", round(eigens[as.numeric(gsub("PC", "", y))], 2), "%)")) +
    guides(color = guide_legend(override.aes = list(size = 1.5))) +
    theme_custom
}

pc_plots <- lapply(2:12, function(i) create_plot(pca_df, "PC1", paste0("PC", i), eigens_df$V1, colors))

# Create grid of PCA plots
legend_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = origin)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(title = "Recruitment\nregion or birthplace", ncol = 2)) +
  theme(legend.position = "bottom")

legend_grob <- ggplotGrob(legend_plot)$grobs[[which(sapply(ggplotGrob(legend_plot)$grobs, function(x) x$name) == "guide-box")]]

empty_plot <- nullGrob()
gg <- grid.arrange(
  arrangeGrob(
    arrangeGrob(pc_plots[[1]], pc_plots[[4]], pc_plots[[7]], pc_plots[[10]], ncol = 1),
    arrangeGrob(pc_plots[[2]], pc_plots[[5]], pc_plots[[8]], pc_plots[[11]], ncol = 1),
    arrangeGrob(pc_plots[[3]], pc_plots[[6]], pc_plots[[9]], empty_plot, ncol = 1),
    ncol = 3
  ),
  legend_grob,
  nrow = 2,
  heights = c(10, 1.5)
)
ggsave("figS3_PCA_12PCs.png", gg, width = 8.27, height = 9.27, dpi = 300, units = "in")

# =========================
# 5. UMAP and DBSCAN
# =========================

# UMAP projection
set.seed(42)
umap_result <- umap(pcs[, paste0("PC", 1:12)])
umap_coords <- as.data.frame(umap_result$layout)
colnames(umap_coords) <- c("Coord1", "Coord2")

# Merge with metadata
infos <- data.frame(IID = pca_df$IID, regions = pca_df$regions, cohorte = pca_df$cohorte, origin = pca_df$origin)
umap_df <- cbind(umapdf, infos)

# UMAP plot by origin
gg1 <- ggplot(umap_df, aes(x = Coord1, y = Coord2, color = origin, alpha = 0.5)) +
  geom_point(size = 0.75) +
  scale_color_manual(values = colors) +
  scale_alpha(guide = 'none') +
  labs(x = "Coord1", y = "Coord2", color = "Recruitment region or birthplace") +
  guides(color = guide_legend(title.position = "top", override.aes = list(size = 1.5), ncol = 3)) +
  theme_custom +
  labs(tag = "a") +
  coord_fixed()

dbscan_result <- dbscan(umap_df[, c("Coord1", "Coord2")], eps = 0.3, minPts = 50)
umap_df$cluster <- dbscan_result$cluster
umap_df$IID <- pca_df$IID
# Remove outliers
no0 <- subset(umap_df, cluster != 0)
# UMAP plot by cluster
gg2 <- ggplot(no0, aes(x = Coord1, y = Coord2, color = as.factor(cluster), alpha = 0.5)) +
  geom_point(size = 0.75) +
  scale_color_manual(values = colors2) +
  scale_alpha(guide = 'none') +
  labs(x = "Coord1", y = "Coord2", color = "DBSCAN cluster") +
  guides(color = guide_legend(title.position = "top", override.aes = list(size = 1.5), ncol = 7)) +
  theme_custom +
  labs(tag = "b") +
  coord_fixed()

gg <- grid.arrange(gg1, gg2, ncol = 2)
ggsave("fig5_12PCs_UMAP.jpeg",
       gg, width = 8.27, height = 5.9, dpi = 300, units = "in")

# Identify individuals from Beauce, UrbanQC and SLSJ clusters
get_cluster_inds <- function(df, cluster_number) {
  clust <- subset(df, cluster == cluster_number)
  data.frame(FID = clust$IID, IID = clust$IID)
}
write.table(get_cluster_inds(no0, 1),
            "Beauce_cluster.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(get_cluster_inds(no0, 5),
            "UrbanQc_cluster.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(get_cluster_inds(no0, 2),
            "SLSJ_cluster.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)                                                           

# Get proportion of individuals per region in the clusters
clusters <- subset(no0, cluster %in% c(1, 5))
clusters <- merge(clusters, infos, by = "IID")

cluster_stats <- clusters %>%
  group_by(cluster, origin) %>%
  summarise(nind = n(), .groups = "drop") %>%
  mutate(
    proportion = case_when(
      cluster == 1 ~ nind / 1362 * 100,
      cluster == 5 ~ nind / 11998 * 100
    ),
    cluster = case_when(
      cluster == 1 ~ "Beauce",
      cluster == 5 ~ "UrbanQc"
    )
  )
gg <- ggplot(cluster_stats, aes(x = cluster, y = proportion, fill = origin)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = colors) +
  labs(x = "Cluster", y = "Percentage of individuals (%)", fill = "Recruitment region \nor birthplace") +
  theme_custom

ggsave("figS4_cluster_regional_composition.svg", 
       gg, width = 8, height = 6, dpi = 300, units = "in")                                                           
                                                           
