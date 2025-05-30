### Clean environment and load libraries
rm(list = ls())
library('ggplot2')
library("data.table")
library("dplyr")
library("stringr")
library("tidyr")
library("gridExtra")

### Load IBD sharing results (by position) from files
files <- list.files(path_in, pattern = "ibd_sharing_carrier_BeauceVSUrbanQc_.*_chr.*\\.txt.sharing.by.pos", full.names = TRUE)

# Extract unique variant identifiers and chromosomes from filenames
vars <- unique(gsub("ibd_sharing_carrier_BeauceVSUrbanQc_(.*)_chr.*\\.txt.sharing.by.pos", "\\1", basename(files)))
chr <- unique(gsub("ibd_sharing_carrier_BeauceVSUrbanQc_.*_chr(.*).txt.sharing.by.pos", "\\1", basename(files)))

results_list <- list()
# Load and combine sharing data by variant and chromosome
for (v in vars) {
  chr_list <- list()
  for (c in chr) {
    file_name <- paste0(path_in, "ibd_sharing_carrier_BeauceVSUrbanQc_", v, "_chr", c, ".txt.sharing.by.pos")
    if (file.exists(file_name)) {
      if (file.info(file_name)$size == 0) {
        warning(paste("File is empty:", file_name))
        next
      }
      df <- read.table(file_name)
      df$chrpos <- paste0(df[, 1], ":", df[, 2])
      names(df) <- c("chrom", "pos", "nind", "chrpos")
      df$chrom <- NULL
      chr_list[[c]] <- df
    } else {
      warning(paste("File not found:", file_name))
      next
    }
  }
  results_list[[v]] <- chr_list
}

### Combine chromosome-specific data into one table per variant
combined_results <- list()
for (v in names(results_list)) {
  print(v)
  if (length(results_list[[v]]) > 0) {
    df <- do.call(rbind, results_list[[v]])
    df$v <- v
    combined_results[[v]] <- df
  } else {
    warning(paste("No valid data frames for", v))
  }
}

### Chromosome sizes (GRCh38) used for cumulative positioning
chromosome_lengths <- data.frame(
  chr = seq(1, 22),
  length = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
             159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
             114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
             58617616, 64444167, 46709983, 50818468)
)
chromosome_lengths$cum_length <- cumsum(chromosome_lengths$length)

### Load BIM file containing SNP positions
bim <- read.table(bim_file)
bim$V7 <- paste0(bim$V1, ":", bim$V4)

### Add chromosome and variant information to each variant result
modified_data_list <- list()
for (v in names(combined_results)) {
  df <- combined_results[[v]]
  df <- merge(df, bim, by.x = "chrpos", by.y = "V7", all.x = TRUE, all.y = TRUE)
  names(df) <- c("chrpos", "pos", "nind", "variant", "chr", "rsid", "cm", "posbim", "A1", "A2")
  df$variant <- gsub(":", "_", df$variant)
  df$var_pos <- sub(".*_", "", df$variant)
  df$var_chr <- sub("_.*", "", df$variant)

  df <- merge(df, chromosome_lengths, by = "chr")
  df$position <- df$pos + df$cum_length

  df$var_chr <- as.numeric(df$var_chr)
  df <- merge(df, chromosome_lengths, by.x = "var_chr", by.y = "chr")
  df$var_pos <- as.numeric(df$var_pos) + as.numeric(df$'cum_length.y')

  modified_data_list[[v]] <- df
}

### Load IBD sharing input files (used to count carriers)
ind_files <- list.files(path_in, pattern = "ibd_sharing_carrier_BeauceVSUrbanQc_.*_chr.*\\.txt$", full.names = TRUE)
var <- unique(gsub("ibd_sharing_carrier_BeauceVSUrbanQc_(.*)_chr.*\\.txt", "\\1", basename(ind_files)))
ch <- unique(gsub("ibd_sharing_carrier_BeauceVSUrbanQc_.*_chr(.*).txt", "\\1", basename(ind_files)))

ind_list <- list()
for (v in var) {
  chr_list <- list()
  for (c in ch) {
    file_name <- paste0(path_in, "ibd_sharing_carrier_BeauceVSUrbanQc_", v, "_chr", c, ".txt")
    if (file.exists(file_name)) {
      if (file.info(file_name)$size == 0) {
        warning(paste("File is empty:", file_name))
        next
      }
      df <- read.table(file_name, header = TRUE)
      chr_list[[c]] <- df
    } else {
      warning(paste("File not found:", file_name))
      next
    }
  }
  ind_list[[v]] <- chr_list
}


# Combine all individuals sharing information per variant
list_to_df <- function(lst) {
  do.call(rbind, lst)
}
df_list <- lapply(ind_list, function(inner_list) {
  as.data.frame(list_to_df(inner_list))
})

### Count unique carriers per variant
carrier_list <- list()
for (v in var) {
  df <- as.data.frame(df_list[[v]])
  carriers <- unique(c(df$ind1, df$ind2))
  carrier_list[[v]] <- length(carriers)
}

### Compute sharing proportion for each variant
proportion_df <- list()
for (v in names(modified_data_list)) {
  df <- modified_data_list[[v]]
  carriers <- carrier_list[[v]]
  pairs <- (carriers * (carriers - 1)) / 2
  df$proportion <- df$nind / pairs
  df$ncarriers <- carriers
  proportion_df[[v]] <- df
}

# Print maximum sharing proportion observed
max_values <- sapply(proportion_df, function(df) max(df$proportion, na.rm = TRUE))
print("max proportion ")
print(max_values)

### Generate and save plots of IBD sharing proportions
graph_list <- list()
for (i in seq_along(proportion_df)) {
  df <- proportion_df[[i]]
  xaxis <- c()
  vlines <- c()

  for (ichr in 1:22) {
    mymin <- min(df$position[df$chr == ichr])
    mymax <- max(df$position[df$chr == ichr])
    mypos <- ((mymax - mymin) / 2) + mymin
    xaxis <- c(xaxis, mypos)
    vlines <- c(vlines, mymax)
  }
  graph_list[[i]] <- df
  gg <- ggplot(df, aes(x = position, y = proportion)) +
    geom_line(linewidth = 0.5) +
    ylim(c(0, 1)) +
    labs(title = paste0("Proportion of pairs sharing IBD segment between \n the ", df$ncarriers, " carriers of ", unique(df$variant))) +
    geom_vline(xintercept = unique(df$var_pos), linetype = "longdash", color = "red", linewidth = 0.25, alpha = 0.70) +
    geom_vline(xintercept = vlines, linetype = "dashed", color = "darkgrey", linewidth = 0.25) +
    scale_x_continuous(name = 'Chromosome', breaks = as.numeric(xaxis), labels = seq(1, 22)) +
    labs(x = "Chromosome", y = "Proportion of pairs sharing") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  var <- unique(df$variant)
  plot_title <- paste0(path_in, "IDB_sharing_proportion_variants_BeauceVSUrbanQc_", var, ".png")
  ggsave(plot_title, plot = gg, width = 12, height = 8)
}

# Save all data used for plotting
save(graph_list, file = paste0(path_in, "IDB_sharing_proportion_variants_BeauceVSUrbanQc.RData"))



