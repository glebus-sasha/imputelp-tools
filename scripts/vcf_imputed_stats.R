#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(dplyr)
  library(ggplot2)
  library(optparse)
  library(readr)
})

# ---- DEBUG / DEV PARAMETERS (uncomment for testing) ----
# if (FALSE) {
#   opt <- list(
#     input = "raw/SSP_1378_BB254867491.vcf.gz",
#     output = "vcf_stats_out"
#   )
# }

# ---- define CLI options ----
if (!exists("opt")) {
  option_list <- list(
    make_option("--input", type="character", help="Input VCF file path (.vcf.gz)"),
    make_option("--output", type="character", default="vcf_stats_out", help="Output directory")
  )
  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
}

# ---- log inputs ----
cat("Input parameters:\n")
print(opt)

# ---- check required arguments ----
if (!file.exists(opt$input)) stop("Input VCF file does not exist: ", opt$input)
if (!dir.exists(opt$output)) dir.create(opt$output, recursive = TRUE)

# ---- read VCF ----
cat("Reading VCF...\n")
vcf <- readVcf(opt$input, "genome")  # genome can be replaced with specific assembly

# ---- extract INFO and GP ----
info_vector <- as.numeric(info(vcf)$INFO)

# GP: матрица вариантов x образцы x 3
gp_array <- geno(vcf)$GP
num_samples <- dim(gp_array)[2]

# max GP по вариантам (по всем образцам)
if (num_samples == 1) {
  gp_matrix <- gp_array[,1,]  # матрица вариантов x 3
  max_gp <- apply(gp_matrix, 1, max)
} else {
  max_gp <- apply(gp_array, 1:1, function(x) max(x, na.rm=TRUE))  # для нескольких образцов
}

# ---- summary table ----
summary_df <- data.frame(
  CHROM = as.character(seqnames(rowRanges(vcf))),
  POS   = start(rowRanges(vcf)),
  INFO  = info_vector,
  max_GP = max_gp
)

# ---- summary statistics ----
summary_stats <- summary_df %>%
  summarise(
    total_variants = n(),
    mean_INFO      = mean(INFO, na.rm=TRUE),
    mean_max_GP    = mean(max_GP, na.rm=TRUE)
  )

cat("Summary statistics:\n")
print(summary_stats)

# ---- save summary statistics (metrics only) ----
metrics_csv <- file.path(opt$output, "vcf_metrics_summary.csv")
write_csv(summary_stats, metrics_csv)
cat("Summary metrics saved to:", metrics_csv, "\n")

# ---- save detailed table per variant (optional) ----
variant_csv <- file.path(opt$output, "vcf_metrics_per_variant.csv")
write_csv(summary_df, variant_csv)
cat("Detailed variant table saved to:", variant_csv, "\n")

# ---- plots ----
p1 <- ggplot(summary_df, aes(x=INFO)) +
  geom_histogram(bins=50, fill="steelblue", color="black") +
  theme_minimal() +
  labs(title="Distribution of INFO", x="INFO", y="Count")
ggsave(filename=file.path(opt$output, "INFO_histogram.png"), plot=p1, width=6, height=4)

p2 <- ggplot(summary_df, aes(x=max_GP)) +
  geom_histogram(bins=50, fill="darkgreen", color="black") +
  theme_minimal() +
  labs(title="Distribution of max GP", x="max_GP", y="Count")
ggsave(filename=file.path(opt$output, "maxGP_histogram.png"), plot=p2, width=6, height=4)

cat("Plots saved in:", opt$output, "\n")
cat("VCF statistics script finished successfully.\n")