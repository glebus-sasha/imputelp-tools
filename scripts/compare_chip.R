#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(vcfR)
  library(readxl)
  library(ggplot2)
  library(plotly)
})

# ---------------------------
# DEBUG / PROD параметры
# ---------------------------
DEBUG <- FALSE

if (DEBUG) {
  args <- list(
    vcf_dir = "raw/imputed",
    passport = "raw/20250331_kira_monogen_passports_final_table.xlsx",
    descr = "raw/Описание_образцов_RnD_low_pass_PCRfree_SGGM_18032025.xlsx",
    interpretation = "raw/SNV_53_interpretation.xlsx",
    output_prefix = "results/BCF_vs_Chip_March2026"
  )
  message("Running in DEBUG mode with hardcoded arguments")
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
  if(length(cmd_args) < 5) stop("Provide 5 arguments: vcf_dir passport descr interpretation output_prefix")
  args <- list(
    vcf_dir = cmd_args[1],
    passport = cmd_args[2],
    descr = cmd_args[3],
    interpretation = cmd_args[4],
    output_prefix = cmd_args[5]
  )
}

prefix <- args$output_prefix

# ---------------------------
# Загрузка данных
# ---------------------------
passport_df <- read_excel(args$passport) %>% 
  select(`Probe Code`, ABCG2:TTF1)

descr_df <- read_excel(args$descr) %>% 
  filter(Источник == 'Мираторг') %>% 
  select(
    sid = `Название образца`,
    `Probe Code` = Наименование
  )

interpretation_df <- read_excel(args$interpretation, sheet = 'Chip') %>%
  mutate(
    # Преобразуем Chr: NA остаётся NA, иначе добавляем префикс "chr" и убираем ".0" если есть
    Chr = ifelse(
      is.na(Chr) | Chr == "NA",
      NA,
      paste0("chr", str_remove(Chr, "\\.0$"))
    )
  ) %>% 
  select(
    Short_name = Short_name_chip,
    Symptom:OMIA_ID,
    Chr,
    START:END
  )

# ---------------------------
# Функция обработки VCF
# ---------------------------
process_vcf <- function(file_path) {
  file_name <- basename(file_path)
  sample_name <- tools::file_path_sans_ext(file_name)  # убираем .vcf, используем как sid
  
  vcf <- read.vcfR(file_path, verbose = FALSE)
  vcf_df <- vcfR2tidy(vcf, single_frame = TRUE)$dat
  
  vcf_df %>% mutate(sid = sample_name)
}

# ---------------------------
# Основной пайплайн
# ---------------------------
vcf_files <- list.files(path = args$vcf_dir, pattern = "\\.vcf$", full.names = TRUE)

passport_long <- passport_df %>%
  pivot_longer(cols = -`Probe Code`,
               names_to = "Short_name",
               values_to = "status_chip") %>%
  right_join(descr_df) %>% 
  rename(Probe_Code = `Probe Code`)

combined_df <- map_dfr(vcf_files, process_vcf) %>% 
  left_join(descr_df) %>%
  mutate(
    status = case_when(
      gt_GT == "0|0" ~ "FREE",
      gt_GT %in% c("0|1", "1|0") ~ "CARRIER",
      gt_GT == "1|1" ~ "AFFECTED",
      TRUE ~ NA_character_
    )
  ) %>% 
  select(sid, Probe_Code = `Probe Code`, CHROM, POS, status) %>%
  left_join(interpretation_df, by = c("CHROM" = "Chr", "POS" = "END")) %>%
  left_join(passport_long, by = c("Probe_Code", "Short_name", "sid"))

# ---------------------------
# Wide table
# ---------------------------
wide_df <- combined_df %>%
  right_join(passport_long, by = c("Probe_Code", "Short_name", "sid")) %>%
  mutate(
    # coalesce берёт status из VCF, если есть, иначе оставляет NA
    combined_status = paste(status, status_chip.y, sep = "/")
  ) %>%
  select(sid, Probe_Code, Short_name, combined_status) %>%
  pivot_wider(
    id_cols = c(sid, Probe_Code),
    names_from = Short_name,
    values_from = combined_status
  ) %>%
  rename(`Probe Code` = Probe_Code)

write_csv(wide_df, paste0(prefix, "_wide_table.csv"))

# ---------------------------
# Статистика
# ---------------------------
stat_df <- combined_df %>%
  filter(!is.na(status) & !is.na(status_chip)) %>%
  group_by(status, status_chip) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

total_compared <- sum(stat_df$count)
total_matches <- sum(stat_df$count[stat_df$status == stat_df$status_chip])
percent_matches <- round(total_matches / total_compared * 100, 2)

summary_stats <- tibble(
  total_variants = total_compared,
  total_matches = total_matches,
  percent_matches = percent_matches
)

write_csv(stat_df, paste0(prefix, "_combined_status_counts.csv"))
write_csv(summary_stats, paste0(prefix, "_summary_stats.csv"))

# ---------------------------
# Heatmap (ggplot PNG)
# ---------------------------
heatmap_df <- passport_long %>%
  distinct(sid, Probe_Code, Short_name, status_chip) %>%  # возьмём статус chip сразу
  left_join(
    combined_df %>% 
      select(sid, Probe_Code, Short_name, status),
    by = c("sid", "Probe_Code", "Short_name")
  ) %>%
  mutate(
    # Обрабатываем NA
    status_show = ifelse(is.na(status), "Missing in VCF", status),
    status_chip_show = ifelse(is.na(status_chip), "Missing in chip", status_chip),
    fill_state = case_when(
      is.na(status) & is.na(status_chip) ~ "Both missing",
      is.na(status) ~ "Only chip",
      is.na(status_chip) ~ "Only VCF",
      status == status_chip ~ "Match",
      TRUE ~ "Mismatch"
    ),
    hover_text = paste0(
      "Sample: ", sid,
      "<br>SNP: ", Short_name,
      "<br>Status VCF: ", status_show,
      "<br>Status chip: ", status_chip_show,
      "<br>Result: ", fill_state
    )
  ) %>% 
  filter(status_chip != "Both missing")

p <- ggplot(heatmap_df, aes(
  x = Short_name,
  y = sid,
  fill = fill_state,
  text = hover_text
)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c(
      "Match" = "green",
      "Mismatch" = "red",
      "Only VCF" = "orange",
      "Only chip" = "dodgerblue",
      "Both missing" = "grey80"
    ),
    drop = FALSE
  ) +
  theme_minimal() +
  labs(fill = "Comparison", x = "SNP / Chip", y = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(paste0(prefix, "_heatmap.png"), plot = p, width = 12, height = 8)

# ---------------------------
# HTML график (интерактивный)
# ---------------------------
# интерактив
interactive_plot <- ggplotly(p, tooltip = "text")

htmlwidgets::saveWidget(
  interactive_plot,
  paste0(prefix, "_heatmap.html"),
  selfcontained = TRUE
)

cat("Done! All results saved with prefix:", prefix, "\n")