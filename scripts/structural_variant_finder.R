# structural_variant_finder.R

library(readr)
library(dplyr)
library(stringr)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(httr)
library(jsonlite)
library(ggplot2)

#-------------------------------
# Argument Parsing and Setup
#-------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Two input files required: <tumor_file> <normal_file>")

tumor_file <- args[1]
normal_file <- args[2]

if (!file.exists(tumor_file)) stop("Tumor file not found: ", tumor_file)
if (!file.exists(normal_file)) stop("Normal file not found: ", normal_file)

outdir <- dirname(tumor_file)
sample <- tools::file_path_sans_ext(basename(tumor_file))
log_file <- file.path(outdir, paste0(sample, "_pipeline_log.txt"))

log_message <- function(msg) {
  cat(msg, "\n"); write(msg, file = log_file, append = TRUE)
}

log_message(paste0("Pipeline started for sample: ", sample))

#-------------------------------
# Load and Annotate SVs
#-------------------------------
annot_tumor <- read_tsv(tumor_file, comment = "#", show_col_types = FALSE)
annot_normal <- read_tsv(normal_file, comment = "#", show_col_types = FALSE)
if (nrow(annot_tumor) == 0 || nrow(annot_normal) == 0) stop("Empty or invalid input.")

oncogene_list <- readLines("data/oncogenes.txt")
tumor_suppressor_list <- readLines("data/tumor_suppressors.txt")
cancer_genes <- unique(c(oncogene_list, tumor_suppressor_list))

annot_filtered <- annot_tumor %>%
  anti_join(annot_normal, by = c("Chr", "Start", "End", "Gene_name")) %>%
  filter(!is.na(Gene_name)) %>%
  mutate(Functional_Impact = case_when(
    Gene_name %in% oncogene_list ~ "Oncogene",
    Gene_name %in% tumor_suppressor_list ~ "Tumor Suppressor",
    Gene_name %in% cancer_genes ~ "Cancer Gene",
    TRUE ~ "Uncertain"
  ))

#-------------------------------
# Clinical API + COSMIC Local
#-------------------------------
civic_api_key <- "YOUR_CIVIC_API_KEY"
oncokb_api_key <- "YOUR_ONCOKB_API_KEY"
cosmic_local <- read_tsv("data/cosmic_gene_hits.tsv")

query_api <- function(gene, source) {
  if (source == "oncokb") {
    url <- paste0("https://oncokb.org/api/v1/genes/", gene, "/variants")
    res <- GET(url, add_headers(Authorization = paste("Bearer", oncokb_api_key)))
  } else {
    url <- paste0("https://civicdb.org/api/variants?entrez_id=", gene)
    res <- GET(url, add_headers(Authorization = paste("Bearer", civic_api_key)))
  }
  if (status_code(res) == 200) fromJSON(content(res, as = "text")) else NULL
}

query_with_fallback <- function(gene) {
  civic <- query_api(gene, "civic")
  oncokb <- query_api(gene, "oncokb")
  cosmic <- cosmic_local %>% filter(Gene == gene)
  score <- 0
  if (!is.null(civic) && length(civic$records) > 0) score <- score + 3
  if (!is.null(oncokb) && length(oncokb) > 0) score <- score + 3
  if (nrow(cosmic) > 0) score <- score + 2
  return(score)
}

annot_filtered <- annot_filtered %>%
  rowwise() %>%
  mutate(Clinical_Evidence_Score = query_with_fallback(Gene_name)) %>%
  ungroup()

#-------------------------------
# Final Score and Output
#-------------------------------
annot_filtered <- annot_filtered %>%
  mutate(Combined_Score = Clinical_Evidence_Score +
           ifelse(Functional_Impact %in% c("Oncogene", "Tumor Suppressor"), 5, 2),
         Clinical_Significance = case_when(
           Combined_Score >= 8 ~ "Pathogenic",
           Combined_Score >= 5 ~ "Likely Pathogenic",
           Combined_Score >= 3 ~ "VUS",
           TRUE ~ "Benign"
         ))

write.xlsx(annot_filtered, file.path(outdir, paste0(sample, "_clinical_significance.xlsx")))
log_message("Pipeline completed with scoring and annotation.")
