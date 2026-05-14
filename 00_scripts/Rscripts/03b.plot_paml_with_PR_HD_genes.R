#!/usr/bin/env Rscript

# === Parse command-line arguments ===
args <- commandArgs(trailingOnly = TRUE)
use_raw_data <- FALSE
max_ds <- 1.0

i <- 1
while (i <= length(args)) {
  if (args[i] == "--raw") {
    use_raw_data <- TRUE
  } else if (args[i] == "--max_ds" && i + 1 <= length(args)) {
    max_ds <- as.numeric(args[i + 1])
    i <- i + 1
  } else if (args[i] %in% c("-h", "--help")) {
    cat("\nUsage: Rscript 03b.plot_paml_with_PR_HD_genes.R [OPTIONS]\n\n")
    cat("Options:\n")
    cat("  --raw             Plot all dS values without filtering\n")
    cat("  --max_ds VALUE    Max dS threshold for filtering (default: 1.0)\n")
    cat("  -h, --help        Show this help message\n\n")
    cat("Examples:\n")
    cat("  Rscript 03b.plot_paml_with_PR_HD_genes.R\n")
    cat("  Rscript 03b.plot_paml_with_PR_HD_genes.R --max_ds 0.5\n")
    cat("  Rscript 03b.plot_paml_with_PR_HD_genes.R --raw\n\n")
    quit("no")
  }
  i <- i + 1
}

cat("\n========== PR/HD Annotation for dS Plot ==========\n")
cat("Mode:    ", if (use_raw_data) "raw (no filter)" else paste0("filtered (dS < ", max_ds, ")"), "\n")
cat("max_ds:  ", max_ds, "\n\n")

req_pkgs <- c("ggplot2", "dplyr", "magrittr", "patchwork", "wesanderson")
missing_pkgs <- setdiff(req_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  cat("Installing missing packages:", paste(missing_pkgs, collapse=", "), "\n")
  install.packages(missing_pkgs, repos="https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(magrittr)
  library(patchwork)
  library(wesanderson)
})

# === PR/HD Gene Positions (on Mlag129A1 ancestral scaffolds) ===
pr_hd_data <- tribble(
  ~contig, ~gene, ~start, ~end,
  "Mlag129A1_contig11", "PR",  969702, 970887,
  "Mlag129A1_contig8",  "HD1", 286526, 287848,
  "Mlag129A1_contig8",  "HD2", 288269, 288826
)

cat("PR/HD Marker Positions:\n")
print(as.data.frame(pr_hd_data))
cat("\n")

# === File paths ===
paml_cml_file   <- "02_results/paml/results_codeml.txt"
paml_yn_file    <- "02_results/paml/results_YN.txt"
sco_file        <- "02_results/paml/single.copy.orthologs"
anc_bed_file    <- "ancestral_sp/ancestral_sp.bed"
scaffold_file   <- "scaffold.txt"

for (f in c(paml_cml_file, sco_file, anc_bed_file, scaffold_file)) {
  if (!file.exists(f)) stop("Required file not found: ", f)
}

# === Load data ===

# 1. Raw PAML dS values (geneX = haplo1 gene ID)
if (file.exists(paml_cml_file)) {
  paml_raw <- read.table(paml_cml_file) %>%
    set_colnames(c("Ds", "Dn", "dNdS", "geneX", "geneY"))
  cat("Using codeml results\n")
} else {
  paml_raw <- read.table(paml_yn_file) %>%
    set_colnames(c("Ds", "SEDs", "Dn", "SEDn", "geneX", "geneY"))
  cat("Using yn00 results\n")
}

# 2. Single-copy orthologs: ancestral gene ID -> haplo1/haplo2 gene IDs
sco <- read.table(sco_file, sep = "\t") %>%
  set_colnames(c("ortho", "anc_gene", "geneX", "geneY"))

# 3. Ancestral genome bed file: positions for ALL ancestral genes
anc_bed <- read.table(anc_bed_file) %>%
  set_colnames(c("scaff", "start", "end", "anc_gene"))

# 4. Scaffold orientation info
scaf_orient <- read.table(scaffold_file, sep = "\t") %>%
  set_colnames(c("haplo", "scaff", "order"))

# === Build full position table from ancestral bed ===
target_scaffs <- scaf_orient$scaff

anc_target <- anc_bed %>%
  filter(scaff %in% target_scaffs)

cat("Ancestral genes on target scaffolds:", nrow(anc_target), "\n")

# Scaffold display order: decreasing sort matches original 03.plot_paml.R
# (contig8 first, contig11 second for these scaffolds)
scaff_levels <- sort(target_scaffs, decreasing = TRUE)

# Join: ancestral position -> haplo1 gene ID -> dS
# orderchp is assigned to ALL ancestral genes (including those without dS),
# preserving physical spacing in Panel B - matching original 03.plot_paml.R
all <- anc_target %>%
  left_join(sco %>% select(anc_gene, geneX), by = "anc_gene") %>%
  left_join(paml_raw %>% select(geneX, Ds), by = "geneX") %>%
  left_join(scaf_orient %>% select(scaff, order), by = "scaff") %>%
  mutate(scaff = factor(scaff, levels = scaff_levels)) %>%
  arrange(scaff, start) %>%
  group_by(scaff) %>%
  mutate(St = ifelse(order == "N", start, rev(start))) %>%
  arrange(St, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(orderchp = row_number())

# Subset with dS values for statistics and Panel A
all_ds <- all %>% filter(!is.na(Ds))

# Find orderchp positions of PR/HD markers for Panel B annotation
# (nearest ancestral gene by midpoint position; HD1 and HD2 merged as "HD")
annotation_b <- pr_hd_data %>%
  filter(gene %in% c("PR", "HD1")) %>%
  rowwise() %>%
  mutate(
    x_pos = {
      nearby <- all %>% filter(as.character(scaff) == contig)
      nearby$orderchp[which.min(abs(nearby$start - (start + end) / 2))]
    },
    label = ifelse(gene == "HD1", "HD", gene)
  ) %>%
  ungroup() %>%
  select(x_pos, label)

cat("Total ancestral genes with orderchp:", nrow(all), "\n")
cat("Genes with dS values:", nrow(all_ds), "\n")
cat("Scaffolds:", paste(unique(as.character(all_ds$scaff)), collapse = ", "), "\n")
cat("\ndS statistics (before filtering):\n")
cat("  n =", nrow(all_ds), "\n")
cat("  min =", round(min(all_ds$Ds), 4), "\n")
cat("  median =", round(median(all_ds$Ds), 4), "\n")
cat("  mean =", round(mean(all_ds$Ds), 4), "\n")
cat("  max =", round(max(all_ds$Ds), 4), "\n")
cat("  values > 0.5:", sum(all_ds$Ds > 0.5), "\n")
cat("  values > 1.0:", sum(all_ds$Ds > 1.0), "\n\n")

# === Apply filter ===
# Panel A uses all_ds (all genes with dS), y-axis clipped by coord_cartesian
# Panel B uses filtered data but retains original orderchp for correct x-axis spacing
if (use_raw_data) {
  plot_data_b <- all_ds
  ylim_max_a <- max(all_ds$Ds) * 1.1
  ylim_max_b <- ylim_max_a
  cat("No dS filter applied (--raw mode)\n")
} else {
  plot_data_b <- all_ds %>% filter(Ds < max_ds)
  n_removed <- nrow(all_ds) - nrow(plot_data_b)
  cat("Genes removed by filter (dS >=", max_ds, "):", n_removed, "\n")
  ylim_max_a <- max(plot_data_b$Ds) * 1.15
  ylim_max_b <- max(max(plot_data_b$Ds) * 1.1, 0.5)
}

# === Theme ===
th_plot <- theme(
  axis.title.x  = element_text(size = 14, family = "Helvetica", face = "bold"),
  axis.text.x   = element_text(size = 12, family = "Helvetica"),
  axis.title.y  = element_text(size = 16, family = "Helvetica", face = "bold"),
  axis.text.y   = element_text(size = 12, family = "Helvetica"),
  strip.text.x  = element_text(size = 14, family = "Helvetica"),
  panel.grid.major = element_blank(),
  plot.title    = element_text(family = "Helvetica", face = "bold", size = 18, hjust = 0.02)
)

# === Fig A: dS by physical position ===
# Uses all_ds (all genes with dS values); y-axis clipped by coord_cartesian
# HD1 and HD2 are adjacent - show single "HD" label at HD1 midpoint
annotation_labels <- tribble(
  ~scaff,                   ~x_pos,                    ~label,
  "Mlag129A1_contig11",     (969702  + 970887)  / 2,   "PR",
  "Mlag129A1_contig8",      (286526  + 287848)  / 2,   "HD"
) %>%
  filter(scaff %in% levels(all_ds$scaff)) %>%
  mutate(scaff = factor(scaff, levels = levels(all_ds$scaff)))

p_a <- all_ds %>%
  ggplot(aes(x = start, y = Ds)) +
  facet_wrap(~scaff, scale = "free_x") +
  geom_point(size = 1) +
  geom_text(
    data = annotation_labels,
    aes(x = x_pos, y = ylim_max_a * 0.96, label = label),
    size = 4, color = "#4472C4", fontface = "italic", inherit.aes = FALSE
  ) +
  theme_classic() +
  coord_cartesian(ylim = c(0, ylim_max_a)) +
  xlab("position along chromosome") +
  ylab(expression(italic(d[S]))) +
  th_plot +
  theme(legend.position = "none") +
  ggtitle("A")

# === Fig B: dS by gene order along ancestral reference ===
ncolors <- length(unique(as.character(plot_data_b$scaff)))
color_palette <- colorRampPalette(wes_palette(name = "GrandBudapest1"))(ncolors)

p_b <- plot_data_b %>%
  ggplot(aes(x = orderchp, y = Ds, colour = scaff)) +
  geom_point(size = 1) +
  geom_text(
    data = annotation_b,
    aes(x = x_pos, y = ylim_max_b * 0.96, label = label),
    size = 4, color = "#4472C4", fontface = "italic", inherit.aes = FALSE
  ) +
  theme_classic() +
  coord_cartesian(ylim = c(0, ylim_max_b)) +
  xlab("order along reference") +
  ylab(expression(italic(d[S]))) +
  scale_color_manual(values = color_palette) +
  th_plot +
  theme(legend.position = "none") +
  ggtitle("B")

# === Combine and save ===
combined <- p_a / p_b

if (use_raw_data) {
  output_file <- "02_results/dsplots/dS_with_PROHD_raw.pdf"
} else if (max_ds == 1.0) {
  output_file <- "02_results/dsplots/dS_with_PROHD.pdf"
} else {
  output_file <- sprintf("02_results/dsplots/dS_with_PROHD_maxds%.1f.pdf", max_ds)
}

pdf(file = output_file, width = 14, height = 8)
print(combined)
dev.off()

cat("\nOutput saved to:", output_file, "\n")
cat("=========================================\n\n")
