suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(svglite)
  library(patchwork)
})


obj <- readRDS("../objects/tcell_TCRab_annotated_final04_clean.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)

# ==============================================================================
# GENE GROUPS
# ==============================================================================

gene_groups <- list(
  "B-cell chemoattraction / core Tph" = c("Cxcl13", "Icos", "Il21", "Cd40lg", "Pdcd1"),
  "Exhaustion-associated TFs"          = c("Tox", "Batf", "Prdm1", "Bcl6"),
  "Tfh-exclusion markers"              = c("Cxcr5", "Cxcr6", "Sell", "Ccr7"),
  "Costimulation / survival"           = c("Tnfrsf4", "Tnfsf8", "Cd200")
)

# Filter to genes present in the object
all_genes_requested <- unlist(gene_groups)
missing <- setdiff(all_genes_requested, rownames(obj))
if (length(missing) > 0)
  message("Genes not found in object (skipped): ", paste(missing, collapse = ", "))

gene_groups <- lapply(gene_groups, function(g) intersect(g, rownames(obj)))
gene_groups  <- gene_groups[sapply(gene_groups, length) > 0]

# ==============================================================================
# BINARY GROUPING
# ==============================================================================

obj$tph_group <- ifelse(obj$celltype_merged == "Tfh/Tph-like", "Tfh/Tph-like", "Other T")
obj$tph_group[is.na(obj$celltype_merged)] <- NA

# ==============================================================================
# EXPRESSION — long-format data frame
# ==============================================================================

all_genes <- unlist(gene_groups)
expr_mat   <- GetAssayData(obj, layer = "data")[all_genes, , drop = FALSE]

plot_df <- as.data.frame(t(as.matrix(expr_mat))) %>%
  tibble::rownames_to_column("cell") %>%
  pivot_longer(-cell, names_to = "gene", values_to = "expression") %>%
  left_join(
    data.frame(cell = colnames(obj), tph_group = obj$tph_group,
               stringsAsFactors = FALSE),
    by = "cell"
  ) %>%
  filter(!is.na(tph_group))

# Factor levels: Tfh/Tph-like first so it appears on the left
plot_df$tph_group <- factor(plot_df$tph_group, levels = c("Tfh/Tph-like", "Other T"))

group_colors <- c("Tfh/Tph-like" = "#B2182B", "Other T" = "#AAAAAA")

# ==============================================================================
# HELPER: one row of violins for a set of genes
# ==============================================================================

make_group_violin <- function(genes, group_label) {
  df <- plot_df %>% filter(gene %in% genes) %>%
    mutate(gene = factor(gene, levels = genes))  # preserve biological order

  ggplot(df, aes(x = tph_group, y = expression, fill = tph_group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.8, linewidth = 0.3) +
    geom_jitter(size = 0.15, alpha = 0.2, width = 0.18, color = "grey20") +
    stat_summary(fun = median, geom = "point", size = 2.5, color = "black",
                 position = position_dodge(0)) +
    facet_wrap(~ gene, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = group_colors) +
    labs(title = group_label, x = NULL, y = "Log-norm expression") +
    theme_classic(base_size = 11) +
    theme(
      panel.grid        = element_blank(),
      axis.text.x       = element_text(angle = 35, hjust = 1, size = 9),
      axis.text.y       = element_text(size = 8),
      axis.title.y      = element_text(size = 9),
      strip.text        = element_text(face = "italic", size = 10),
      strip.background  = element_rect(fill = "#EEF2F7", color = NA),
      plot.title        = element_text(face = "bold", size = 11, hjust = 0),
      legend.position   = "none"
    )
}

# ==============================================================================
# BUILD PANELS
# ==============================================================================

panels <- mapply(make_group_violin,
                 genes       = gene_groups,
                 group_label = names(gene_groups),
                 SIMPLIFY    = FALSE)

# ==============================================================================
# COMBINE WITH PATCHWORK
# ==============================================================================

final_plot <- Reduce(`/`, panels) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ==============================================================================
# SAVE
# ==============================================================================

out <- "../figures/Tph_marker_violin_combined.svg"
svglite(out, width = 20, height = 16)
print(final_plot)
dev.off()
cat("Saved:", out, "\n")
