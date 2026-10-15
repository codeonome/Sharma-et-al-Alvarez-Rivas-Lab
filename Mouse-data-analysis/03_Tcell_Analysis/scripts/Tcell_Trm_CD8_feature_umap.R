suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(svglite)
  library(patchwork)
})

set.seed(123)

# ==============================================================================
# FEATURE PLOTS — CD8/CD4 module scores + raw lineage category for cluster 6,
# on (A) the main T-cell UMAP and (B) cluster 6's own local re-embedding.
# Purpose: see WHERE the double-negative ("dropout") cells sit relative to the
# Cd8+ and Cd4+ cells — do they colocalize with CD8 (=> likely CD8 dropout) or
# form their own territory?
# ==============================================================================

tcr_pat <- "(^Trav|^Trbv|^Trgv|^Trdv|^Trac$|^Trbc1$|^Trbc2$|^Traj|^Trbj)"
obj <- readRDS("../objects/tcell_TCRab_annotated_final04_clean.rds")
DefaultAssay(obj) <- "RNA"; obj <- JoinLayers(obj)

# module scores (same as identity check)
obj <- AddModuleScore(obj, features = list(intersect(c("Cd8a","Cd8b1"), rownames(obj))),
                      name = "CD8score", seed = 123)
obj <- AddModuleScore(obj, features = list(intersect(c("Cd4","Cd40lg"), rownames(obj))),
                      name = "CD4score", seed = 123)

trm_cells <- colnames(obj)[which(obj$celltype_merged == "Trm type mixed CD8+ CD4")]

# raw lineage category (the >0 4-way used for the "dropout" discussion)
raw <- FetchData(obj, vars = c("Cd8a","Cd8b1","Cd4"), cells = trm_cells, layer = "counts")
lineage <- ifelse((raw$Cd8a>0 | raw$Cd8b1>0) & raw$Cd4==0, "CD8+ / CD4-",
           ifelse(raw$Cd4>0 & raw$Cd8a==0 & raw$Cd8b1==0, "CD4+ / CD8-",
           ifelse((raw$Cd8a>0 | raw$Cd8b1>0) & raw$Cd4>0, "double-positive",
                  "double-negative (dropout)")))
names(lineage) <- trm_cells
obj$c6_lineage <- NA_character_
obj$c6_lineage[trm_cells] <- lineage

pal_lin <- c("CD8+ / CD4-"="#4E79A7","CD4+ / CD8-"="#E15759",
             "double-positive"="#B07AA1","double-negative (dropout)"="grey55")
grad <- c("grey90","#FEE5D9","#FC9272","#DE2D26","#67000D")

base_theme <- theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(), axis.line = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5))

# ------------------------------------------------------------------------------
# A) MAIN T-CELL UMAP (context)
# ------------------------------------------------------------------------------
emb <- Embeddings(obj, "umap.pca_TCRab") %>% as.data.frame() %>% rownames_to_column("cell")
xg <- colnames(emb)[2]; yg <- colnames(emb)[3]
mdf <- emb %>%
  left_join(obj@meta.data %>% rownames_to_column("cell") %>%
              select(cell, CD8score1, CD4score1, c6_lineage), by="cell") %>%
  mutate(UMAP2 = -.data[[yg]])

feat_main <- function(col, ttl) {
  ggplot(mdf[order(mdf[[col]]),], aes(.data[[xg]], UMAP2, color=.data[[col]])) +
    geom_point(size=0.35, alpha=0.8) +
    scale_color_gradientn(colors=grad, name="score") +
    base_theme + labs(title=ttl, x="UMAP_1", y="UMAP_2")
}
p_c6 <- ggplot(mdf, aes(.data[[xg]], UMAP2)) +
  geom_point(data=filter(mdf, is.na(c6_lineage)), color="grey88", size=0.3, alpha=0.6) +
  geom_point(data=filter(mdf, !is.na(c6_lineage)), aes(color=c6_lineage), size=0.7) +
  scale_color_manual(values=pal_lin, name="Cluster 6") +
  base_theme + labs(title="Cluster 6 cells on main UMAP (by lineage)", x="UMAP_1", y="UMAP_2")

main_fig <- p_c6 + feat_main("CD8score1","CD8 module score (all cells)") +
  feat_main("CD4score1","CD4 module score (all cells)") + plot_layout(nrow=1)
svglite("../figures/Trm_CD8_featureplot_mainUMAP.svg", width = 16, height = 4.8)
print(main_fig); dev.off()
cat("Saved: Trm_CD8_featureplot_mainUMAP.svg\n")

# ------------------------------------------------------------------------------
# B) CLUSTER-6 LOCAL RE-EMBEDDING
# ------------------------------------------------------------------------------
cl6 <- subset(obj, cells = trm_cells)
cl6 <- JoinLayers(cl6)
cl6 <- NormalizeData(cl6, verbose=FALSE)
cl6 <- FindVariableFeatures(cl6, nfeatures=2000, verbose=FALSE)
VariableFeatures(cl6) <- setdiff(VariableFeatures(cl6),
                                 grep(tcr_pat, VariableFeatures(cl6), value=TRUE, ignore.case=TRUE))
cl6 <- ScaleData(cl6, vars.to.regress=c("percent.mt","nCount_RNA"), verbose=FALSE)
cl6 <- RunPCA(cl6, npcs=30, seed.use=123, verbose=FALSE)
cl6 <- RunUMAP(cl6, reduction="pca", dims=1:10, reduction.name="umap.cl6", seed.use=123, verbose=FALSE)

lemb <- Embeddings(cl6, "umap.cl6") %>% as.data.frame() %>% rownames_to_column("cell")
lx <- colnames(lemb)[2]; ly <- colnames(lemb)[3]
ldf <- lemb %>%
  left_join(cl6@meta.data %>% rownames_to_column("cell") %>%
              select(cell, CD8score1, CD4score1, c6_lineage), by="cell") %>%
  left_join(FetchData(cl6, vars=c("Cd8b1","Cd4"), layer="data") %>% rownames_to_column("cell"), by="cell")

p_lin <- ggplot(ldf, aes(.data[[lx]], .data[[ly]], color=c6_lineage)) +
  geom_point(size=1.1) + scale_color_manual(values=pal_lin, name=NULL) +
  base_theme + labs(title="Lineage category (raw >0)", x="UMAP1", y="UMAP2")
feat_loc <- function(col, ttl) ggplot(ldf[order(ldf[[col]]),], aes(.data[[lx]], .data[[ly]], color=.data[[col]])) +
  geom_point(size=1.1) + scale_color_gradientn(colors=grad, name="") +
  base_theme + labs(title=ttl, x="UMAP1", y="UMAP2")

local_fig <- p_lin + feat_loc("CD8score1","CD8 module score") + feat_loc("CD4score1","CD4 module score") +
  feat_loc("Cd8b1","Cd8b1 (norm)") + feat_loc("Cd4","Cd4 (norm)") + plot_layout(nrow=2)
svglite("../figures/Trm_CD8_featureplot_localUMAP.svg", width = 13, height = 9)
print(local_fig); dev.off()
cat("Saved: Trm_CD8_featureplot_localUMAP.svg\n")

# quick quantitative aid: where do dropout cells sit? mean scores by lineage category
cat("\n=== Mean module scores by cluster-6 lineage category ===\n")
print(ldf %>% group_by(c6_lineage) %>%
        summarise(n=n(), mean_CD8score=round(mean(CD8score1),3),
                  mean_CD4score=round(mean(CD4score1),3), .groups="drop"))
cat("\n=== Cluster complete ===\n")
