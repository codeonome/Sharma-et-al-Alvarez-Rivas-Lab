suppressPackageStartupMessages(library(Seurat))


t_cells_TCRab <- readRDS("../objects/tcell_TCRab_annotated_final04.rds")
cat("Loaded:", ncol(t_cells_TCRab), "cells\n")

mait_cells <- colnames(t_cells_TCRab)[is.na(t_cells_TCRab$celltype_merged)]
cat("MAIT/innate-like cells to remove (celltype_merged == NA):", length(mait_cells), "\n")
print(table(t_cells_TCRab$condition[mait_cells]))

t_cells_TCRab_clean <- subset(t_cells_TCRab, cells = setdiff(colnames(t_cells_TCRab), mait_cells))
cat("\nRemaining cells:", ncol(t_cells_TCRab_clean), "\n")
cat("celltype_merged breakdown after removal:\n")
print(table(t_cells_TCRab_clean$celltype_merged, useNA = "always"))

saveRDS(t_cells_TCRab_clean, "../objects/tcell_TCRab_annotated_final04_clean.rds")
cat("\nSaved: ../objects/tcell_TCRab_annotated_final04_clean.rds (original tcell_TCRab_annotated_final04.rds left untouched)\n")
