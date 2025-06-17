library(tidyverse)
library(Seurat)

sample_id <- "SK350_Pool1_NBB"
# sample_id <- "SK350_Pool2_PDUK"

seurat_obj <- paste0(
  "cellranger-output/count/",
  sample_id,
  "/filtered_feature_bc_matrix"
) |>
  Read10X() |>
  CreateSeuratObject(counts = _)

clusters <- paste0("souporcell-output/", sample_id, "/clusters.tsv") |>
  read.table(header = TRUE, sep = "\t")

# Rename columns to match what Seurat expects
clusters$donor_cluster <- clusters$assignment
clusters$doublet <- clusters$status == "doublet" # TRUE/FALSE

# Set barcode as rownames
rownames(clusters) <- clusters$barcode

# Subset to barcodes present in Seurat
matched_clusters <- clusters[colnames(seurat_obj), ]

# Add to Seurat metadata
seurat_obj$donor_cluster <- matched_clusters$donor_cluster
seurat_obj$doublet <- matched_clusters$doublet

# Normalize (based on UMI counts) and scale (based on gene expression) data
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)

seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

seurat_obj <- RunTSNE(seurat_obj)

# Visualize
singlets_only <- seurat_obj[, nchar(seurat_obj$donor_cluster) == 1] |>
  RunTSNE()

dimplt_all <- DimPlot(seurat_obj,
  group.by = "donor_cluster",
  label = TRUE,
  reduction = "tsne"
) +
  labs(
    title = paste("tSNE plot for", sample_id, ("donor clusters")),
    subtitle = paste(
      "Total cells:",
      ncol(seurat_obj),
      "\nTotal singlets:",
      ncol(singlets_only)
    ),
  )
dimplt_all
ggplot2::ggsave(
  dimplt_all,
  filename = paste0("out/", sample_id, "_souporcell_tsne.png"),
  width = 10,
  height = 10
)

DimPlot(
  singlets_only,
  group.by = "donor_cluster",
  label = TRUE,
  reduction = "tsne"
)



seurat_obj$expression.IRF8 <- FetchData(seurat_obj, vars = "IRF8") > 0
dimplt_irf8 <- DimPlot(seurat_obj,
  group.by = "expression.IRF8",
  label = TRUE,
  reduction = "tsne"
) +
  labs(
    title = paste("tSNE plot for", sample_id, ("IRF8 expression")),
    subtitle = paste(
      "Total cells:",
      ncol(seurat_obj),
      "\nIRF8 positive cells:",
      sum(seurat_obj$expression.IRF8)
    ),
  )
dimplt_irf8
ggplot2::ggsave(
  dimplt_irf8,
  filename = paste0("out/", sample_id, "_souporcell_IRF8_expression.png"),
  width = 10,
  height = 10
)


# Table - IRF8 cell expression count per donor cluster
irf8_counts <- table(seurat_obj$donor_cluster, seurat_obj$expression.IRF8) |>
  as.data.frame.matrix() |>
  rownames_to_column("donor_cluster") |>
  # pivot_wider(names_from = expression.IRF8, values_from = n) |>
  rename("IRF8_negative" = `FALSE`, "IRF8_positive" = `TRUE`) |>
  mutate(
    "total_cells" = IRF8_negative + IRF8_positive,
    "fraction_IRF8_positive" = IRF8_positive / total_cells
  ) |>
  relocate("total_cells", .after = "donor_cluster")
write_tsv(irf8_counts,
  file = paste0("out/", sample_id, "_souporcell_IRF8_expression_counts.tsv")
)
cat(
  "IRF8 expression counts saved to out/",
  sample_id,
  "souporcell_IRF8_expression_counts.tsv\n", 
  sep = ""
)

# Visualize - IRF8 expression by donor cluster, bar plot
# Only for singlets
irf8_counts_plt <- irf8_counts[!grepl("/", irf8_counts$donor_cluster), ] |>
  pivot_longer(
    cols = c("IRF8_positive", "IRF8_negative"),
    names_to = "IRF8_expression",
    values_to = "count"
  ) |>
  mutate(IRF8_expression = recode(IRF8_expression,
    "IRF8_negative" = "No",
    "IRF8_positive" = "Yes"
  )) |>
  ggplot(aes(x = donor_cluster, y = count, fill = IRF8_expression)) +
  geom_bar(stat = "identity") +
  labs(
    title = "IRF8 Expression by Donor Cluster",
    subtitle = paste("Pool:", sample_id),
    x = "Donor Cluster",
    y = "Number of Cells",
    fill = "IRF8 Expression?"
  )

irf8_counts_plt
ggplot2::ggsave(
  irf8_counts_plt,
  filename = paste0("out/", sample_id, "_souporcell_IRF8_expression_counts.png"),
  width = 10,
  height = 10
)
