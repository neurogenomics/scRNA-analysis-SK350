# Run view-soupercell.R

sce <- singlets_only

sce[["RNA"]]$counts
sce[["RNA"]]$counts |> ncol()

sex.classified <- speckle::classifySex(
    sce[["RNA"]]$counts,
    genome = "Hs"   
)

sex.classified
sex.classified |> nrow()


# Combine with Seurat metadata
sce$sex <- sex.classified
sex.table <- table(sce$sex, sce$donor_cluster) %>% as.data.frame.matrix()


# Write output to file
out_dir <- "out"
file <- file.path(out_dir, paste0(sample_id, "_sex_classification.txt"))
write.table(sex.table, file = file, quote = FALSE)







