# Load mg.combined using load.R first
library(tidyverse)
library(SingleR)
library(scRNAseq)

mg.combined

sce.mg.PDUK <- as.SingleCellExperiment(mgPDUK)
sce.mg.PDUK <- scuttle::logNormCounts(sce.mg.PDUK)
# sce.mg.NBB <- as.SingleCellExperiment(mgNBB)
# sce.mg.NBB <- scuttle::logNormCounts(sce.mg.NBB)


# Search for human reference
searchDatasets(
    defineTextQuery("brain", partial = TRUE) &
    defineTextQuery("GRCh38", field = "genome") &
    defineTextQuery("9606", field = "taxonomy_id")
)[,c("name", "title", "genome", "taxonomy_id", "version")]

# Load reference
ref.dataset.name <- "zhong-prefrontal-2018"
ref.dataset.version <- "2023-12-22"
listPaths(ref.dataset.name, ref.dataset.version)
ref.dataset.path <- NA
sce.ref <- fetchDataset(ref.dataset.name, ref.dataset.version, ref.dataset.path)
# assay(sce.ref, "logcounts") <- assay(sce.ref, "counts")
sce.ref <- sce.ref[, colSums(counts(sce.ref)) > 0] # https://github.com/edward130603/BayesSpace/issues/28#issuecomment-717547233
sce.ref <- scuttle::logNormCounts(sce.ref)

# Run SingleR
sce.labeled <- SingleR(
    test = sce.mg.PDUK,
    # test = sce.mg.NBB,
    ref = sce.ref,
    labels = factor(sce.ref$cell_types)
)

cell.counts <- table(sce.labeled$labels) %>%
    data.frame() %>%
    mutate(
        perc = paste0(round(Freq / sum(Freq) * 100, 1), "%"),
    )
cell.counts

sample.name <- "PDUK"
# sample.name <- "NBB"
results.file <- file(file.path(paste0("out/labels_", sample.name, ".txt")), open = "wt")
writeLines(paste("Sample:", sample.name), results.file)
writeLines(paste("Total cells:", nrow(sce.labeled), "\n"), results.file)
write.table(cell.counts, results.file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
close(results.file)

