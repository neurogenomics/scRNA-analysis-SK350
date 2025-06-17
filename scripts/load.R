library(Seurat)
library(ggplot2)

### Set vars ###################################################################
Pool1.data.dir <- "cellranger-output/count/SK350_Pool1_NBB/filtered_feature_bc_matrix"
Pool2.data.dir <- "cellranger-output/count/SK350_Pool2_PDUK/filtered_feature_bc_matrix"
out.dir <- "out"

################################################################################

# Prep combined SeuratObject
mgNBB.data <- Read10X(data.dir = Pool1.data.dir)
mgNBB <- CreateSeuratObject(counts = mgNBB.data, project = "SK350_Pool1_NBB", min.features = 200)
mgNBB

mgPDUK.data <- Read10X(data.dir = Pool2.data.dir)
mgPDUK <- CreateSeuratObject(counts = mgPDUK.data, project = "SK350_Pool2_PDUK", min.features = 200)
mgPDUK

mg.combined <- merge(mgNBB, y = mgPDUK, add.cell.ids = c("NBB", "PDUK"), project = "SK350_combined")
mg.combined

# QC
mg.combined[["percent.mt"]] <- PercentageFeatureSet(mg.combined, pattern = "^MT-")
VlnPlt <- VlnPlot(mg.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlt
ggsave(file.path(out.dir, "QC_vlnplot_minfeat200.png"), VlnPlt)


