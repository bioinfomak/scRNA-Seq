library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(sctransform)

# Load the PBMC dataset
dbm_control1 <- Read10X(data.dir = "/home/gskbioinfo122/Desktop/dbm_control1/filtered_gene_bc_matrices/mm10/")
# Initialize the Seurat object with the raw (non-normalized data).
dbm_control1 <- CreateSeuratObject(counts = dbm_control1, project = "ctrl", min.cells = 3, min.features = 200)
dbm_control1

# Load the PBMC dataset
dbm_control2 <- Read10X(data.dir = "/home/gskbioinfo122/Desktop/dbm_control2/filtered_gene_bc_matrices/mm10/")
# Initialize the Seurat object with the raw (non-normalized data).
dbm_control2 <- CreateSeuratObject(counts = dbm_control2, project = "ctrl", min.cells = 3, min.features = 200)
dbm_control2

# Load the PBMC dataset
dbm_control3 <- Read10X(data.dir = "/home/gskbioinfo122/Desktop/dbm_control3/filtered_gene_bc_matrices/mm10/")
# Initialize the Seurat object with the raw (non-normalized data).
dbm_control3 <- CreateSeuratObject(counts = dbm_control3, project = "ctrl", min.cells = 3, min.features = 200)
dbm_control3


dbdb_diabetic1 <- Read10X(data.dir =  "/home/gskbioinfo122/Desktop/dbdb_diabetic1/filtered_gene_bc_matrices/mm10/")
dbdb_diabetic1 <- CreateSeuratObject(counts = dbdb_diabetic1, project = "diabetic", min.cells = 3, min.features = 200)
dbdb_diabetic1


dbdb_diabetic2 <- Read10X(data.dir =  "/home/gskbioinfo122/Desktop/dbdb_diabetic2/filtered_gene_bc_matrices/mm10/")
dbdb_diabetic2 <- CreateSeuratObject(counts = dbdb_diabetic2, project = "diabetic", min.cells = 3, min.features = 200)
dbdb_diabetic2


dbdb_diabetic3 <- Read10X(data.dir =  "/home/gskbioinfo122/Desktop/dbdb_diabetic3/filtered_gene_bc_matrices/mm10/")
dbdb_diabetic3 <- CreateSeuratObject(counts = dbdb_diabetic3, project = "diabetic", min.cells = 3, min.features = 200)
dbdb_diabetic3

#Merge Two Datasets
dbmdbdb_combined <- merge(dbm_control1, y = c (dbm_control2, dbm_control3, dbdb_diabetic1, dbdb_diabetic2, dbdb_diabetic3), add.cell.ids = c("dbm1", "dbm2", "dbm3", "dbdb1", "dbdb2", "dbdb3"), project = "DR")
dbmdbdb_combined


# notice the cell names now have an added identifier
head(colnames(dbmdbdb_combined))

table(dbmdbdb_combined$orig.ident)


#QC Control
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
dbmdbdb_combined[["percent.mt"]] <- PercentageFeatureSet(dbmdbdb_combined, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(dbmdbdb_combined@meta.data, 5)
mt.genes <- rownames(dbmdbdb_combined )[grep("^MT-",rownames(dbmdbdb_combined ))]
C<-GetAssayData(object = dbmdbdb_combined , slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
dbmdbdb_combined  <- AddMetaData(dbmdbdb_combined , percent.mito, col.name = "percent.mito")
rb.genes <- rownames(dbmdbdb_combined )[grep("^RP[SL]",rownames(dbmdbdb_combined ))]
percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
dbmdbdb_combined <- AddMetaData(dbmdbdb_combined , percent.ribo, col.name = "percent.ribo")
dbmdbdb_combined[["percent.mt"]] <- PercentageFeatureSet(dbmdbdb_combined, pattern = "^MT-")
dbmdbdb_combined[["percent.rb"]] <- PercentageFeatureSet(dbmdbdb_combined, pattern = "^RP[SL]")
VlnPlot(dbmdbdb_combined, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
head(dbmdbdb_combined@meta.data, 5)


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
dbmdbdb_combined[["percent.mt"]] <- PercentageFeatureSet(dbmdbdb_combined, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(dbmdbdb_combined@meta.data, 5)


#plot QC
VlnPlot(dbmdbdb_combined , features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(dbmdbdb_combined , features = "nCount_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(dbmdbdb_combined, features =c('nFeature_RNA','nCount_RNA'))
VlnPlot(dbmdbdb_combined, features =c("nFeature_RNA", "nCount_RNA"), ncol = 2)
RidgePlot(dbmdbdb_combined, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


#plot scatter
FeatureScatter(dbmdbdb_combined , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(dbmdbdb_combined , feature1 = "nFeature_RNA", feature2 = "percent.mito")
FeatureScatter(dbmdbdb_combined , feature1="percent.ribo", feature2="nFeature_RNA")


#filter_out
data.filter1 <- subset(dbmdbdb_combined , subset = nFeature_RNA > 200  & nFeature_RNA < 8000 & percent.mito < 20)
data.filter2 <- subset(dbmdbdb_combined , subset = nFeature_RNA > 200  & nFeature_RNA < 2500 & percent.mito < 5)
VlnPlot(data.filter1, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(data.filter1, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(data.filter2, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(data.filter2, features = "nCount_RNA", pt.size = 0.1) + NoLegend()


table(Idents(dbmdbdb_combined ))
table(Idents(data.filter1))
table(Idents(data.filter2))


#Normalization of data
dbmdbdb_combined <- NormalizeData(dbmdbdb_combined, normalization.method = "LogNormalize", scale.factor = 10000)
dbmdbdb_combined <- NormalizeData(dbmdbdb_combined)
NormalizeData(dbmdbdb_combined)

data.filter2 <- NormalizeData(data.filter2, normalization.method = "LogNormalize", scale.factor = 10000)
data.filter2 <- NormalizeData(data.filter2)
NormalizeData(data.filter2)

#FindVariableFeatures
dbmdbdb_combined <- FindVariableFeatures(dbmdbdb_combined , nfeatures = 8000, mean.cutoff = c(0.0125,3), dispersion.cutoff = c(0.5, Inf))

length(dbmdbdb_combined @assays$RNA@var.features)
summary(dbmdbdb_combined @assays$RNA@meta.features)

# Identify the 20 most highly variable genes
db.markers  <- FindVariableFeatures(dbmdbdb_combined , selection.method = "vst", nfeatures = 8000)
dim(db.markers)
top20 <- head(VariableFeatures(dbmdbdb_combined ), 20)
top10 <- head(VariableFeatures(dbmdbdb_combined ), 10)



# plot variable features with and without labels
plot1 <- VariableFeaturePlot(dbmdbdb_combined )
plot2 <-LabelPoints(plot = plot1, points = top10, repel = T, xnudge = 0, ynudge = 0)
plot1+plot2



#Scaling the data
all.genes <- rownames(dbmdbdb_combined )
dbmdbdb_combined  <- ScaleData(dbmdbdb_combined , features = all.genes)
dbmdbdb_combined  <- ScaleData(dbmdbdb_combined )
dbmdbdb_combined  <- ScaleData(dbmdbdb_combined , vars.to.regress = "percent.mt")
dbmdbdb_combined  <- SCTransform(dbmdbdb_combined , features = all.genes)
dbmdbdb_combined  <- SCTransform(dbmdbdb_combined )
dbmdbdb_combined  <- SCTransform(dbmdbdb_combined , vars.to.regress = "percent.mt")




#Dimensanality Reduction
dbmdbdb_combined  <- RunPCA(dbmdbdb_combined , features = VariableFeatures(object = dbmdbdb_combined ))
print(dbmdbdb_combined [["pca"]], dims = 1:10, nfeatures = 10)
VizDimLoadings(dbmdbdb_combined , dims = 1:10, reduction = "pca")
DimPlot(dbmdbdb_combined , reduction = "pca")
DimHeatmap(dbmdbdb_combined , dims = 1:10, cells = 50, balanced = TRUE)
DimHeatmap(dbmdbdb_combined , dims = 1:10, cells = 500, balanced = TRUE)
PCHeatmap(dbmdbdb_combined ,  dim = 1:10, cells = 1000, balanced=TRUE)
DimHeatmap(dbmdbdb_combined, dims = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))


# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
dbmdbdb_combined <- JackStraw(dbmdbdb_combined, num.replicate = 100)
dbmdbdb_combined <- ScoreJackStraw(dbmdbdb_combined, dims = 1:20)
dbmdbdb_combined <- JackStraw(dbmdbdb_combined, reduction = "pca", assay = NULL, dims = 20, num.replicate = 100, prop.freq = 0.01, verbose = TRUE, maxit = 1000)

JackStrawPlot(dbmdbdb_combined, dims = 1:15)
JackStrawPlot(dbmdbdb_combined, dims = 1:15, cols = NULL, reduction = "pca", xmax = 0.1, ymax = 0.3)




#Dimensationality
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
ElbowPlot(dbmdbdb_combined)


#cluster
dbmdbdb_combined <- FindNeighbors(dbmdbdb_combined, dims = 1:14)
dbmdbdb_combined <- FindClusters(dbmdbdb_combined, resolution = 0.5)


# Look at cluster IDs of the first 5 cells
head(Idents(dbmdbdb_combined), 5)


#UMAP/tSNE
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
dbmdbdb_combined <- RunUMAP(dbmdbdb_combined
                        , dims = 1:10)
dbmdbdb_combined <- RunTSNE(dbmdbdb_combined, dims = 1:10 , perplexity= 30)



# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(dbmdbdb_combined, reduction = "umap")
DimPlot(dbmdbdb_combined, reduction = "tsne")


# find all markers of cluster 2
cluster3.markers <- FindMarkers(dbmdbdb_combined, ident.1 = 2, min.pct = 0.25)
head(cluster3.markers, n = 5)



# find markers for every cluster compared to all remaining cells, report only the positive ones
db.markers <- FindAllMarkers(dbmdbdb_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
db1.markers <- FindAllMarkers(dbmdbdb_combined, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.26)
db.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
db1.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
cluster0.markers <- FindMarkers(dbmdbdb_combined, ident.1 = 2, logfc.threshold = 0.26, test.use = "roc", only.pos = TRUE)

write.csv(db.markers, file = "db.markers.csv")


# you can plot raw counts as well
VlnPlot(dbmdbdb_combined, features = c("Clu"))
VlnPlot(dbmdbdb_combined, features = c("Clu", "Glul"), slot = "counts", log = TRUE)
VlnPlot(dbmdbdb_combined, features = c("Clu", "Glul", "Gng13",  "Sag", "Pdc", "Rho", "Gnat1", "Pcp2"))
FeaturePlot(dbmdbdb_combined, features = c("Clu", "Glul", "Gng13", "Sag", "Pdc", "Rho", "Pcp2", "Pde6g"))
FeaturePlot(dbmdbdb_combined, features = c("Rho", "Gngt1"), blend = TRUE)
RidgePlot(dbmdbdb_combined, features = c("Rho", "Gngt1", "Sag"), ncol = 2)



akimba.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(dbmdbdb_combined, features = top10$gene) + NoLegend()
DoHeatmap(subset(dbmdbdb_combined, downsample = 100), features = c("Clu", "Glul", "Gng13", "Sag", "Pdc", "Rho", "Pcp2", "Pde6g"), size = 3)



#cluster.ids
new.cluster.ids <- c("Muller Glial Cell","Cardiomyocytes", "Photo Receptor Cells", "Photo Receptor Cells", "Muller Glial Cell", "Bipolar Cells", "Horizontal Cells","Bipolar Cells", "Endothelial Cells", "Muller Glial Cell", "Photo Receptor Cells", "Muller Glial Cell", "Bipolar Cells", "Endothelial Cells", "Photo Receptor Cells", "Photo Receptor Cells", "Endothelial Cells")
names(new.cluster.ids) <- levels(dbmdbdb_combined)
data.filter2 <- RenameIdents(data.filter2, new.cluster.ids)
DimPlot(dbmdbdb_combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(dbmdbdb_combined, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()


#CellCycleScoring
dbmdbdb_combined <- CellCycleScoring (object = dbmdbdb_combined, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
VlnPlot(dbmdbdb_combined, features = c("S.Score","G2M.Score"))
FeaturePlot(data.filter2,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))

dbmdbdb_combined <- CellCycleScoring (object = dbmdbdb_combined, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
VlnPlot(dbmdbdb_combined,features = c("S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))

dbmdbdb_combined <- CellCycleScoring(dbmdbdb_combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(dbmdbdb_combined[[]])


# Visualize the distribution of cell cycle markers across
RidgePlot(dbmdbdb_combined, features = c("Clu", "Glul", "Gng13", "Sag", "Pdc", "Rho", "Pcp2", "Pde6g"), ncol = 2)


# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
dbmdbdb_combined <- RunPCA(dbmdbdb_combined, features = c(s.genes, g2m.genes))
DimPlot(dbmdbdb_combined)



dbmdbdb_combined <- ScaleData(dbmdbdb_combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(dbmdbdb_combined))


# Now, a PCA on the variable genes no longer returns components associated with cell cycle
dbmdbdb_combined <- RunPCA(dbmdbdb_combined, features = VariableFeatures(dbmdbdb_combined), nfeatures.print = 10)


# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
dbmdbdb_combined <- RunPCA(dbmdbdb_combined, features = c(s.genes, g2m.genes))
DimPlot(dbmdbdb_combined)





dbmdbdb_combined <- CellCycleScoring(dbmdbdb_combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(dbmdbdb_combined[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(dbmdbdb_combined, features = c("Rho", "Gngt1", "Pde6g", "Gnat1"), ncol = 2)




# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
dbmdbdb_combined <- RunPCA(dbmdbdb_combined, features = c(s.genes, g2m.genes))
DimPlot(dbmdbdb_combined)


dbmdbdb_combined <- ScaleData(dbmdbdb_combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(dbmdbdb_combined))


# Now, a PCA on the variable genes no longer returns components associated with cell cycle
dbmdbdb_combined <- RunPCA(dbmdbdb_combined, features = VariableFeatures(dbmdbdb_combined), nfeatures.print = 10)


# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
dbmdbdb_combined <- RunPCA(dbmdbdb_combined, features = c(s.genes, g2m.genes))
DimPlot(dbmdbdb_combined)



