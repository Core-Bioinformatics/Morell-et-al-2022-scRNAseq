library(Seurat)
set.seed(777)

min.cells = 0
min.features = 0

# read the cellranger output
mat = Read10X('10x-counts')

so = CreateSeuratObject(counts=mat, min.cells=min.cells, min.features=min.features)

mt.genes=grep("^MT-", rownames(so), value=FALSE)
rp.genes=grep("^RP[SL]", rownames(so), value=FALSE)
so[['percent.mt']] = PercentageFeatureSet(so, features=rownames(so)[mt.genes])
so[['percent.rp']] = PercentageFeatureSet(so, features=rownames(so)[rp.genes])
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4, pt.size=0)

# filter away cellsbased on QC distributions
so = subset(so, subset = nCount_RNA < 80000 & nFeature_RNA > 1000)

so = so[-c(mt.genes, rp.genes)]
so = SCTransform(so, return.only.var.genes=FALSE, verbose=FALSE)

# get list of abundant genes
n.abundant = 2000
matrix.so = GetAssayData(so, assay='SCT', slot='counts')
abundant.genes=rownames(matrix.so)[order(Matrix::rowSums(matrix.so), decreasing=TRUE)[1:n.abundant]]
length(abundant.genes)
abundant.genes[1:200]
rm(matrix.so)

# dim reductions
so <- RunPCA(so, features = abundant.genes, verbose = FALSE, npcs=30)
so <- RunUMAP(so, dims = 1:30, verbose = FALSE)
so = CellCycleScoring(so, cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)

# clustering using community detection on NN graph
so = FindNeighbors(so, reduction='pca', dims=1:30, k.param=5, verbose=FALSE)
so = FindClusters(so, algorithm=3, resolution=0.8, verbose=FALSE)

# get marker genes
lfc.thresh = log(2^1)
Idents(so) = so@meta.data$seurat_clusters
seurat.markers=FindAllMarkers(so, logfc.threshold=lfc.thresh, min.pct=0.0, test.use='roc', verbose=FALSE)
