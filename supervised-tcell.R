library(STACAS)

#count gene per cell
sum(GetAssayData(object = tcell.pca, slot = "data")["CD8A",]>0)


tcell.pca <- subset(x = tcell.subset.int, idents = c("NKT"), invert = TRUE) 

pca.genes = c('CD3D','CD3G','CD4','CD8A','CD8B','GZMA','GZMK','GZMH','GZMB','GNLY',
              'PRF1','CCR7','SELL','CD69','CD28','CTLA4','PDCD1','LAG3','CD244','TBX21',
              'TIGIT','HAVCR2','IFNG','KLRG1','KLRD1','IL7R','FOS','JUN','TCF7',"CD44","CXCR6","MKI67","TCF7","EOMES","TOX")

tcell.pca <- DietSeurat(tcell.supervised.sub, assays = "RNA")


tcell.pca <- CellCycleScoring(tcell.pca, 
                                  g2m.features =  cc.genes$g2m.genes, 
                                  s.features =  cc.genes$s.genes)

tcell.pca <- PercentageFeatureSet(tcell.pca, "^MT-", col.name = "percent_mito")
tcell.pca <- PercentageFeatureSet(tcell.pca, "^RP[SL]", col.name = "percent_ribo")

Seurat_object.list <- SplitObject(tcell.pca, split.by = "Timepoint") %>%
  lapply(SCTransform, verbose = T, vars.to.regress = c("nFeature_RNA","nCount_RNA", "percent_mito","percent_ribo","S.Score", "G2M.Score"))

Seurat_object.features <- SelectIntegrationFeatures(object.list = Seurat_object.list,
                                                    nfeatures = 3000)


#REMOVE TCR/CC/RPS/MITO
var.genes <-Seurat_object.features

#TCR
tra.remove <- var.genes[grep("^TRA[VDJC]", var.genes)]
trb.remove <- var.genes[grep("^TRB[VDJC]", var.genes)]
trd.remove <- var.genes[grep("^TRD[VDJC]", var.genes)]
trg.remove <- var.genes[grep("^TRG[VDJC]", var.genes)]

# IG
ig.remove <- var.genes[grep("^IG[HKL]", var.genes)]

# mito/rb
mt.remove <- var.genes[grep("^MT-", var.genes)]
rb.remove <-var.genes[grep("^RP[SL]",var.genes)]

#unknown pseudogenes
mtr.remove <- var.genes[grep("^MTR", var.genes)]
rpl <-var.genes[grep("^RP[0-9]{1,}-",var.genes)]



combined.genes.remove <- c(tra.remove,trb.remove,trd.remove,trg.remove,ig.remove,cellCycle.symbol,mt.remove,rb.remove,mtr.remove,rpl)

#remove from variable genes so we dont cluster on ig/tcr
Seurat_object.features <- setdiff(var.genes, combined.genes.remove)



# some data wrangling
Seurat_object.list <- PrepSCTIntegration(object.list = Seurat_object.list,
                                         anchor.features = Seurat_object.features)

Seurat_object.list <- lapply(X = Seurat_object.list, FUN = RunPCA, features = Seurat_object.features)

# identify anchors shared by the datasets
Seurat_object.anchors <- FindIntegrationAnchors(object.list = Seurat_object.list,
                                                normalization.method = "SCT", 
                                                anchor.features = Seurat_object.features,
                                                dims = 1:30, reduction = "rpca")
# 
# Seurat_object.anchors <- FindIntegrationAnchors(object.list = Seurat_object.list,
#                                                 normalization.method = "SCT", 
#                                                 anchor.features = Seurat_object.features)

# proceed with integration
tcell.pca <- IntegrateData(anchorset = Seurat_object.anchors,
                                  normalization.method = "SCT", dims = 1:30)
# 
# tcell.pca <- IntegrateData(anchorset = Seurat_object.anchors,
#                            normalization.method = "SCT")
#clean up
rm(Seurat_object.list)
# rm(Seurat_object.features)
rm(Seurat_object.anchors)


### STACAS 

var.genes.n <- 800
var.genes.integrated.n <-500

cellCycle.symbol <- read.csv("cellCycle.symbol.DE.specific.170120.csv",as.is=T)$x
cellCycle.symbol = toupper(cellCycle.symbol)

for (i in 1:length(Seurat_object.list)) {
  Seurat_object.list[[i]] <- NormalizeData(Seurat_object.list[[i]], verbose = T)
  
  Seurat_object.list[[i]] <- FindVariableFeatures(Seurat_object.list[[i]], selection.method = "vst", 
                                        nfeatures = var.genes.n*2, verbose = T)
  
  mito.genes <- grep(pattern = "^mt-", rownames(Seurat_object.list[[i]]), value = TRUE)
  ribo.genes <- grep(pattern = "^Rp[ls]", rownames(Seurat_object.list[[i]]), value = TRUE)
  
  Seurat_object.list[[i]]@assays$RNA@var.features <- setdiff(Seurat_object.list[[i]]@assays$RNA@var.features, cellCycle.symbol)
  Seurat_object.list[[i]]@assays$RNA@var.features <- setdiff(Seurat_object.list[[i]]@assays$RNA@var.features, mito.genes)
  Seurat_object.list[[i]]@assays$RNA@var.features <- setdiff(Seurat_object.list[[i]]@assays$RNA@var.features, ribo.genes)
  Seurat_object.list[[i]]@assays$RNA@var.features <- head( Seurat_object.list[[i]]@assays$RNA@var.features, var.genes.n)
  
  
  
}


ndim=10
ref.anchors <- FindAnchors.STACAS(Seurat_object.list, dims=1:ndim, anchor.features=var.genes.integrated.n)

names <- names(Seurat_object.list)

plots <- PlotAnchors.STACAS(ref.anchors, obj.names=names)

g.cols <- 2
g.rows <- as.integer((length(plots)+2)/g.cols)
g <- do.call("arrangeGrob", c(plots, ncol=g.cols, nrow=g.rows))

plot(g)

ref.anchors.filtered <- FilterAnchors.STACAS(ref.anchors,dist.pct=0.8)

all.genes <- row.names(Seurat_object.list[[1]])
for (i in 2:length(Seurat_object.list)) {
  all.genes <- intersect(all.genes, row.names(Seurat_object.list[[i]]))
}

mySampleTree <- SampleTree.STACAS(ref.anchors.filtered)
print(mySampleTree)

ref.integrated <- IntegrateData(anchorset=ref.anchors.filtered, dims=1:ndim, features.to.integrate=all.genes,
                                sample.tree=mySampleTree, preserve.order=T)


ref.integrated <- ScaleData(ref.integrated, verbose = TRUE)

ref.integrated <- RunPCA(object = ref.integrated,features=pca.genes,approx=F)

ElbowPlot(ref.integrated,ndims=40)

ref.integrated <- RunUMAP(ref.integrated, 
                     dims = 1:10,
                     reduction = "pca")
DimPlot(ref.integrated)

DefaultAssay(ref.integrated) <-"RNA"
FeaturePlot(ref.integrated,features="CD8A",order=T)




DefaultAssay(tcell.pca) <-"RNA"

tcell.pca<-NormalizeData(tcell.pca)
tcell.pca<-ScaleData(tcell.pca)
tcell.pca <- RunPCA(object = tcell.pca,features=pca.genes,approx=F)



DimHeatmap(tcell.pca, dims = 1:15, cells = 500, balanced = T)
PCAPlot(tcell.pca,group.by="Disease_status")

ElbowPlot(tcell.pca,ndims=40)


tcell.pca <- RunUMAP(tcell.pca, 
                    dims = 1:10,
                    reduction = "pca")


tcell.pca <- FindNeighbors(object = tcell.pca, 
                          dims = 1:10,k.param=25)


tcell.pca <- FindClusters(object = tcell.pca,
                         resolution = 0.4)




# Plot UMAP                             
DimPlot(tcell.pca,label=T,repel=T)


DefaultAssay(tcell.pca) <-"RNA"
FeaturePlot(tcell.pca,features=c("PDCD1","GZMA","GZMH","GZMB","GZMK","PRF1","CTLA4","HAVCR2","LAG3","CD244","EOMES"),order=T)

FeaturePlot(tcell.pca,features=c("CD4"),order=T)

TCR_6 <- filter(tcell.pca@meta.data, clonotype == "6_TCR")$tcell_barcode
TCR_3 <- filter(tcell.pca@meta.data, clonotype == "3_TCR")$tcell_barcode
TCR_139 <- filter(tcell.pca@meta.data, clonotype == "139_TCR")$tcell_barcode
TCR_46 <- filter(tcell.pca@meta.data, clonotype == "46_TCR")$tcell_barcode
TCR_17 <- filter(tcell.pca@meta.data, clonotype == "17_TCR")$tcell_barcode
TCR_86 <- filter(tcell.pca@meta.data, clonotype == "86_TCR")$tcell_barcode
TCR_48 <- filter(tcell.pca@meta.data, clonotype == "48_TCR")$tcell_barcode
TCR_77 <- filter(tcell.pca@meta.data, clonotype == "77_TCR")$tcell_barcode
TCR_19 <- filter(tcell.pca@meta.data, clonotype == "19_TCR")$tcell_barcode
TCR_589 <- filter(tcell.pca@meta.data, clonotype == "589_TCR")$tcell_barcode
TCR_99 <- filter(tcell.pca@meta.data, clonotype == "99_TCR")$tcell_barcode
TCR_199 <- filter(tcell.pca@meta.data, clonotype == "199_TCR")$tcell_barcode
TCR_187 <- filter(tcell.pca@meta.data, clonotype == "187_TCR")$tcell_barcode
TCR_593 <- filter(tcell.pca@meta.data, clonotype == "593_TCR")$tcell_barcode
TCR_4021 <- filter(tcell.pca@meta.data, clonotype == "4021_TCR")$tcell_barcode
TCR_2663 <- filter(tcell.pca@meta.data, clonotype == "2663_TCR")$tcell_barcode
TCR_4087 <- filter(tcell.pca@meta.data, clonotype == "4087_TCR")$tcell_barcode

TCR_1068 <- filter(tcell.pca@meta.data, clonotype == "1068_TCR")$tcell_barcode
TCR_506 <- filter(tcell.pca@meta.data, clonotype == "506_TCR")$tcell_barcode
TCR_423 <- filter(tcell.pca@meta.data, clonotype == "423_TCR")$tcell_barcode

DimPlot(tcell.pca,label=F, repel=T,cells.highlight= list(TCR_6,TCR_3,TCR_139,TCR_46,TCR_17,TCR_86,TCR_48,TCR_77,TCR_19,
                                                              TCR_589,TCR_99,TCR_199,TCR_187,TCR_593,TCR_4021,TCR_2663,TCR_4087,
                                                              TCR_506),split.by="Timepoint")  + 
  scale_color_manual(labels = c("unselected","TCR_19","TCR_77","TCR_48","TCR_86","TCR_17","TCR_46","TCR_139","TCR_3","TCR_506","TCR_4087","TCR_2663","TCR_4021","TCR_593","TCR_187","TCR_199","TCR_99","TCR_589","TCR_6"), 
                     values = c("darkgrey", "blue","darkolivegreen3","brown4","blueviolet","deepskyblue","deeppink2","red","orange","cyan","black","darkgoldenrod3","cadetblue4","coral","darkseagreen","burlywood","chartreuse3","darkslateblue","coral3")) +
  labs(color = "Top TCR's Across Timepoint")  + theme(legend.position="bottom")

tcell.pca.markers <- FindAllMarkers(tcell.pca, only.pos = F,
                                   min.pct = 0.05, logfc.threshold = 0.15)


tcell.pca.markers.top10 <- tcell.pca.markers %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()
