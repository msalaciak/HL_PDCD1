library(Seurat)	
library(patchwork)	
library(Matrix)	
library(tidyverse)
library(plyr)
library(RCurl)	
library(cowplot)
library(clusterProfiler)
library("org.Dm.eg.db",character.only = TRUE)
library(org.Dm.eg.db)
library(DOSE)
library(plotly)
theme_set(theme_cowplot())
library(enrichplot)
library(reshape2)
library(EnhancedVolcano)
library(Nebulosa)
options(ggrepel.max.overlaps = Inf)
library(msigdbr)
library(magrittr)
library(GSVA)
library(scCATCH)


## load 10x datasets
#load datasets of each time point
timepoint_1.data <- Read10X(data.dir ="/home/matthew/datatransfer/mercier/151231/outs/filtered_gene_bc_matrices/GRCh38")
timepoint_2.data <- Read10X(data.dir ="/home/matthew/datatransfer/mercier/160411/outs/filtered_gene_bc_matrices/GRCh38")
timepoint_3.data <- Read10X(data.dir ="/home/matthew/datatransfer/mercier/161962/outs/filtered_gene_bc_matrices/GRCh38")
timepoint_4.data <- Read10X(data.dir ="/home/matthew/datatransfer/mercier/171094/outs/filtered_gene_bc_matrices/GRCh38")
timepoint_5.data <- Read10X(data.dir ="/home/matthew/datatransfer/mercier/171642/outs/filtered_gene_bc_matrices/GRCh38")
timepoint_6.data <- Read10X(data.dir ="/home/matthew/datatransfer/mercier/180251/outs/filtered_gene_bc_matrices/GRCh38")

#create seurat objects for each
timepoint_1 <-CreateSeuratObject(counts = timepoint_1.data, project = "pbmc timepoint 1", min.cells =3,min.features=200)
timepoint_2 <-CreateSeuratObject(counts = timepoint_2.data, project = "pbmc timepoint 2", min.cells =3,min.features=200)
timepoint_3 <-CreateSeuratObject(counts = timepoint_3.data, project = "pbmc timepoint 3", min.cells =3,min.features=200)
timepoint_4 <-CreateSeuratObject(counts = timepoint_4.data, project = "pbmc timepoint 4", min.cells =3,min.features=200)
timepoint_5 <-CreateSeuratObject(counts = timepoint_5.data, project = "pbmc timepoint 5", min.cells =3,min.features=200)
timepoint_6 <-CreateSeuratObject(counts = timepoint_6.data, project = "pbmc timepoint 6", min.cells =3,min.features=200)

#add metadata
#add metadata for each object
timepoint_1@meta.data$Timepoint <- "1"
timepoint_1@meta.data$Disease_status <- "HL"
timepoint_1@meta.data$Pembrolizumab <- "No"
timepoint_1@meta.data$iRAE <- "No"

timepoint_2@meta.data$Timepoint <- "2"
timepoint_2@meta.data$Disease_status <- "Remission"
timepoint_2@meta.data$Pembrolizumab <- "Yes"
timepoint_2@meta.data$iRAE <- "Yes"

timepoint_3@meta.data$Timepoint <- "3"
timepoint_3@meta.data$Disease_status <- "Remission"
timepoint_3@meta.data$Pembrolizumab <- "Yes"
timepoint_3@meta.data$iRAE <- "No"

timepoint_4@meta.data$Timepoint <- "4"
timepoint_4@meta.data$Disease_status <- "Relapse"
timepoint_4@meta.data$Pembrolizumab <- "No"
timepoint_4@meta.data$iRAE <- "No"

timepoint_5@meta.data$Timepoint <- "5"
timepoint_5@meta.data$Disease_status <- "Relapse"
timepoint_5@meta.data$Pembrolizumab <- "Yes"
timepoint_5@meta.data$iRAE <- "Yes"

timepoint_6@meta.data$Timepoint <- "6"
timepoint_6@meta.data$Disease_status <- "Relapse"
timepoint_6@meta.data$Pembrolizumab <- "No"
timepoint_6@meta.data$iRAE <- "No"


#remove doublets (used scrublet SEE python code if needed)

#tp1

tp1doublet <- read.table('/home/matthew/Research/HL_doublet/tp1-doublet.txt',  sep=',', header = TRUE)
tp1doublet <- tp1doublet[tp1doublet$predicted_doublets=='True',]
print("before doublet removal")
tibble(
  total_cells  = nrow(timepoint_1@meta.data),
  total_doublets = nrow(tp1doublet)
) 
timepoint_1cells <- colnames(timepoint_1)[!(colnames(timepoint_1) %in% tp1doublet$barcode)]
timepoint_1 <- subset(timepoint_1, cells=timepoint_1cells)
print("after doublet removal")
tibble(
  total_cells  = nrow(timepoint_1@meta.data),
  total_doublets = nrow(tp1doublet)
) 


#tp2
tp2doublet <- read.table('/home/matthew/Research/HL_doublet/tp2-doublet.txt',  sep=',', header = TRUE)
tp2doublet <- tp2doublet[tp2doublet$predicted_doublets=='True',]
print("before doublet removal")
tibble(
  total_cells  = nrow(timepoint_2@meta.data),
  total_doublets = nrow(tp2doublet)
) 

timepoint_2cells <- colnames(timepoint_2)[!(colnames(timepoint_2) %in% tp2doublet$barcode)]
timepoint_2 <- subset(timepoint_2, cells=timepoint_2cells)

print("after doublet removal")
tibble(
  total_cells  = nrow(timepoint_2@meta.data),
  total_doublets = nrow(tp2doublet)
) 


#tp3
tp3doublet <- read.table('/home/matthew/Research/HL_doublet/tp3-doublet.txt',  sep=',', header = TRUE)
tp3doublet <- tp3doublet[tp3doublet$predicted_doublets=='True',]
print("before doublet removal")
tibble(
  total_cells  = nrow(timepoint_3@meta.data),
  total_doublets = nrow(tp3doublet)
) 
timepoint_3cells <- colnames(timepoint_3)[!(colnames(timepoint_3) %in% tp3doublet$barcode)]
timepoint_3 <- subset(timepoint_3, cells=timepoint_3cells)

print("after doublet removal")
tibble(
  total_cells  = nrow(timepoint_3@meta.data),
  total_doublets = nrow(tp3doublet)
) 


#tp4
tp4doublet <- read.table('/home/matthew/Research/HL_doublet/tp4-doublet.txt',  sep=',', header = TRUE)
tp4doublet <- tp4doublet[tp4doublet$predicted_doublets=='True',]
print("before doublet removal")
tibble(
  total_cells  = nrow(timepoint_4@meta.data),
  total_doublets = nrow(tp4doublet)
) 
timepoint_4cells <- colnames(timepoint_4)[!(colnames(timepoint_4) %in% tp4doublet$barcode)]
timepoint_4 <- subset(timepoint_4, cells=timepoint_4cells)
print("after doublet removal")
tibble(
  total_cells  = nrow(timepoint_4@meta.data),
  total_doublets = nrow(tp4doublet)
) 

#tp5
tp5doublet <- read.table('/home/matthew/Research/HL_doublet/tp5-doublet.txt',  sep=',', header = TRUE)
tp5doublet <- tp5doublet[tp5doublet$predicted_doublets=='True',]
print("before doublet removal")
tibble(
  total_cells  = nrow(timepoint_5@meta.data),
  total_doublets = nrow(tp5doublet)
) 
timepoint_5cells <- colnames(timepoint_5)[!(colnames(timepoint_5) %in% tp5doublet$barcode)]
timepoint_5 <- subset(timepoint_5, cells=timepoint_5cells)
print("after doublet removal")
tibble(
  total_cells  = nrow(timepoint_5@meta.data),
  total_doublets = nrow(tp5doublet)
) 

#tp6
print("amount of doublets")
tp6doublet <- read.table('/home/matthew/Research/HL_doublet/tp6-doublet.txt', sep=',', header = TRUE)
tp6doublet <- tp6doublet[tp6doublet$predicted_doublets=='True',]
print("before doublet removal")
tibble(
  total_cells  = nrow(timepoint_6@meta.data),
  total_doublets = nrow(tp6doublet)
) 
timepoint_6cells <- colnames(timepoint_6)[!(colnames(timepoint_6) %in% tp6doublet$barcode)]
timepoint_6 <- subset(timepoint_6, cells=timepoint_6cells)
print("after doublet removal")
tibble(
  total_cells  = nrow(timepoint_6@meta.data),
  total_doublets = nrow(tp6doublet)
) 




#merge objects together

pbmc <- merge(timepoint_1, c(timepoint_2, timepoint_3, timepoint_4, 
                             timepoint_5,timepoint_6), add.cell.ids = c("Timepoint_1", "Timepoint_2", "Timepoint_3", 
                                                                "Timepoint_4", "Timepoint_5","Timepoint_6"))

#delete single objects and run clean up
rm(timepoint_1.data,timepoint_2.data,timepoint_3.data,timepoint_4.data,timepoint_5.data,timepoint_6.data,
   timepoint_1,timepoint_2,timepoint_3,timepoint_4,timepoint_5,timepoint_6, tp1doublet,tp2doublet,tp3doublet,tp4doublet,tp5doublet,
   tp6doublet,timepoint_1cells,timepoint_2cells,timepoint_3cells,timepoint_4cells,timepoint_5cells,timepoint_6cells)
gc()

#calculate qc

#umi
pbmc@meta.data$log10GenesPerUMI <- log10(pbmc$nFeature_RNA) / log10(pbmc$nCount_RNA)

#mitochondrial genes
pbmc <- PercentageFeatureSet(pbmc, "^MT-", col.name = "percent_mito")
mito_genes <- rownames(pbmc)[grep("^MT-", rownames(pbmc))]
head(mito_genes, 10)

#ribosomal proteins
pbmc <- PercentageFeatureSet(pbmc, "^RP[SL]", col.name = "percent_ribo")
ribo_genes <- rownames(pbmc)[grep("^RP[SL]", rownames(pbmc))]
head(ribo_genes, 10)


#hemoglobin genes
pbmc <- PercentageFeatureSet(pbmc, "^HB[^(P)]", col.name = "percent_hb")
hemoglob_genes <- rownames(pbmc)[grep("^HB[^(P)]", rownames(pbmc))]
head(hemoglob_genes, 10)

#plot qc
qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(pbmc, features = qc_feats, pt.size = 0.1, ncol=3) + 
  NoLegend() 

#make metadata df for plotting
metadata <- pbmc@meta.data
#cells column
metadata$cells <- rownames(metadata)
#rename column for plotting
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata$sample <-metadata$seq_folder

#num of cells
metadata %>% 
  ggplot(aes(x=Timepoint, fill=Timepoint)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# num of counts
metadata %>% 
  ggplot(aes(color=Timepoint, x=nUMI, fill= Timepoint)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=Timepoint, x=nGene, fill= Timepoint)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=Timepoint, y=log10(nGene), fill=Timepoint)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent_mito)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 300) +
  facet_wrap(~Timepoint)

#percent mito
metadata %>% 
  ggplot(aes(color=Timepoint, x=percent_mito, fill=Timepoint)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)

#percent hb
metadata %>% 
  ggplot(aes(color=Timepoint, x=percent_hb, fill=Timepoint)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = Timepoint, fill=Timepoint)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


# add some line to see our potential filter cutoff

#nFeature_RNA
VlnPlot(pbmc, features = "nFeature_RNA", pt.size = 0.1) + 
  NoLegend() + geom_hline(yintercept=4100, linetype="dashed", 
             color = "red", size=1)
#nCount_RNA
VlnPlot(pbmc, features = "nCount_RNA", pt.size = 0.1) + 
  NoLegend() + geom_hline(yintercept=40000, linetype="dashed", 
                          color = "red", size=1)

#percent_mito
VlnPlot(pbmc, features = "percent_mito", pt.size = 0.1) + 
  NoLegend() + geom_hline(yintercept=20, linetype="dashed", 
                          color = "red", size=1)

#percent_hb
VlnPlot(pbmc, features = "percent_hb", pt.size = 0.1) + 
  NoLegend() + geom_hline(yintercept=50, linetype="dashed", 
                          color = "red", size=1)

FeatureScatter(pbmc, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent_mito")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent_ribo")

counts_per_cell <- Matrix::colSums(pbmc)
mito_gene_read_counts <- Matrix::colSums(pbmc[mito_genes,])
pct_mito <- mito_gene_read_counts / counts_per_cell * 100
plot(sort(pct_mito), xlab = "cells sorted by percentage mitochondrial counts", ylab = "percentage mitochondrial counts")
abline(h=10, col='red')
abline(h=15, col='red')
abline(h=20, col='red')
#pre-filter stats
# Get some summary stats for the sample:
summary_stats_before_filtering <- tibble(
  total_cells  = nrow(pbmc@meta.data),
  mean_n_genes = mean(pbmc@meta.data$nFeature_RNA),
  sd_n_genes   = sd(pbmc@meta.data$nFeature_RNA),
  max_n_genes  = max(pbmc@meta.data$nFeature_RNA),
  min_n_genes  = min(pbmc@meta.data$nFeature_RNA),
  mean_UMI     = mean(pbmc@meta.data$nCount_RNA),
  sd_UMI       = sd(pbmc@meta.data$nCount_RNA),
  max_UMI      = max(pbmc@meta.data$nCount_RNA),
  min_UMI      = min(pbmc@meta.data$nCount_RNA)
) %>% mutate_all(function(x) round(x, 2))

summary_stats_before_filtering

# filter
pbmc <- subset(pbmc,  subset = nFeature_RNA > 200 & nFeature_RNA < 4100 & 
                 percent_hb < 20 & percent_mito < 10)


# Get some summary stats for the sample:
summary_stats_after_filtering <- tibble(
  total_cells  = nrow(pbmc@meta.data),
  mean_n_genes = mean(pbmc@meta.data$nFeature_RNA),
  sd_n_genes   = sd(pbmc@meta.data$nFeature_RNA),
  max_n_genes  = max(pbmc@meta.data$nFeature_RNA),
  min_n_genes  = min(pbmc@meta.data$nFeature_RNA),
  mean_UMI     = mean(pbmc@meta.data$nCount_RNA),
  sd_UMI       = sd(pbmc@meta.data$nCount_RNA),
  max_UMI      = max(pbmc@meta.data$nCount_RNA),
  min_UMI      = min(pbmc@meta.data$nCount_RNA)
) %>% mutate_all(function(x) round(x, 2))

summary_stats_after_filtering

#remove zero genes
counts <- GetAssayData(object = pbmc, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

pbmc.filtered <- CreateSeuratObject(filtered_counts, meta.data = pbmc@meta.data)

#check cell cycle
# Normalize the counts
seurat_phase <- NormalizeData(pbmc.filtered)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features =  cc.genes$g2m.genes, 
                                 s.features =  cc.genes$s.genes)

              

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase")
DimPlot(seurat_phase,
        reduction = "pca",
        split.by = "Phase")
#no large differences no need to regress out cell cycle but can do it to be safe
pbmc.filtered <- NormalizeData(pbmc.filtered)
pbmc.filtered <- ScaleData(pbmc.filtered)
pbmc.filtered <- CellCycleScoring(pbmc.filtered, 
                                 g2m.features =  cc.genes$g2m.genes, 
                                 s.features =  cc.genes$s.genes)


#clean up
rm(seurat_phase)
rm(nonzero)
rm(filtered_counts)
rm(counts)
rm(metadata)




#integrate data
Seurat_object.list <- SplitObject(pbmc.filtered, split.by = "Timepoint") %>%
  lapply(SCTransform, verbose = T, vars.to.regress = c("nFeature_RNA","nCount_RNA", "percent_mito","percent_ribo","S.Score", "G2M.Score"))

Seurat_object.features <- SelectIntegrationFeatures(object.list = Seurat_object.list,
                                                    nfeatures = 3000)

# some data wrangling
Seurat_object.list <- PrepSCTIntegration(object.list = Seurat_object.list,
                                         anchor.features = Seurat_object.features)

# identify anchors shared by the datasets
Seurat_object.anchors <- FindIntegrationAnchors(object.list = Seurat_object.list,
                                                normalization.method = "SCT", 
                                                anchor.features = Seurat_object.features)

# proceed with integration
pbmc.int <- IntegrateData(anchorset = Seurat_object.anchors,
                               normalization.method = "SCT")

#clean up
rm(Seurat_object.list)
rm(Seurat_object.features)
rm(Seurat_object.anchors)
rm(pbmc)
rm(pbmc.filtered)



#variable genes
var.genes <-VariableFeatures(pbmc.int)

#TCR
tra.remove <- var.genes[grep("^TRA[VDJC]", var.genes)]
trb.remove <- var.genes[grep("^TRB[VDJC]", var.genes)]
trd.remove <- var.genes[grep("^TRD[VDJC]", var.genes)]
trg.remove <- var.genes[grep("^TRG[VDJC]", var.genes)]

# IG
ig.remove <- var.genes[grep("^IG[HKL]", var.genes)]

combined.genes.remove <- c(tra.remove,trb.remove,trd.remove,trg.remove,ig.remove)

#remove from variable genes so we dont cluster on ig/tcr
VariableFeatures(pbmc.int) <- setdiff(var.genes, combined.genes.remove)

#PCA
pbmc.int <- RunPCA(object = pbmc.int)


DimHeatmap(pbmc.int, dims = 1:15, cells = 500, balanced = T)
PCAPlot(pbmc.int,
        split.by = "Timepoint")  

ElbowPlot(pbmc.int,ndims=40)
#more in depth to look for number of pc's
pct <- pbmc.int[["pca"]]@stdev / sum(pbmc.int[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

pbmc.int <- RunUMAP(pbmc.int, 
                             dims = 1:19,
                             reduction = "pca",n.neighbors=30L,min.dist = 0.3,spread=1,metric = 'manhattan')

# Plot UMAP                             
DimPlot(pbmc.int,group.by="Timepoint")   



pbmc.int <- FindNeighbors(object = pbmc.int, 
                                   dims = 1:19,k.param=20)


pbmc.int <- FindClusters(object = pbmc.int,
                                  resolution = c(0.5, 0.6, 0.8, 1,1.2),algorithm = 1)



#go through resolutions to see
Idents(object = pbmc.int) <- "integrated_snn_res.0.5"
DimPlot(pbmc.int,
        reduction = "umap",
        label = TRUE)

Idents(object = pbmc.int) <- "integrated_snn_res.0.6"
DimPlot(pbmc.int,
        reduction = "umap",
        label = TRUE)

Idents(object = pbmc.int) <- "integrated_snn_res.0.8"
DimPlot(pbmc.int,
        reduction = "umap",
        label = TRUE)

Idents(object = pbmc.int) <- "integrated_snn_res.1"
DimPlot(pbmc.int,
        reduction = "umap",
        label = TRUE)

Idents(object = pbmc.int) <- "integrated_snn_res.1.2"
DimPlot(pbmc.int,
        reduction = "umap",
        label = TRUE)


#res 0.8 is the winner. lets remove the other resolutions
pbmc.int[['integrated_snn_res.0.5']] <- NULL
pbmc.int[['integrated_snn_res.0.6']] <- NULL
pbmc.int[['integrated_snn_res.1']] <- NULL
pbmc.int[['integrated_snn_res.1.2']] <- NULL




DimPlot(pbmc.int,
        reduction = "umap",
        label = TRUE,split.by="Timepoint")

# clustering metrics

#heatshock score

# heatshock_features <- list(c('AAAS','AKT1S1','ANO1','ARPP21','ATM','ATP2A2','ATR','ATXN3','BAG1','BAG2','BAG3','BAG4','BAG5','CAMK2A','CAMK2B','CAMK2D','CAMK2G','CCAR2','CD34','CDKN1A','CHORDC1','CLPB','CREBBP','CRYAB','CXCL10','DAXX','DHX36','DNAJB1','DNAJB6','DNAJC2','DNAJC7','EIF2S1','EP300','FGF1','FKBP4','GSK3B','HDAC2','HIKESHI','HMOX1','HSBP1','HSBP1L1','HSF1','HSP90AA1','HSP90AA4P','HSP90AB1','HSP90AB2P','HSP90AB3P','HSP90AB4P','HSPA1A','HSPA1B','HSPA1L','HSPA6','HSPA8','HSPB8','HSPH1','HTRA2','IER5','IL1A','IRAK1','LYN','MAPK1','MAPK3','MAPKAPK2','MAPT','MLST8','MTOR','NDC1','NUP107','NUP133','NUP153','NUP155','NUP160','NUP188','NUP205','NUP210','NUP214','NUP35','NUP37','NUP42','NUP43','NUP50','NUP54','NUP58','NUP62','NUP85','NUP88','NUP93','NUP98','PDCD6','PDCL3','POLR2D','POM121','POM121C','PRKACA','PTGES3','RAE1','RANBP2','RBBP7','RPA1','RPA2','RPA3','RPTOR','SCARA5','SEC13','SEH1L','SIRT1','SLC52A3','SLU7','ST8SIA1','STAC','STUB1','SUMO1','TCIM','TFEC','THBS1','TPR','TRPV1','TRPV4','VCP','YWHAE'))
# 
# pbmc.int <- AddModuleScore(
#   object = pbmc.int,
#   features = heatshock_features,
#   ctrl = 5,
#   name = 'heatshock_features'
# )


metrics <-  c("nCount_RNA", "nFeature_RNA",  "percent_mito", "percent_ribo", "percent_hb")

FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

#look at pca's on umap
# Adding cluster label to center of cluster on UMAP
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(pbmc.int, 
                     vars = columns)
umap_label <- FetchData(pbmc.int, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

print(x = pbmc.int[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)



# look at genes
DefaultAssay(pbmc.int) <- "RNA"

pbmc.int <- NormalizeData(pbmc.int)

#T-cells
FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("CD3D","CD3E","CD3G","CD2"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#B-CELLS
FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("CD19", "MS4A1"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#CD4+ T-cells

FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("CD3D", "CD4"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("ZNF683", "FCGR3A"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#CD8+ T-cells
FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("CD3D", "CD8A"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#cDC
FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("FCER1A", "CST3"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


#pDC
FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("IL3RA", "SERPINF1"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#nk cells
FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("NCAM1", "NKG7","KLRB1","FCGR3A"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#megakaryocytes
FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("PPBP"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#CD14+
FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("CD14","LYZ"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#CD16+
FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("FCGR3A","MS4A7"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#macro
FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("MARCO", "ITGAM", "ADGRE1"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


FeaturePlot(pbmc.int, 
            reduction = "umap", 
            features = c("PDCD1", "HAVCR2", "CD244","EOMES","LAG3","CTLA4"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


#find DEGS
pbmc.int.markers <- FindAllMarkers(pbmc.int, only.pos = TRUE,
                                        min.pct = 0.1, logfc.threshold = 0.25)

pbmc.int.markers.full <-FindAllMarkers(pbmc.int, only.pos = FALSE,
                                       min.pct = 0.1, logfc.threshold = 0.05)


top10 <- pbmc.int.markers %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()


#to find out what cluster 17 and 23 are
cluster17 <- FindMarkers(pbmc.int, ident.1 = 17,ident.2=14,min.pct = .1,only.pos=TRUE,logfc.threshold = 0.15)
cluster17$gene<- rownames(cluster17)

cluster23 <- FindMarkers(pbmc.int,ident.1 = 23,min.pct = .1,only.pos=TRUE,logfc.threshold = 0.25)
cluster23$gene<- rownames(cluster23)

cluster12 <- FindMarkers(pbmc.int,ident.1 = 12,ident.2=c(5,6,8,15,19),min.pct = .1,only.pos=TRUE,logfc.threshold = 0.25)
cluster12$gene<- rownames(cluster12)
                                      

cluster19.conserved <-FindConservedMarkers(pbmc.int, ident.1 = 19,grouping.var = "iRAE", verbose = TRUE)
cluster19.conserved$gene <-rownames(cluster19.conserved)

#rename clusters
new.idents <- c('CD14 Monocyte-1','CD14 Monocyte-2','CD4 TREG',
                'CD14 Monocyte-3','CD4 Naive','CD8 TEM-3','CD8 TEFF-1',
                'CD14 Monocyte-4','CD8 TEM-2','CD16 Monocyte','NK CD56+',
                'B Cell','CD8 TEFF-2','CD4 TH17','CD14 Monocyte-5','CD8 TEM',
                'CD14 Monocyte-6','CD4 TEM','cDC2','CD8 TEX','NKT',
                'CD14 Monocyte-7','pDC','Platelet')

names(new.idents) <- levels(pbmc.int)
pbmc.int <- RenameIdents(pbmc.int, new.idents)
DimPlot(pbmc.int, reduction = "umap", label = TRUE,label.size = 3,repel=TRUE) + NoLegend()
pbmc.int[["new.ident"]] <-Idents(object = pbmc.int)

## CD8 subset
cd8.subset <- subset(pbmc.int, idents = c('CD8 TEM-2','CD8 TEFF-1','CD8 TEM','CD8 TEFF-2','CD8 TEM-3','CD8 TEX'))
cd8.subset <- DietSeurat(cd8.subset, assays = "RNA")

#redo integration
#integrate data
Seurat_object.list <- SplitObject(cd8.subset, split.by = "Timepoint") %>%
  lapply(SCTransform, verbose = T, vars.to.regress = c("nFeature_RNA","nCount_RNA", "percent_mito","percent_ribo","S.Score", "G2M.Score"))

Seurat_object.features <- SelectIntegrationFeatures(object.list = Seurat_object.list,
                                                    nfeatures = 3000)

# some data wrangling
Seurat_object.list <- PrepSCTIntegration(object.list = Seurat_object.list,
                                         anchor.features = Seurat_object.features)

# identify anchors shared by the datasets
Seurat_object.anchors <- FindIntegrationAnchors(object.list = Seurat_object.list,
                                                normalization.method = "SCT", 
                                                anchor.features = Seurat_object.features)

# proceed with integration
cd8.subset.int <- IntegrateData(anchorset = Seurat_object.anchors,
                          normalization.method = "SCT")

#clean up
rm(Seurat_object.list)
rm(Seurat_object.features)
rm(Seurat_object.anchors)
rm(cd8.subset)

# remove tcr genes for variable features
#variable genes
var.genes.cd8subset <-VariableFeatures(cd8.subset.int)

#TCR
tra.remove <- var.genes.cd8subset[grep("^TRA[VDJC]", var.genes.cd8subset)]
trb.remove <- var.genes.cd8subset[grep("^TRB[VDJC]", var.genes.cd8subset)]
trd.remove <- var.genes.cd8subset[grep("^TRD[VDJC]", var.genes.cd8subset)]
trg.remove <- var.genes.cd8subset[grep("^TRG[VDJC]", var.genes.cd8subset)]

# IG
ig.remove <- var.genes.cd8subset[grep("^IG[HKL]", var.genes.cd8subset)]

combined.genes.remove <- c(tra.remove,trb.remove,trd.remove,trg.remove,ig.remove)

#remove from variable genes so we dont cluster on ig/tcr
VariableFeatures(cd8.subset.int) <- setdiff(var.genes.cd8subset, combined.genes.remove)

cd8.subset.int.save<-cd8.subset.int

DefaultAssay(cd8.subset.int) <- "integrated"
#PCA
cd8.subset.int <- RunPCA(object = cd8.subset.int)


ElbowPlot(cd8.subset.int,ndims=40)
#more in depth to look for number of pc's
pct <- cd8.subset.int[["pca"]]@stdev / sum(cd8.subset.int[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

cd8.subset.int <- RunUMAP(cd8.subset.int, 
                          dims = 1:9,
                          reduction = "pca",n.neighbors=30L,min.dist = 0.1,spread=1,metric = 'manhattan')

# Plot UMAP                             
DimPlot(cd8.subset.int,group.by="Timepoint")   



cd8.subset.int <- FindNeighbors(object = cd8.subset.int, 
                                dims = 1:9,k.param=10)


cd8.subset.int <- FindClusters(object = cd8.subset.int,
                               resolution = c(0.5, 0.6, 0.8, 1,1.2),algorithm = 1)



#go through resolutions to see
Idents(object = cd8.subset.int) <- "integrated_snn_res.0.5"
DimPlot(cd8.subset.int,
        reduction = "umap",
        label = TRUE)

Idents(object = cd8.subset.int) <- "integrated_snn_res.0.6"
DimPlot(cd8.subset.int,
        reduction = "umap",
        label = TRUE)

Idents(object = cd8.subset.int) <- "integrated_snn_res.0.8"
DimPlot(cd8.subset.int,
        reduction = "umap",
        label = TRUE)

Idents(object = cd8.subset.int) <- "integrated_snn_res.1"
DimPlot(cd8.subset.int,
        reduction = "umap",
        label = TRUE)

Idents(object = cd8.subset.int) <- "integrated_snn_res.1.2"
DimPlot(cd8.subset.int,
        reduction = "umap",
        label = TRUE)

DimPlot(cd8.subset.int,
        reduction = "umap",
        label = TRUE,group.by="new.ident")

cd8.subset.int[['integrated_snn_res.0.5']] <- NULL
cd8.subset.int[['integrated_snn_res.1']] <- NULL
cd8.subset.int[['integrated_snn_res.1.2']] <- NULL

#look at markers to help decide
DefaultAssay(cd8.subset.int) <- "RNA"

cd8.subset.int <- NormalizeData(cd8.subset.int)

#T-cells
FeaturePlot(cd8.subset.int, 
            reduction = "umap", 
            features = c("PRF1"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

cd8.subset.int.markers <- FindAllMarkers(cd8.subset.int, only.pos = TRUE,
                                   min.pct = 0.1, logfc.threshold = 0.25)
cd8.subset.int.markers.full <- FindAllMarkers(cd8.subset.int, only.pos = TRUE,
                                         min.pct = 0.1, logfc.threshold = 0.05)

cluster7v10 <- FindMarkers(cd8.subset.int, ident.1 = 7,ident.2=10,min.pct = .1,only.pos=FALSE,logfc.threshold = 0.05)
cluster7v10$gene<- rownames(cluster7v10)

EnhancedVolcano(cluster7v10,
                lab = rownames(cluster7v10),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'cluster 7 vs clsuter 10',
                drawConnectors = TRUE,
                widthConnectors = 0.75)

DimPlot(pbmc.int,
        reduction = "umap",
        label = T,group.by="integrated_snn_res.0.8", repel = TRUE) + NoLegend()



### cd4/cd8/nk subset

## subset
tcell.subset <- subset(pbmc.int, idents = c('CD8 TEM-2','CD8 TEFF-1','CD8 TEM','CD8 TEFF-2','CD8 TEM-3','CD8 TEX','CD4 TREG','CD4 Naive' ,'NKT' , 'CD4 TH17' , 'CD4 TEM'))
tcell.subset <- DietSeurat(tcell.subset, assays = "RNA")

#redo integration
#integrate data
Seurat_object.list <- SplitObject(tcell.subset, split.by = "Timepoint") %>%
  lapply(SCTransform, verbose = T, vars.to.regress = c("nFeature_RNA","nCount_RNA", "percent_mito","percent_ribo","S.Score", "G2M.Score"))

Seurat_object.features <- SelectIntegrationFeatures(object.list = Seurat_object.list,
                                                    nfeatures = 4000)

# some data wrangling
Seurat_object.list <- PrepSCTIntegration(object.list = Seurat_object.list,
                                         anchor.features = Seurat_object.features)

# identify anchors shared by the datasets
Seurat_object.anchors <- FindIntegrationAnchors(object.list = Seurat_object.list,
                                                normalization.method = "SCT", 
                                                anchor.features = Seurat_object.features)

# proceed with integration
tcell.subset.int <- IntegrateData(anchorset = Seurat_object.anchors,
                                  normalization.method = "SCT")

#clean up
rm(Seurat_object.list)
rm(Seurat_object.features)
rm(Seurat_object.anchors)
rm(tcell.subset)

# remove tcr genes for variable features
#variable genes
var.genes.cd8subset <-VariableFeatures(tcell.subset.int)

#TCR
tra.remove <- var.genes.cd8subset[grep("^TRA[VDJC]", var.genes.cd8subset)]
trb.remove <- var.genes.cd8subset[grep("^TRB[VDJC]", var.genes.cd8subset)]
trd.remove <- var.genes.cd8subset[grep("^TRD[VDJC]", var.genes.cd8subset)]
trg.remove <- var.genes.cd8subset[grep("^TRG[VDJC]", var.genes.cd8subset)]

# IG
ig.remove <- var.genes.cd8subset[grep("^IG[HKL]", var.genes.cd8subset)]

combined.genes.remove <- c(tra.remove,trb.remove,trd.remove,trg.remove,ig.remove)

#remove from variable genes so we dont cluster on ig/tcr
VariableFeatures(tcell.subset.int) <- setdiff(var.genes.cd8subset, combined.genes.remove)



DefaultAssay(tcell.subset.int) <- "integrated"
#PCA

tcell.subset.int <- RunPCA(object = tcell.subset.int)


ElbowPlot(tcell.subset.int,ndims=40)
#more in depth to look for number of pc's
pct <- tcell.subset.int[["pca"]]@stdev / sum(tcell.subset.int[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

tcell.subset.int <- RunUMAP(tcell.subset.int, 
                            dims = 1:14,
                            reduction = "pca",n.neighbors=30L,min.dist = 0.1,spread=1,metric = 'manhattan',
                            n.epochs=1000
                            )

# Plot UMAP                             
DimPlot(tcell.subset.int,group.by="new.ident",label=T)   
DimPlot(tcell.subset.int,split.by="new.ident",label=T)   


tcell.subset.int <- FindNeighbors(object = tcell.subset.int, 
                                  dims = 1:14,k.param=20)


tcell.subset.int <- FindClusters(object = tcell.subset.int,
                                 resolution = c(0.4,0.5, 0.6, 0.8, 1,1.2,1.4),algorithm = 4)


#go through resolutions to see

Idents(object = tcell.subset.int) <- "integrated_snn_res.0.4"
DimPlot(tcell.subset.int,
        reduction = "umap",
        label = TRUE)

Idents(object = tcell.subset.int) <- "integrated_snn_res.0.5"
DimPlot(tcell.subset.int,
        reduction = "umap",
        label = TRUE)

Idents(object = tcell.subset.int) <- "integrated_snn_res.0.6"
DimPlot(tcell.subset.int,
        reduction = "umap",
        label = TRUE)

Idents(object = tcell.subset.int) <- "integrated_snn_res.0.8"
DimPlot(tcell.subset.int,
        reduction = "umap",
        label = TRUE)

Idents(object = tcell.subset.int) <- "integrated_snn_res.1"
DimPlot(tcell.subset.int,
        reduction = "umap",
        label = TRUE)

Idents(object = tcell.subset.int) <- "integrated_snn_res.1.2"
DimPlot(tcell.subset.int,
        reduction = "umap",
        label = TRUE)

Idents(object = tcell.subset.int) <- "integrated_snn_res.1.4"
DimPlot(tcell.subset.int,
        reduction = "umap",
        label = TRUE)

DimPlot(tcell.subset.int,
        reduction = "umap",
        label = TRUE,group.by="new.ident")
DimPlot(tcell.subset.int,
        reduction = "umap",
        label = TRUE,split.by="Timepoint")

tcell.subset.int[['integrated_snn_res.0.4']] <- NULL
tcell.subset.int[['integrated_snn_res.0.5']] <- NULL
tcell.subset.int[['integrated_snn_res.0.6']] <- NULL
tcell.subset.int[['integrated_snn_res.1.2']] <- NULL
tcell.subset.int[['integrated_snn_res.1.4']] <- NULL


#look at markers to help decide
DefaultAssay(tcell.subset.int) <- "RNA"

Idents(object = tcell.subset.int) <- "integrated_snn_res.0.8"
tcell.subset.int.markers_res0.8 <- FindAllMarkers(tcell.subset.int, only.pos = TRUE,
                                         min.pct = 0.1, logfc.threshold = 0.25)

tcell.subset.top10_res0.8  <- tcell.subset.int.markers_res0.8 %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()


Idents(object = tcell.subset.int) <- "integrated_snn_res.1"
tcell.subset.int.markers_res1.full <- FindAllMarkers(tcell.subset.int, only.pos = FALSE,
                                                  min.pct = 0.1, logfc.threshold = 0.05)

tcell.subset.top10_res1  <- tcell.subset.int.markers_res1 %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()





FeaturePlot(subset(tcell.subset.int,ident=c(5,10)), 
            reduction = "umap", 
            features = c("CD4","CD8A","CD8B","nCount_RNA"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

subset.expression<-AverageExpression(object = tcell.subset.int)

subset.rna<-as.data.frame(subset.expression$RNA)
subset.rna$gene <-rownames(subset.rna)


# going with res 1
DotPlot(
  tcell.subset.int,
  assay = NULL,
  unique(tcell.subset.top10_res1$gene),
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))


tcell.17 <- FindMarkers(tcell.subset.int, ident.1=17,only.pos = TRUE,
                                                     min.pct = 0.1, logfc.threshold = 0.10)
tcell.17$gene <-rownames(tcell.17)

new.idents.subset <- c('CD4 TH2', 'CD8 TEM-1','CD8 TEM-2','CD4 NAIVE','CD8 TCM-1','CD4 TH17','CD8 TEM-3','CD8 TEFF-1',
                       'CD4 TREG','CD4 TCM','NKT','CD8 TEX-1','CD8 TEFF-2','CD8 TEX-2','CD8 TCM-2','CD8 TEM-4','Unknown Cell')

names(new.idents.subset) <- levels(tcell.subset.int)
tcell.subset.int <- RenameIdents(tcell.subset.int, new.idents.subset)
tcell.subset.int[["new.ident"]] <-Idents(object = tcell.subset.int)
DimPlot(tcell.subset.int, reduction = "umap", label = TRUE,label.size = 3,repel=TRUE) + NoLegend()
