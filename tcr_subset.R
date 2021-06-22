#subset for timepoints for tcr analysis in scanpy/scripy


t1<- subset(x = tcell.subset.int, subset = Timepoint == "1")
t2<- subset(x = tcell.subset.int, subset = Timepoint == "2")
t3<- subset(x = tcell.subset.int, subset = Timepoint == "3")
t4<- subset(x = tcell.subset.int, subset = Timepoint == "4")
t5<- subset(x = tcell.subset.int, subset = Timepoint == "5")
t6<- subset(x = tcell.subset.int, subset = Timepoint == "6")


sceasy::convertFormat(t1, from="seurat", to="anndata",
                      outFile='t1.h5ad')
sceasy::convertFormat(t2, from="seurat", to="anndata",
                      outFile='t2.h5ad')
sceasy::convertFormat(t3, from="seurat", to="anndata",
                      outFile='t3.h5ad')
sceasy::convertFormat(t4, from="seurat", to="anndata",
                      outFile='t4.h5ad')
sceasy::convertFormat(t5, from="seurat", to="anndata",
                      outFile='t5.h5ad')
sceasy::convertFormat(t6, from="seurat", to="anndata",
                      outFile='t6.h5ad')
rm(t1)
rm(t2)
rm(t3)
rm(t4)
rm(t5)
rm(t6)


#import csv from scanpy/scirpy and merge to seurat obj

tcrInfo <- read.csv(file = 'tcell_subsets/tcr_analyzed.csv')
clonotypes <-select(tcrInfo,clonotype, X)

#this orders so it matches the 
myBarcode = rownames(tcell.subset.int@meta.data) #get barcode from seurat
clonotypes = clonotypes[match(myBarcode, clonotypes$X), ] #match the order of seurat barcode with you data
tcell.subset.int$clonotype = clonotypes$clonotype #put mutation into metadata
tcell.subset.int$tcell_barcode = clonotypes$X #put mutation into metadata





# 
# DimPlot(cd8.subset.cc, label=F,  cells.highlight= list(TCR_4,TCR_20,TCR_9,TCR_36,TCR_22,TCR_32,TCR_11,TCR_43,TCR_210,TCR_75,TCR_1475,TCR_374,TCR_59,TCR_3801,TCR_187,TCR_3873,TCR_881),split.by="Timepoint")  + 
#   scale_color_manual(labels = c("unselected","210_TCR","43_TCR","11_TCR","32_TCR","22_TCR","36_TCR","9_TCR","20_TCR","881_TCR","3873_TCR","187_TCR","3801_TCR","59_TCR","374_TCR","1475_TCR","75_TCR","4_TCR"), 
#                      values = c("darkgrey", "blue","darkolivegreen3","brown4","black","deeppink2","blueviolet","red","orange","darkgoldenrod3","darkseagreen","deepskyblue","cyan","coral","cadetblue4","burlywood","chartreuse1","darkslateblue")) +
#   labs(color = "Top TCR's Across Timepoint") 


#top everything
TCR_6 <- filter(tcell.subset.int@meta.data, clonotype == "6_TCR")$tcell_barcode
TCR_3 <- filter(tcell.subset.int@meta.data, clonotype == "3_TCR")$tcell_barcode
TCR_139 <- filter(tcell.subset.int@meta.data, clonotype == "139_TCR")$tcell_barcode
TCR_46 <- filter(tcell.subset.int@meta.data, clonotype == "46_TCR")$tcell_barcode
TCR_17 <- filter(tcell.subset.int@meta.data, clonotype == "17_TCR")$tcell_barcode
TCR_86 <- filter(tcell.subset.int@meta.data, clonotype == "86_TCR")$tcell_barcode
TCR_48 <- filter(tcell.subset.int@meta.data, clonotype == "48_TCR")$tcell_barcode
TCR_77 <- filter(tcell.subset.int@meta.data, clonotype == "77_TCR")$tcell_barcode
TCR_19 <- filter(tcell.subset.int@meta.data, clonotype == "19_TCR")$tcell_barcode
TCR_589 <- filter(tcell.subset.int@meta.data, clonotype == "589_TCR")$tcell_barcode
TCR_99 <- filter(tcell.subset.int@meta.data, clonotype == "99_TCR")$tcell_barcode
TCR_199 <- filter(tcell.subset.int@meta.data, clonotype == "199_TCR")$tcell_barcode
TCR_187 <- filter(tcell.subset.int@meta.data, clonotype == "187_TCR")$tcell_barcode
TCR_593 <- filter(tcell.subset.int@meta.data, clonotype == "593_TCR")$tcell_barcode
TCR_4021 <- filter(tcell.subset.int@meta.data, clonotype == "4021_TCR")$tcell_barcode
TCR_2663 <- filter(tcell.subset.int@meta.data, clonotype == "2663_TCR")$tcell_barcode
TCR_4087 <- filter(tcell.subset.int@meta.data, clonotype == "4087_TCR")$tcell_barcode

TCR_1068 <- filter(tcell.subset.int@meta.data, clonotype == "1068_TCR")$tcell_barcode
TCR_506 <- filter(tcell.subset.int@meta.data, clonotype == "506_TCR")$tcell_barcode
TCR_423 <- filter(tcell.subset.int@meta.data, clonotype == "423_TCR")$tcell_barcode

DimPlot(tcr.subset.int,label=F, repel=T,cells.highlight= list(TCR_6,TCR_3,TCR_139,TCR_46,TCR_17,TCR_86,TCR_48,TCR_77,TCR_19,
                                                                                  TCR_589,TCR_99,TCR_199,TCR_187,TCR_593,TCR_4021,TCR_2663,TCR_4087,
                                                                                  TCR_506),split.by="Timepoint")  + 
  scale_color_manual(labels = c("unselected","TCR_19","TCR_77","TCR_48","TCR_86","TCR_17","TCR_46","TCR_139","TCR_3","TCR_506","TCR_4087","TCR_2663","TCR_4021","TCR_593","TCR_187","TCR_199","TCR_99","TCR_589","TCR_6"), 
                     values = c("darkgrey", "blue","darkolivegreen3","brown4","blueviolet","deepskyblue","deeppink2","red","orange","cyan","black","darkgoldenrod3","cadetblue4","coral","darkseagreen","burlywood","chartreuse3","darkslateblue","coral3")) +
  labs(color = "Top TCR's Across Timepoint")  + theme(legend.position="bottom")




#exhausted tcr
TCR_4087 <- filter(tcell.subset.int@meta.data, clonotype == "4087_TCR")$tcell_barcode
TCR_1068 <- filter(tcell.subset.int@meta.data, clonotype == "1068_TCR")$tcell_barcode
TCR_3909 <- filter(tcell.subset.int@meta.data, clonotype == "3909_TCR")$tcell_barcode
TCR_7747 <- filter(tcell.subset.int@meta.data, clonotype == "7747_TCR")$tcell_barcode
TCR_8870 <- filter(tcell.subset.int@meta.data, clonotype == "8870_TCR")$tcell_barcode
TCR_1104 <- filter(tcell.subset.int@meta.data, clonotype == "1104_TCR")$tcell_barcode
TCR_8149 <- filter(tcell.subset.int@meta.data, clonotype == "8149_TCR")$tcell_barcode
TCR_7613 <- filter(tcell.subset.int@meta.data, clonotype == "7613_TCR")$tcell_barcode



 DimPlot(tcr.subset.int,label=F,cells.highlight=list(TCR_7613,TCR_8149,TCR_1104,TCR_8870,TCR_3909,TCR_1068,TCR_4087,TCR_7747),
         split.by="Timepoint",sizes.highlight=1.2) +
scale_color_manual(labels = c("unselected","TCR_7747","TCR_4087","TCR_1068","TCR_3909","TCR_8870","TCR_1104","TCR_8149","TCR_7613"),
                   values = c("lightgrey", "blue","chartreuse3","brown4","black","deeppink2","blueviolet","red","orange","darkgoldenrod3","darkseagreen","deepskyblue","cyan","coral","cadetblue4","burlywood","chartreuse1","darkslateblue")) +
 labs(color = "CD8 TEX Top TCR's Across Timepoint")  + theme(legend.position="bottom")


 ## markers
 cd8tem2.1v4 <- FindMarkers(subset(tcell.subset.int,idents=c('CD8 TEM-1','CD8 TEM-2','CD8 TCM-1','CD8 TEM-3',
                                                             'CD8 TEFF-1','CD8 TEX-1','CD8 TEFF-2','CD8 TEX-2',
                                                             'CD8 TCM-2','CD8 TEM-4')), ident.1 = "3",ident.2="4",group.by="Timepoint",
                         only.pos = FALSE, min.pct = .10,logfc.threshold = 0.05)
 
 cd8tem2.1v4$gene <-rownames(cd8tem2.1v4)
 
 EnhancedVolcano(cd8tem2.1v4,
                 lab = rownames(cd8tem2.1v4),
                 x = 'avg_log2FC',
                 y = 'p_val_adj',
                 pCutoff = 0.05,
                 FCcutoff = 0.75,
                 title = 'CD8 TEM-2 Timepoint 1 vs Timepoint 4',
                 drawConnectors = TRUE,
                 widthConnectors = 0.75)
 
 cd8teff2.1v4 <- FindMarkers(subset(tcell.subset.int,new.ident=="CD8 TEM-2"), ident.1 = 1,ident.2=4,group.by="Timepoint",
                            only.pos = FALSE, min.pct = .10,logfc.threshold = 0.05)
 
 cd8teff2.1v4$gene <-rownames(cd8teff2.1v4)
 
 EnhancedVolcano(cd8teff2.1v4,
                 lab = rownames(cd8teff2.1v4),
                 x = 'avg_log2FC',
                 y = 'p_val_adj',
                 pCutoff = 0.05,
                 FCcutoff = 0.75,
                 title = 'CD8 TEM-2 Timepoint 1 vs Timepoint 4',
                 drawConnectors = FALSE,
                 widthConnectors = 0.75)
 
 Idents(object = tcell.subset.int) <- "new.ident"
 
 tp1v4 <- FindAllMarkers(tcell.subset.int,
                             only.pos = TRUE, min.pct = .10,logfc.threshold = 0.25)
 

 tp1v4.top10  <- tp1v4 %>% 
   group_by(cluster) %>% 
   slice_max(avg_log2FC, n = 10) %>%
   ungroup()

 DotPlot(
   tcell.subset.int,
   assay = NULL,
   unique(tp1v4.top10$gene),
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
 
 DoHeatmap(subset(tcell.subset.int, downsample = 1000),  features= unique(tp1v4.top10$gene), size = 3)

 
 DimPlot(tcell.subset.int,split.by="Timepoint",label=F) + NoLegend()
 
 ## TCR DIVERSITY 
 clonotypes.div <- filter(tcrInfo, has_ir=="True")
 clonotypes.div <-select(clonotypes.div,clonotype, new.ident, Timepoint )
 
 tp1.clonotype.count <- nrow(filter(clonotypes.div, Timepoint==1))
 tp2.clonotype.count <- nrow(filter(clonotypes.div, Timepoint==2))
 tp3.clonotype.count <- nrow(filter(clonotypes.div, Timepoint==3))
 tp4.clonotype.count <- nrow(filter(clonotypes.div, Timepoint==4))
 tp5.clonotype.count <- nrow(filter(clonotypes.div, Timepoint==5))
 tp6.clonotype.count <- nrow(filter(clonotypes.div, Timepoint==6))
 
 
## cluster tcr only
 tcr.cells<-filter(tcrInfo, has_ir=="True")
 tcr.subset <-subset(tcell.subset.int, cells=tcr.cells$X)
 tcr.subset <-subset(tcr.subset, idents = c("Unknown Cell", "NKT"), invert = TRUE)

 tcr.subset <- DietSeurat(tcr.subset, assays = "RNA")
 
 #redo integration
 #integrate data
 Seurat_object.list <- SplitObject(tcr.subset, split.by = "Timepoint") %>%
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
 tcr.subset.int <- IntegrateData(anchorset = Seurat_object.anchors,
                                   normalization.method = "SCT")
 
 #clean up
 rm(Seurat_object.list)
 rm(Seurat_object.features)
 rm(Seurat_object.anchors)
 rm(tcr.subset)
 
 # remove tcr genes for variable features
 #variable genes
 var.genes.cd8subset <-VariableFeatures(tcr.subset.int)
 
 #TCR
 tra.remove <- var.genes.cd8subset[grep("^TRA[VDJC]", var.genes.cd8subset)]
 trb.remove <- var.genes.cd8subset[grep("^TRB[VDJC]", var.genes.cd8subset)]
 trd.remove <- var.genes.cd8subset[grep("^TRD[VDJC]", var.genes.cd8subset)]
 trg.remove <- var.genes.cd8subset[grep("^TRG[VDJC]", var.genes.cd8subset)]
 
 # IG
 ig.remove <- var.genes.cd8subset[grep("^IG[HKL]", var.genes.cd8subset)]
 
 combined.genes.remove <- c(tra.remove,trb.remove,trd.remove,trg.remove,ig.remove)
 
 #remove from variable genes so we dont cluster on ig/tcr
 VariableFeatures(tcr.subset.int) <- setdiff(var.genes.cd8subset, combined.genes.remove)
 
 
 
 DefaultAssay(tcr.subset.int) <- "integrated"
 #PCA
 
 tcr.subset.int <- RunPCA(object = tcr.subset.int)
 
 ElbowPlot(tcr.subset.int,ndims=40)

 tcr.subset.int <- RunUMAP(tcr.subset.int, 
                             dims = 1:20,
                             reduction = "pca")
 
 DimPlot(tcr.subset.int,group.by="new.ident",label=T,repel=T)  
 
 tcr.subset.int <- FindNeighbors(object = tcr.subset.int, 
                                   dims = 1:20)
 
 
 tcr.subset.int <- FindClusters(object = tcr.subset.int,
                                  resolution = 0.8,algorithm = 4)
 
 DimPlot(tcr.subset.int,label=T,repel=T)
 tcr.subset.int$celltype.cond <- paste(Idents(tcr.subset.int), tcr.subset.int$Timepoint, sep = "_tp_")
 
 DefaultAssay(tcr.subset.int) <- "SCT"
 

 
 
 DefaultAssay(tcr.subset.int) <- "RNA"
 
 #MODULE SCORE

 tcr.subset.int <- AddModuleScore(object = tcr.subset.int, features = OXIDATIVE_PHOSPHORYLATION, name = "OXIDATIVE_PHOSPHORYLATION",ctrl=20)
 tcr.subset.int <- AddModuleScore(object = tcr.subset.int, features = GLYCOLYSIS, name = "GLYCOLYSIS",ctrl=20)
 tcr.subset.int <- AddModuleScore(object = tcr.subset.int, features = MTOR, name = "MTOR",ctrl=20)
 tcr.subset.int <- AddModuleScore(object = tcr.subset.int, features = PI3K_AKT_MTOR, name = "PI3K_AKT_MTOR",ctrl=20)
 tcr.subset.int <- AddModuleScore(object = tcr.subset.int, features = FATTY_ACID_MET, name = "FATTY_ACID_METABOLISM",ctrl=20)
 


 tcr.subset.markers<- FindAllMarkers(tcr.subset.int, only.pos = TRUE,
                                                   min.pct = 0.1, logfc.threshold = 0.15,features=all.genes.nobcr)
 
 tcr.subset.markers.top15  <- tcr.subset.markers %>% 
    group_by(cluster) %>% 
    slice_max(avg_log2FC, n = 10) %>%
    ungroup()
 
 tcr.subset.int <- NormalizeData(tcr.subset.int)
 tcr.subset.int <- ScaleData(tcr.subset.int, features = rownames(tcr.subset.int))
 Idents(tcr.subset.int) <-"seurat_clusters"
 Idents(tcr.subset.int) <-"celltype.cond"
 Idents(object = tcr.subset.int) <- "new.ident.tcr"
 Idents(object = tcr.subset.int) <- "Timepoint"
 Idents(object = tcr.subset.int) <- "Disease_status"
 
 
 new.idents.subset <- c('CD4 TH2', "CD4 NAIVE","CD8 TEFF-3 ","CD8 TEFF-2","CD8 TEM-2","CD8 TEFF-1","CD4 TREG","CD4 TH17","CD8 TEFF-2","CD8 TEM-3","CD8 TEX-1")
 
 names(new.idents.subset) <- levels(tcr.subset.int)
 tcr.subset.int <- RenameIdents(tcr.subset.int, new.idents.subset)
 tcr.subset.int[["new.ident.tcr"]] <-Idents(object = tcr.subset.int)
 
 
tcr.subset.cluster1 <-FindMarkers(tcr.subset.int, ident.1="1",only.pos = FALSE,
                           min.pct = 0.1, logfc.threshold = 0.10)

tcr.subset.cluster1$gene <-rownames(tcr.subset.cluster1)


EnhancedVolcano(tcr.subset.cluster5.tp4v1,
                lab = rownames(tcr.subset.cluster5.tp4v1),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.75,
                title = 'CD8 TEM-2 Timepoint 4 vs all clusters',
                drawConnectors = TRUE,
                widthConnectors = 0.75)




#### TCR DIVERISTY
tcr.div <- data.frame(clonotype= tcr.subset.int$clonotype,ID=tcr.subset.int$new.ident.tcr,Timepoint=tcr.subset.int$Timepoint)

tcr.diversity
tcr.div.count.id <- tcr.div %>%                              # Applying group_by & summarise
   group_by(ID) %>%
   summarise(count = n_distinct(clonotype))

tcr.div.clonotype <- tcr.diversity  %>% group_by(ID, Timepoint,clonotype)
tcr.div.clonotype <- tcr.div.clonotype %>% tally()


