
t1<- subset(x = pbmc.int, idents="B Cell", subset = Timepoint == "1")
t2<- subset(x = pbmc.int, idents="B Cell", subset = Timepoint == "2")
t3<- subset(x = pbmc.int, idents="B Cell", subset = Timepoint == "3")
t4<- subset(x = pbmc.int, idents="B Cell", subset = Timepoint == "4")
t5<- subset(x = pbmc.int, idents="B Cell", subset = Timepoint == "5")
t6<- subset(x = pbmc.int, idents="B Cell", subset = Timepoint == "6")


sceasy::convertFormat(t1, from="seurat", to="anndata",
                      outFile='b11.h5ad')
sceasy::convertFormat(t2, from="seurat", to="anndata",
                      outFile='b22.h5ad')
sceasy::convertFormat(t3, from="seurat", to="anndata",
                      outFile='b3.h5ad')
sceasy::convertFormat(t4, from="seurat", to="anndata",
                      outFile='b4.h5ad')
sceasy::convertFormat(t5, from="seurat", to="anndata",
                      outFile='b5.h5ad')
sceasy::convertFormat(t6, from="seurat", to="anndata",
                      outFile='b6.h5ad')
rm(t1)
rm(t2)
rm(t3)
rm(t4)
rm(t5)
rm(t6)

bcell.subset<-subset(x = pbmc.int, idents="B Cell")

bcrInfo <- read.csv(file = 'bcell_subsets/bcr_analyzed.csv')
clonotypes <-select(bcrInfo,clonotype, X)

#this orders so it matches the 
myBarcode = rownames(bcell.subset@meta.data) #get barcode from seurat
clonotypes = clonotypes[match(myBarcode, clonotypes$X), ] #match the order of seurat barcode with you data
bcell.subset$clonotype = clonotypes$clonotype #put mutation into metadata
bcell.subset$bcell_barcode = clonotypes$X #put mutation into metadata

BCR_143 <- filter(bcell.subset@meta.data, clonotype == "143_BCR")$bcell_barcode
BCR_65 <- filter(bcell.subset@meta.data, clonotype == "65_BCR")$bcell_barcode
BCR_313 <- filter(bcell.subset@meta.data, clonotype == "313_BCR")$bcell_barcode
BCR_255 <- filter(bcell.subset@meta.data, clonotype == "255_BCR")$bcell_barcode
BCR_19 <- filter(bcell.subset@meta.data, clonotype == "19_BCR")$bcell_barcode

DimPlot(bcell.subset, label=F,  cells.highlight=list(BCR_143,BCR_65),split.by="Timepoint")
DimPlot(bcell.subset, label=F,  split.by="Timepoint")


