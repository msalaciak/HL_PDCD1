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
# TCR_4 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSQGTGYTNTEAFF")$t_barcode
# TCR_20 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRA:CASLNQAGTALIF")$t_barcode
# TCR_9 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRA:CASLNQAGTALIF;TRB:CASSQGTGYTNTEAFF")$t_barcode
# TCR_36 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSLTGTTYNEQFF")$t_barcode
# TCR_22 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CSARDSAASTDTQYF")$t_barcode
# TCR_32 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSPKDYNNEQFF")$t_barcode
# TCR_11 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSVAGEHEQFF")$t_barcode
# TCR_43 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSLGLRESEQFF")$t_barcode
# TCR_210 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSLTTGNEQFF")$t_barcode
# TCR_75 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSLLQGKINEQFF")$t_barcode
# TCR_1475 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASGPGTYGYTF")$t_barcode
# TCR_374 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSHSATGESYEQYF")$t_barcode
# TCR_59 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSLGGATDTQYF")$t_barcode
# TCR_3801 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRA:CAVSPMNTGFQKLVF;TRB:CASSPAEMNTEAFF")$t_barcode
# TCR_187 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CSARDPGLAGKWDTQYF")$t_barcode
# TCR_3873 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSPAEMNTEAFF")$t_barcode
# TCR_881 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSHGGGNYEQYF")$t_barcode
# 
# 
# DimPlot(cd8.subset.cc, label=F,  cells.highlight= list(TCR_4,TCR_20,TCR_9,TCR_36,TCR_22,TCR_32,TCR_11,TCR_43,TCR_210,TCR_75,TCR_1475,TCR_374,TCR_59,TCR_3801,TCR_187,TCR_3873,TCR_881),split.by="Timepoint")  + 
#   scale_color_manual(labels = c("unselected","210_TCR","43_TCR","11_TCR","32_TCR","22_TCR","36_TCR","9_TCR","20_TCR","881_TCR","3873_TCR","187_TCR","3801_TCR","59_TCR","374_TCR","1475_TCR","75_TCR","4_TCR"), 
#                      values = c("darkgrey", "blue","darkolivegreen3","brown4","black","deeppink2","blueviolet","red","orange","darkgoldenrod3","darkseagreen","deepskyblue","cyan","coral","cadetblue4","burlywood","chartreuse1","darkslateblue")) +
#   labs(color = "Top TCR's Across Timepoint") 


TCR_6 <- filter(tcell.subset.int@meta.data, clonotype == "7747_TCR")$tcell_barcode

 DimPlot(tcell.subset.int, label=F,  cells.highlight=c(TCR_6),split.by="Timepoint")




