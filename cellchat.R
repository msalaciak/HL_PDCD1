library(CellChat)

#get info from seurat for cellchat object
timepoint_1.subset <- subset(x = pbmc.int, subset = Timepoint == "1")
timepoint_2.subset <- subset(x = pbmc.int, subset = Timepoint == "2")
timepoint_3.subset <- subset(x = pbmc.int, subset = Timepoint == "3")
timepoint_4.subset <- subset(x = pbmc.int, subset = Timepoint == "4")
timepoint_5.subset <- subset(x = pbmc.int, subset = Timepoint == "5")
timepoint_6.subset <- subset(x = pbmc.int, subset = Timepoint == "6")



count_raw<- timepoint_6.subset@assays$RNA@data
count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)

write.table(count_norm, '/home/matthew/Research/HL_PD1/cellchatdb/cellphonedb_count.txt', sep='\t', quote=F)

meta_data <- cbind(rownames(timepoint_6.subset@meta.data), timepoint_6.subset@meta.data[,'new.ident', drop=F])
write.table(meta_data, '/home/matthew/Research/HL_PD1/cellchatdb/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F) 
