#degs for tcell subsets

# cd8 tem-2
tcell.subset.int$celltype.cond <- paste(Idents(tcell.subset.int), tcell.subset.int$Timepoint, sep = "_")

DefaultAssay(tcell.subset.int)<- "RNA"
DefaultAssay(tcell.subset.int)<- "SCT"
Idents(object = tcell.subset.int) <- "Timepoint"
Idents(object = tcell.subset.int) <- "new.ident"

cdteff.t6 <- FindMarkers(subset(tcell.subset.int,Timepoint==6), ident.1 = "CD8 TEFF-1", verbose = TRUE,logfc.threshold = 0.05,features=all.genes.notcr)
cdteff.t6$gene <- rownames(cdteff.t6)

top.cdteff1.t6 <- cdteff.t6 %>% 
  slice_max(avg_log2FC, n = 35)


#combined all top and bottoms

gene.expressions<- data.frame(top_teff1_t1 =top.teff1.t1$gene,
                              bottom_teff_t1=bottom.teff1.t1$gene,
                              top_teff1_t3 =top.cdteff1.t3$gene,
                              bottom_teff_t3=bottom.cdteff1.t3$gene,
                              top_teff1_t4=top.teff1.t4$gene,
                              bottom_teff1_t4=bottom.teff1.t4$gene,
                              top_teff1_t6=top.cdteff1.t6$gene,
                              bottom_teff1_t6=bottom.cdteff1.t6$gene,
                              top_tem2_t1=top.tem2.t1$gene,
                              bottom_tem2_t1=bottom.tem2.t1$gene,
                              top_tem2_t3=top.cdtem2.t3$gene,
                              bottom_tem2_t3=bottom.cdtem2.t3$gene,
                              top_tem2_t4=top.tem2.t4$gene,
                              bottom_tem2_t4=bottom.tem2.t4$gene,
                              top_tem2_t6=top.cdtem2.t6$gene,
                              bottom_tem2_t6=bottom.cdtem2.t6$gene,
                              top_tex1_t1=top.tex1.tp1$gene,
                              bottom_tex1_t1=bottom.tex1.tp1$gene,
                              top_tex1_t3=top.cdtex1.t3$gene,
                              bottom_tex1_t3=bottom.cdtex1.t3$gene,
                              top_tex1_t4=top.tex1.tp4$gene,
                              bottom_tex1_t4=bottom.tex1.tp4$gene,
                              top_tex1_t6=top.cdtex1.t6$gene,
                              bottom_tex1_t6=bottom.cdtex1.t6$gene,
                              top_tex2_t1=top.tex2.tp1$gene,
                              bottom_tex2_t1=bottom.tex2.tp1$gene,
                              top_tex2_t3=top.cdtex2.t3$gene,
                              bottom_tex_t3=bottom.cdtex2.t3$gene,
                              top_tex2_t4=top.tex2.tp4$gene,
                              bottom_tex2_t4=bottom.tex2.tp4$gene,
                              top_tex2_t6=top.cdtex2.t6$gene,
                              bottom_tex2_t6=bottom.cdtex2.t6$gene)





Idents(object = tcell.subset.int) <- "new.ident"
cdtem1<- FindMarkers(tcell.subset.int, ident.1 = "CD8 TEM-1", verbose = TRUE,logfc.threshold = 0.05,features=all.genes.notcr)
cdtem1$gene <- rownames(cdtem1)





Idents(object = tcell.subset.int) <- "celltype.cond"
cdtem2.t4.t1 <- FindMarkers(tcell.subset.int, ident.1 = "CD8 TEM-2_4", ident.2 = "CD8 TEM-2_1", 
                              verbose = TRUE,logfc.threshold = 0.05,min.pct = .10,features=all.genes.notcr)

cdtem2.t4.t1$gene <- rownames(cdtem2.t4.t1)



sub<-subset(tcell.subset.int,new.ident=="CD8 TEM-2")
Idents(sub) <- "Timepoint"
 

cd8tex2avg <- log1p(AverageExpression(sub)$RNA)
cd8tex2avg<-as.data.frame(cd8tex2avg)
cd8tex2avg$gene <- rownames(cd8tex2avg)




cd8teff1avg<- mutate(cd8teff1avg,
    tp1v4 = cd8teff1avg[,1] - cd8teff1avg[,4],
    tp4v1 = cd8teff1avg[,4] - cd8teff1avg[,1]
  )



#exhausted
cdtex1v2 <- FindMarkers(tcell.subset.int, ident.1 = "CD8 TEX-1" ,ident.2="CD8 TEX-2"
                        ,logfc.threshold = 0.05,features=all.genes.notcr)
head(cdtex1v2, n = 15)
cdtex1v2$gene <- rownames(cdtex1v2)


#conserved
Idents(tcell.subset.int) <- "Timepoint"
sub <-subset(tcell.subset.int, idents = c(1, 4))
Idents(sub) <- "new.ident"
cdtex2.conserved <- FindConservedMarkers(sub, ident.1 = "CD8 TEX-2", grouping.var = "Timepoint", verbose = TRUE,logfc.threshold=0.05,feautres=genes.notcr)
cdtex2.conserved$genes <-rownames(cdtex2.conserved)




all.genes<-rownames(pbmc.int)
#TCR
tra.remove <- all.genes[grep("^TRA[VDJC]", all.genes)]
trb.remove <- all.genes[grep("^TRB[VDJC]", all.genes)]
trd.remove <- all.genes[grep("^TRD[VDJC]", all.genes)]
trg.remove <- all.genes[grep("^TRG[VDJC]", all.genes)]

# IG
ig.remove <- all.genes[grep("^IG[HKL]", all.genes)]

combined.genes.remove <- c(tra.remove,trb.remove,trd.remove,trg.remove,ig.remove)
all.genes.nobcr <- setdiff(all.genes, combined.genes.remove)

### BCELLS
Idents(pbmc.int)<- "new.ident"
pbmc.int$celltype.cond <- paste(Idents(pbmc.int), pbmc.int$Timepoint, sep = "_")
Idents(pbmc.int)<- "celltype.cond"
DefaultAssay(pbmc.int) <-"RNA"


bcell.t1 <- FindMarkers(pbmc.int, ident.1 = "B Cell_1" ,logfc.threshold = 0.25,features = all.genes.nobcr)

head(bcell.t1, n = 15)
bcell.t1$gene <- rownames(bcell.t1)




bcell.t3 <- FindMarkers(pbmc.int, ident.1 = "B Cell_3" ,logfc.threshold = 0.25,features = all.genes.nobcr)

head(bcell.t3, n = 15)
bcell.t3$gene <- rownames(bcell.t3)


bcell.t4 <- FindMarkers(pbmc.int, ident.1 = "B Cell_4" ,logfc.threshold = 0.25,features = all.genes.nobcr)

head(bcell.t4, n = 15)
bcell.t4$gene <- rownames(bcell.t4)


bcell.t6 <- FindMarkers(pbmc.int, ident.1 = "B Cell_6" ,logfc.threshold = 0.25,features = all.genes.nobcr)

head(bcell.t6, n = 15)
bcell.t6$gene <- rownames(bcell.t6)

bcell.t1.t3 <- FindMarkers(pbmc.int, ident.1 = "B Cell_1",ident.2="B Cell_3" ,logfc.threshold = 0.25,features = all.genes.nobcr)

head(bcell.t1.t3, n = 15)
bcell.t1.t3$gene <- rownames(bcell.t1.t3)



bcell.t3.t4 <- FindMarkers(pbmc.int, ident.1 = "B Cell_3",ident.2="B Cell_4" ,logfc.threshold = 0.25,features = all.genes.nobcr)

head(bcell.t3.t4, n = 15)
bcell.t3.t4$gene <- rownames(bcell.t3.t4)



bcell.t4.t6 <- FindMarkers(pbmc.int, ident.1 = "B Cell_4",ident.2="B Cell_6" ,logfc.threshold = 0.25,features = all.genes.nobcr)

head(bcell.t4.t6, n = 15)
bcell.t4.t6$gene <- rownames(bcell.t4.t6)


#cd8 tem-2
cdtem2.t3.t6 <- FindMarkers(pbmc.int, ident.1 = "CD8 TEM-2_3",ident.2="CD8 TEM-2_6" ,logfc.threshold = 0.05,features = all.genes.notcr)

head(cdtem2.t3.t6, n = 15)
cdtem2.t3.t6$gene <- rownames(cdtem2.t3.t6)

## cd8 teff-1
cdteff1.t3.t6 <- FindMarkers(pbmc.int, ident.1 = "CD8 TEFF-1_3",ident.2="CD8 TEFF-1_6" ,logfc.threshold = 0.05,features = all.genes.notcr)

head(cdteff1.t3.t6, n = 15)
cdteff1.t3.t6$gene <- rownames(cdteff1.t3.t6)


##Pdc 
Idents(pbmc.int)<- "celltype.cond"
DefaultAssay(pbmc.int) <-"RNA"


pdc.t1.t2 <- FindMarkers(pbmc.int, ident.1 = "pDC_1",ident.2="pDC_2" ,logfc.threshold = 0.25,features = all.genes.nobcr)

head(pdc.t1.t2, n = 15)
pdc.t1.t2$gene <- rownames(pdc.t1.t2)



pdcd1.subset <-subset(x = tcell.subset.int, subset = PDCD1 > 1)
DimPlot(pdcd1.subset,label=T) + NoLegend()
table(Idents(pdcd1.subset))
table(Idents(pdcd1.subset), pdcd1.subset$Timepoint)


sub<-subset(tcr.subset.int,idents=c(3))
Idents(sub) <- "Timepoint"

cd8tem2.avg <- AverageExpression(sub, return.seurat = TRUE)

cd8teff3.tcr.heat<-DoHeatmap(cd8tem2.avg, features=c("GZMA","GZMH","GZMB","GZMK","PRF1","GNLY","KLRG1",
                                  "KLRD1","IFNG","IL7R","SELL","CD27","CD28","CD69","FOS","JUN",
                                  "PDCD1","EOMES","TIGIT","CD244","LAG3","CTLA4",'HAVCR2','TBX21'
), size = 6, 
          draw.lines = F) + scale_fill_viridis_c() +ggtitle("CD8 TEFF-3")

cd8tem2.heat <- cd8tem2.heat+ ggtitle("CD8 TEM-2")
cd8teff1.heat <- cd8teff1.heat+ ggtitle("CD8 TEFF-1")

cd8tem2.tcr.heat <- cd8tem2.tcr.heat+ ggtitle("CD8 TEM-2 TCR SUBSET")
cd8teff1.tcr.heat <- cd8teff1.tcr.heat+ ggtitle("CD8 TEFF-1 TCR SUBSET")


cd8tem2.heat  +
cd8teff1.heat 

cd8tem2.tcr.heat +
cd8teff1.tcr.heat

(cd8tem2.heat + cd8tem2.tcr.heat) /
(cd8teff1.heat +cd8teff1.tcr.heat )

cd8tex1.heat + cd8tex2.heat

cd8tex1.tcr.heat + cd8teff3.tcr.heat
