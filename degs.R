#degs for tcell subsets

# cd8 tem-2
tcell.subset.int$celltype.cond <- paste(Idents(tcell.subset.int), tcell.subset.int$Timepoint, sep = "_")

DefaultAssay(tcell.subset.int)<- "RNA"
Idents(object = tcell.subset.int) <- "celltype.cond"
cdtem2.4v1 <- FindMarkers(tcell.subset.int, ident.1 = "CD8 TEM-2_4", ident.2 = "CD8 TEM-2_1", verbose = TRUE,logfc.threshold = 0.05)
head(cdtem2.4v1, n = 15)
cdtem2.4v1$gene <- rownames(cdtem2.4v1)

Idents(object = tcell.subset.int) <- "Timepoint"
Idents(object = tcell.subset.int) <- "new.ident"

cdtem6 <- FindMarkers(subset(tcell.subset.int,Timepoint==6), ident.1 = "CD8 TEM-2", verbose = TRUE,logfc.threshold = 0.05)
head(cdtem6, n = 15)
cdtem6$gene <- rownames(cdtem6)


sub<-subset(tcell.subset.int,new.ident=="CD8 TEM-2")
Idents(sub) <- "Timepoint"
 

cd8tem2avg <- log1p(AverageExpression(sub)$RNA)
cd8tem2avg<-as.data.frame(cd8tem2avg)
cd8tem2avg$gene <- rownames(cd8tem2avg)


cd8tem2avg<- mutate(cd8tem2avg,
    tp1v4 = cd8tem2avg[,1] - cd8tem2avg[,4],
    tp4v1 = cd8tem2avg[,4] - cd8tem2avg[,1]
  )



#exhausted
cdtem1.test <- FindMarkers(subset(tcell.subset.int,Timepoint==1), ident.1 = "CD8 TEM-2",
                           verbose = TRUE,logfc.threshold = 0.05,features=all.genes.notcr)
head(cdtem1.test, n = 15)
cdtem1.test$gene <- rownames(cdtem1.test)

all.genes<-rownames(tcell.subset.int)
#TCR
tra.remove <- all.genes[grep("^TRA[VDJC]", all.genes)]
trb.remove <- all.genes[grep("^TRB[VDJC]", all.genes)]
trd.remove <- all.genes[grep("^TRD[VDJC]", all.genes)]
trg.remove <- all.genes[grep("^TRG[VDJC]", all.genes)]

# IG
ig.remove <- all.genes[grep("^IG[HKL]", all.genes)]

combined.genes.remove <- c(tra.remove,trb.remove,trd.remove,trg.remove,ig.remove)
all.genes.notcr <- setdiff(all.genes, combined.genes.remove)

