#get db ready for gsva


m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)

msigdbr_list = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name)


# function to take seurat object and make gsva object
RunscGSVA <- function(object = object, assay = "RNA", slot = "data", geneset = geneset, ...) {
  object <- as.matrix(GetAssayData(object = object, assay = assay, slot = slot))
  object <- GSVA::gsva(expr = object, gset.idx.list = geneset, ...)
  return(object)
}

tcell.subset.int

gsva_results <- gsva(
  cd8tex2avg,
  msigdbr_list,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = 15,
  # Maximum gene set size
  max.sz = 500,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = TRUE
)
head(gsva_results[, 1:6])

gsva_df <-gsva_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway")

pheatmap(gsva_results,
                                      annotation_col = c(1,2,3,4,5,6), # Add metadata labels!
                                      show_colnames = FALSE, # Don't show sample labels
                                      fontsize_row = 6 # Shrink the pathway labels a tad
)

GetAssayData(object = subset(tcell.subset.int,subset=Timepoint==1), assay = 'RNA', slot = 'data')
gsva_tcell_full<-RunscGSVA(object=tcell.subset.int,assay = 'SCT',slot='data', geneset=msigdbr_list)

##add to assay in tcell

myBarcode = rownames(test@meta.data) #get barcode from seurat
gsva_tcell_full_transpose = gsva_tcell_full_transpose[match(myBarcode, rownames(gsva_tcell_full_transpose)), ] #match the order of seurat barcode with you data
gsva_tcell_full_transpose$TimePoint =test$Timepoint #put mutation into metadata
gsva_tcell_full_transpose$new.ident =test$new.ident 



# test[['gsva']] <- CreateAssayObject(counts = gsva_tcell_full)
# DefaultAssay(test) <- "gsva"


test<-tcell.subset.int



test <- performGeneSetEnrichmentAnalysis(
  object = test,
  GMT_file = "/home/matthew/Research/HL_PD1/genesets/c6.all.v7.4.symbols.gmt",
  groups = 'new.ident',
  thresh_p_val = 0.05,
  thresh_q_val = 0.1,
)

gsva.results <- as.data.frame(test@misc)


#### GSEA
#set up msigdbr
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)

## if DF of DEGS is from findallmarkers use this
## assume that 1st column is ID
## 2nd column is fold change

## feature 1: numeric vector
cluster0.genelist <- filter(cdtem2.t1, cluster==7 , pct.1>=.10)[,2]

## feature 2: named vector
names(cluster0.genelist) <- as.character(filter(cd8.subset.markers, cluster==7 , pct.1>=.10 )[,7])

## feature 3: decreasing order
cluster0.genelist <- sort(cluster0.genelist, decreasing = TRUE)

## if DF of DEGS is from findmarkers use this
## assume that 1st column is ID
## 2nd column is fold change

## feature 1: numeric vector
cluster0.genelist <- cdtem2.4v1[,2]
cluster0.genelist <- filter(cdtem2.4v1, pct.1>=.10, p_val_adj >=0.05)[,2]

## feature 2: named vector
names(cluster0.genelist) <- as.character(cdtem2.4v1[,6])
names(cluster0.genelist) <- as.character(filter(cdtem2.4v1, pct.1>=.10,p_val_adj >=0.05)[,6])

## feature 3: decreasing order
cluster0.genelist <- sort(cluster0.genelist, decreasing = TRUE)

# gsea


em2 <- GSEA(cluster0.genelist, TERM2GENE = m_t2g)

# ridgeplot(em2,showCategory = 20) + ggtitle("CD8 Exhausted 2v5")

dotplot(em2,showCategory = 20,split=".sign")+ facet_grid(.~.sign)  + ggtitle("CD 8 TEM2 Timepoint 4v1")

#GO
gse <- gseGO(geneList=cluster0.genelist, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb =org.Hs.eg.db, 
             pAdjustMethod = "none")

dotplot(gse, showCategory=15, split=".sign") + facet_grid(.~.sign)  + ggtitle("CD 8 TEM-1 Timepoint 4 vs Timpoint 1 GSE-GO")

# gene.test <- names(cluster0.genelist)[abs(cluster0.genelist) >1]
gene.test.up <- names(cluster0.genelist)[cluster0.genelist >= 1]
gene.test.down <- names(cluster0.genelist)[cluster0.genelist <= (-1)]
gene.test.2 <-names(cluster0.genelist)

gene.test.up<-bitr(gene.test.up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene.test.down<-bitr(gene.test.down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene.test.2<-bitr(gene.test.2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

gse.up <- enrichGO(gene          = gene.test.up$ENTREZID,
                   universe      = gene.test.2$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
)
dotplot(gse.up, showCategory =30)  + ggtitle("CD 8 TEM-1 Timepoint 1 vs Timpoint 4 GO UP")

gse.down <- enrichGO(gene          = gene.test.down$ENTREZID,
                     universe      = gene.test.2$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
)
dotplot(gse.down, showCategory =30)  + ggtitle("CD 8 TEM-1 Timepoint 1 vs Timpoint 4 GO DOWN")


##### escape
GS.hallmark <- getGeneSets(library = "C2")

ES.seurat <- enrichIt(obj = sub, gene.sets = GS.hallmark, groups = 1000, cores = 2)

tcr.subset.int_test <- Seurat::AddMetaData(sub, ES.seurat)


multi_dittoPlot(tcr.subset.int_test, vars = c("REACTOME_GLYCOLYSIS",
                                              "PID_PI3KCI_AKT_PATHWAY",
                                              "PID_CD8_TCR_PATHWAY"), 
                group.by = "Timepoint", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

ES2 <- data.frame(tcr.subset.int_test[[]], Idents(tcr.subset.int_test))
colnames(ES2)[ncol(ES2)] <- "cluster"

test<-getSignificance(ES2, group = "cluster", fit = "ANOVA")

filter(test,FDR<=0.05)
