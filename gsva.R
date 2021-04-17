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





