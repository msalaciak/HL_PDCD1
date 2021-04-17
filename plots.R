#plots


DefaultAssay(tcell.subset.int) <-"RNA"
DefaultAssay(tcell.subset.int) <-"SCT"
DefaultAssay(tcell.subset.int) <-"integrated"

DimPlot(tcell.subset.int,label=T)

#overview
FeaturePlot(subset(tcell.subset.int,Timepoint==4), 
            features = c("GZMA","GZMH","GZMB","GZMK","PRF1"),
            cols = c("grey", "blue"))

plot_density(tcell.subset.int, c("GZMA","GZMH","GZMB","GZMK","PRF1","GNLY","KLRG1","KLRF1","KLRC1",
                              "KLRD1","IL7R","IFNG","SELL","CCR7","CD28","CD27","CD69","FOS","JUN","CD44",
                              "LAMP1","HLA-DRA","ICOS"))


DotPlot(tcell.subset.int, features = c("GZMA","GZMH","GZMB","GZMK","PRF1","GNLY","KLRG1","KLRF1","KLRC1",
                                       "KLRD1","IL7R","IFNG","SELL","CCR7","CD28","CD27","CD69","FOS","JUN","CD44",
                                       "LAMP1","HLA-DRA","ICOS")) + RotatedAxis()


FeaturePlot(tcell.subset.int, features=c("PDCD1"))

DoHeatmap(tcell.subset.int, features=c("PDCD1","EOMES","TOX","CD244","LAG3","CTLA4",'HAVCR2','TBX21','TCF7')) 


#genesets
OXIDATIVE_PHOSPHORYLATION<-list(c('ABCB7','ACAA1','ACAA2','ACADM','ACADSB','ACADVL','ACAT1','ACO2','AFG3L2','AIFM1','ALAS1','ALDH6A1','ATP1B1','ATP5F1A','ATP5F1B','ATP5F1C','ATP5F1D','ATP5F1E','ATP5MC1','ATP5MC2','ATP5MC3','ATP5ME','ATP5MF','ATP5MG','ATP5PB','ATP5PD','ATP5PF','ATP5PO','ATP6AP1','ATP6V0B','ATP6V0C','ATP6V0E1','ATP6V1C1','ATP6V1D','ATP6V1E1','ATP6V1F','ATP6V1G1','ATP6V1H','BAX','BCKDHA','BDH2','CASP7','COX10','COX11','COX15','COX17','COX4I1','COX5A','COX5B','COX6A1','COX6B1','COX6C','COX7A2','COX7A2L','COX7B','COX7C','COX8A','CPT1A','CS','CYB5A','CYB5R3','CYC1','CYCS','DECR1','DLAT','DLD','DLST','ECH1','ECHS1','ECI1','ETFA','ETFB','ETFDH','FDX1','FH','FXN','GLUD1','GOT2','GPI','GPX4','GRPEL1','HADHA','HADHB','HCCS','HSD17B10','HSPA9','HTRA2','IDH1','IDH2','IDH3A','IDH3B','IDH3G','IMMT','ISCA1','ISCU','LDHA','LDHB','LRPPRC','MAOB','MDH1','MDH2','MFN2','MGST3','MPC1','MRPL11','MRPL15','MRPL34','MRPL35','MRPS11','MRPS12','MRPS15','MRPS22','MRPS30','MTRF1','MTRR','MTX2','NDUFA1','NDUFA2','NDUFA3','NDUFA4','NDUFA5','NDUFA6','NDUFA7','NDUFA8','NDUFA9','NDUFAB1','NDUFB1','NDUFB2','NDUFB3','NDUFB4','NDUFB5','NDUFB6','NDUFB7','NDUFB8','NDUFC1','NDUFC2','NDUFS1','NDUFS2','NDUFS3','NDUFS4','NDUFS6','NDUFS7','NDUFS8','NDUFV1','NDUFV2','NNT','NQO2','OAT','OGDH','OPA1','OXA1L','PDHA1','PDHB','PDHX','PDK4','PDP1','PHB2','PHYH','PMPCA','POLR2F','POR','PRDX3','RETSAT','RHOT1','RHOT2','SDHA','SDHB','SDHC','SDHD','SLC25A11','SLC25A12','SLC25A20','SLC25A3','SLC25A4','SLC25A5','SLC25A6','SUCLA2','SUCLG1','SUPV3L1','SURF1','TCIRG1','TIMM10','TIMM13','TIMM17A','TIMM50','TIMM8B','TIMM9','TOMM22','TOMM70','UQCR10','UQCR11','UQCRB','UQCRC1','UQCRC2','UQCRFS1','UQCRH','UQCRQ','VDAC1','VDAC2','VDAC3'))
GLYCOLYSIS<-list(c('ABCB6','ADORA2B','AGL','AGRN','AK3','AK4','AKR1A1','ALDH7A1','ALDH9A1','ALDOA','ALDOB','ALG1','ANG','ANGPTL4','ANKZF1','ARPP19','ARTN','AURKA','B3GALT6','B3GAT1','B3GAT3','B3GNT3','B4GALT1','B4GALT2','B4GALT4','B4GALT7','BIK','BPNT1','CACNA1H','CAPN5','CASP6','CD44','CDK1','CENPA','CHPF','CHPF2','CHST1','CHST12','CHST2','CHST4','CHST6','CITED2','CLDN3','CLDN9','CLN6','COG2','COL5A1','COPB2','CTH','CXCR4','CYB5A','DCN','DDIT4','DEPDC1','DLD','DPYSL4','DSC2','ECD','EFNA3','EGFR','EGLN3','ELF3','ENO1','ENO2','ERO1A','EXT1','EXT2','FAM162A','FBP2','FKBP4','FUT8','G6PD','GAL3ST1','GALE','GALK1','GALK2','GAPDHS','GCLC','GFPT1','GLCE','GLRX','GMPPA','GMPPB','GNE','GNPDA1','GOT1','GOT2','GPC1','GPC3','GPC4','GPR87','GUSB','GYS1','GYS2','HAX1','HDLBP','HK2','HMMR','HOMER1','HS2ST1','HS6ST2','HSPA5','IDH1','IDUA','IER3','IGFBP3','IL13RA1','IRS2','ISG20','KDELR3','KIF20A','KIF2A','LCT','LDHA','LDHC','LHPP','LHX9','MDH1','MDH2','ME1','ME2','MED24','MERTK','MET','MIF','MIOX','MPI','MXI1','NANP','NASP','NDST3','NDUFV3','NOL3','NSDHL','NT5E','P4HA1','P4HA2','PAM','PAXIP1','PC','PDK3','PFKFB1','PFKP','PGAM1','PGAM2','PGK1','PGLS','PGM2','PHKA2','PKM','PKP2','PLOD1','PLOD2','PMM2','POLR3K','PPFIA4','PPIA','PPP2CB','PRPS1','PSMC4','PYGB','PYGL','QSOX1','RARS1','RBCK1','RPE','RRAGD','SAP30','SDC1','SDC2','SDC3','SDHC','SLC16A3','SLC25A10','SLC25A13','SLC35A3','SLC37A4','SOD1','SOX9','SPAG4','SRD5A3','STC1','STC2','STMN1','TALDO1','TFF3','TGFA','TGFBI','TKTL1','TPBG','TPI1','TPST1','TSTA3','TXN','UGP2','VCAN','VEGFA','VLDLR','XYLT2','ZNF292'))
MTOR<-list(c('ABCF2','ACACA','ACLY','ACSL3','ACTR2','ACTR3','ADD3','ADIPOR2','AK4','ALDOA','ARPC5L','ASNS','ATP2A2','ATP5MC1','ATP6V1D','AURKA','BCAT1','BHLHE40','BTG2','BUB1','CACYBP','CALR','CANX','CCNF','CCNG1','CCT6A','CD9','CDC25A','CDKN1A','CFP','COPS5','CORO1A','CTH','CTSC','CXCR4','CYB5B','CYP51A1','DAPP1','DDIT3','DDIT4','DDX39A','DHCR24','DHCR7','DHFR','EBP','EDEM1','EEF1E1','EGLN3','EIF2S2','ELOVL5','ELOVL6','ENO1','EPRS1','ERO1A','ETF1','FADS1','FADS2','FDXR','FGL2','FKBP2','G6PD','GAPDH','GBE1','GCLC','GGA2','GLA','GLRX','GMPS','GOT1','GPI','GSK3B','GSR','GTF2H1','HK2','HMBS','HMGCR','HMGCS1','HPRT1','HSP90B1','HSPA4','HSPA5','HSPA9','HSPD1','HSPE1','IDH1','IDI1','IFI30','IFRD1','IGFBP5','IMMT','INSIG1','ITGB2','LDHA','LDLR','LGMN','LTA4H','M6PR','MAP2K3','MCM2','MCM4','ME1','MLLT11','MTHFD2','MTHFD2L','NAMPT','NFIL3','NFKBIB','NFYC','NIBAN1','NMT1','NUFIP1','NUP205','NUPR1','P4HA1','PDAP1','PDK1','PFKL','PGK1','PGM1','PHGDH','PIK3R3','PITPNB','PLK1','PLOD2','PNO1','PNP','POLR3G','PPA1','PPIA','PPP1R15A','PRDX1','PSAT1','PSMA3','PSMA4','PSMB5','PSMC2','PSMC4','PSMC6','PSMD12','PSMD13','PSMD14','PSME3','PSMG1','PSPH','QDPR','RAB1A','RDH11','RIT1','RPA1','RPN1','RRM2','RRP9','SC5D','SCD','SDF2L1','SEC11A','SERP1','SERPINH1','SHMT2','SKAP2','SLA','SLC1A4','SLC1A5','SLC2A1','SLC2A3','SLC37A4','SLC6A6','SLC7A11','SLC7A5','SLC9A3R1','SORD','SQLE','SQSTM1','SRD5A1','SSR1','STARD4','STC1','STIP1','SYTL2','TBK1','TCEA1','TES','TFRC','TM7SF2','TMEM97','TOMM40','TPI1','TRIB3','TUBA4A','TUBG1','TXNRD1','UBE2D3','UCHL5','UFM1','UNG','USO1','VLDLR','WARS1','XBP1','YKT6'))
PI3K_AKT_MTOR <-list(c('ACACA','ACTR2','ACTR3','ADCY2','AKT1','AKT1S1','AP2M1','ARF1','ARHGDIA','ARPC3','ATF1','CAB39','CAB39L','CALR','CAMK4','CDK1','CDK2','CDK4','CDKN1A','CDKN1B','CFL1','CLTC','CSNK2B','CXCR4','DAPP1','DDIT3','DUSP3','E2F1','ECSIT','EGFR','EIF4E','FASLG','FGF17','FGF22','FGF6','GNA14','GNGT1','GRB2','GRK2','GSK3B','HRAS','HSP90B1','IL2RG','IL4','IRAK4','ITPR2','LCK','MAP2K3','MAP2K6','MAP3K7','MAPK1','MAPK10','MAPK8','MAPK9','MAPKAP1','MKNK1','MKNK2','MYD88','NCK1','NFKBIB','NGF','NOD1','PAK4','PDK1','PFN1','PIK3R3','PIKFYVE','PIN1','PITX2','PLA2G12A','PLCB1','PLCG1','PPP1CA','PPP2R1B','PRKAA2','PRKAG1','PRKAR2A','PRKCB','PTEN','PTPN11','RAC1','RAF1','RALB','RIPK1','RIT1','RPS6KA1','RPS6KA3','RPTOR','SFN','SLA','SLC2A1','SMAD2','SQSTM1','STAT2','TBK1','THEM4','TIAM1','TNFRSF1A','TRAF2','TRIB3','TSC2','UBE2D3','UBE2N','VAV3','YWHAB'))


tcell.subset.int <- AddModuleScore(object = tcell.subset.int, features = OXIDATIVE_PHOSPHORYLATION, name = "OXIDATIVE_PHOSPHORYLATION")
tcell.subset.int <- AddModuleScore(object = tcell.subset.int, features = GLYCOLYSIS, name = "GLYCOLYSIS")
tcell.subset.int <- AddModuleScore(object = tcell.subset.int, features = MTOR, name = "MTOR")
tcell.subset.int <- AddModuleScore(object = tcell.subset.int, features = PI3K_AKT_MTOR, name = "PI3K_AKT_MTOR")




plot_density(subset(tcell.subset.int,Timepoint==4), features=c("MTOR1","PI3K_AKT_MTOR1","OXIDATIVE_PHOSPHORYLATION1","GLYCOLYSIS1"))


p<-plot_density(subset(tcell.subset.int,Timepoint==4),features=c("MTOR1","PI3K_AKT_MTOR1","OXIDATIVE_PHOSPHORYLATION1","GLYCOLYSIS1"))
p1<-plot_density(subset(tcell.subset.int,Timepoint==4),features=c("GZMA","GZMH","GZMB","GZMK","PRF1","GNLY","KLRG1","KLRF1","KLRC1","KLRD1","IFNG",'MKI67'))
p2<-plot_density(subset(tcell.subset.int,Timepoint==4),features=c("IL7R","SELL","CCR7","CD28","CD27","CD44"))
p3<-plot_density(subset(tcell.subset.int,Timepoint==4),features=c("CD69","FOS","JUN","LAMP1","HLA-DRA","ICOS")) 
p4<-plot_density(subset(tcell.subset.int,Timepoint==4),features=c("PDCD1","EOMES","TOX","CD244","LAG3","CTLA4",'HAVCR2','TBX21','TCF7')) 
p + plot_annotation(
  title = 'Timepoint 4')
p1 + plot_annotation(
  title = 'Timepoint 4')
p2 + plot_annotation(
  title = 'Timepoint 4')
p3 + plot_annotation(
  title = 'Timepoint 4')
p4 + plot_annotation(
  title = 'Timepoint 4')

#dimplots saved
pub<-DimPlot(tcell.subset.int,reduction = "umap",label = TRUE ,repel = TRUE, label.size = 2.9) + theme_classic() + NoLegend() 


ggsave("dim_tcell.png",width = 8, height = 5,dpi = 300)



# SaveH5Seurat(tcell.subset.int, filename = "tcells.h5Seurat")
# Convert("tcells.h5Seurat", dest = "h5ad")


#table

tcellprop <- read.csv(file = 'tcellprop.csv')
pbmcprop <- read.csv(file = 'pbmccount.csv')


table1<- gt(data = tcellprop) %>% tab_header(
  title = "Amount of Cells Per Timepoint"
  
)%>%
  tab_stubhead(label = "Cell Identity") %>%
  cols_label(
    N.Cells = html("N Cells"),
    X..Cells = html("% Cells"),
    N.Cells.1 = html("N Cells"),
    X..Cells.1 = html("% Cells"),
    N.Cells.2 = html("N Cells"),
    X..Cells.2 = html("% Cells"),
    N.Cells.3 = html("N Cells"),
    X..Cells.3 = html("% Cells"),
    N.Cells.4 = html("N Cells"),
    X..Cells.4 = html("% Cells"),
    N.Cells.5 = html("N Cells"),
    X..Cells.5 = html("% Cells"),
  ) %>%
  tab_spanner(
    label = "1",  
    columns = vars(N.Cells, X..Cells),
    
  )  %>%
  tab_spanner(
    label = "3" ,
    columns = vars(N.Cells.2, X..Cells.2)
    
  ) %>%
  tab_spanner(
    label = "4", 
    columns = vars(N.Cells.3, X..Cells.3)
    
  )  %>%
  tab_spanner(
    label = "5",
    columns = vars(N.Cells.4, X..Cells.4)
    
  ) %>%
  tab_spanner(
    label = "6",
    columns = vars(N.Cells.5, X..Cells.5)
    
  )  %>%
  tab_spanner(
    label = "2",
    columns = vars(N.Cells.1, X..Cells.1)
    
  )  %>%

  
  tab_row_group(
    group = "NKT and Unknown Cells ",
    rows = 16:17) %>%
  
  tab_row_group(
    group = "CD4 T Cells",
    rows = 1:5) %>%
  
  tab_row_group(
    group = "CD8 T Cell ",
    rows = 6:15) 
  
  tab_style(
    style = list(
      
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      
      rows = 6)
  )

table1 %>% gt_theme_538()

table1 %>% gt_theme_538() %>%
  gtsave(
    "tcelltable.png",expand=10,vwidth = 1628,
    vheight = 882
    
  )


table1<- gt(data = pbmcprop) %>% tab_header(
  title = "Amount of Cells Per Timepoint"

)%>%
  tab_stubhead(label = "Cell Identity") %>%
  cols_label(
    N.Cells = html("N Cells"),
    X..Cells = html("% Cells"),
    N.Cells.1 = html("N Cells"),
    X..Cells.1 = html("% Cells"),
    N.Cells.2 = html("N Cells"),
    X..Cells.2 = html("% Cells"),
    N.Cells.3 = html("N Cells"),
    X..Cells.3 = html("% Cells"),
    N.Cells.4 = html("N Cells"),
    X..Cells.4 = html("% Cells"),
    N.Cells.5 = html("N Cells"),
    X..Cells.5 = html("% Cells"),
  ) %>%
  tab_spanner(
    label = "1",
    columns = vars(N.Cells, X..Cells),

  )  %>%
  tab_spanner(
    label = "3" ,
    columns = vars(N.Cells.2, X..Cells.2)

  ) %>%
  tab_spanner(
    label = "4",
    columns = vars(N.Cells.3, X..Cells.3)

  )  %>%
  tab_spanner(
    label = "5",
    columns = vars(N.Cells.4, X..Cells.4)

  ) %>%
  tab_spanner(
    label = "6",
    columns = vars(N.Cells.5, X..Cells.5)

  )  %>%
  tab_spanner(
    label = "2",
    columns = vars(N.Cells.1, X..Cells.1)

  )  %>%
  tab_row_group(
    group = "PBMC with T Cell Subsets",
    rows = 3:20) %>%
  tab_row_group(
    group = "T Cell Summary",
    rows = 1:2)%>%
  tab_style(
    style = list(
      cell_fill(color = "#ACEACE"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(

      rows = c(2)
      )
  )%>%
  tab_style(
    style = list(

      cell_text(weight = "bold")
    ),
    locations = cells_body(

      rows = 21)
  )



table1 %>% gt_theme_538()

table1 %>% gt_theme_538() %>%
  gtsave(
    "pbmctable.png",expand=10,vwidth = 1628,
    vheight = 882
    
  )





gt_theme_538 <- function(data,...) {
  data %>%
    opt_all_caps()  %>%
    opt_table_font(
      font = list(
        google_font("Myriad Pro")
      )
    ) %>%
    tab_style(
      style = cell_borders(
        sides = "bottom", color = "transparent", weight = px(2)
      ),
      locations = cells_body(
        columns = TRUE,
        # This is a relatively sneaky way of changing the bottom border
        # Regardless of data size
        rows = nrow(data$`_data`)
      )
    )  %>% 
    tab_options(
      column_labels.background.color = "white",
      table.border.top.width = px(3),
      table.border.top.color = "transparent",
      table.border.bottom.color = "transparent",
      table.border.bottom.width = px(3),
      column_labels.border.top.width = px(3),
      column_labels.border.top.color = "transparent",
      column_labels.border.bottom.width = px(3),
      column_labels.border.bottom.color = "black",
      data_row.padding = px(3),
      source_notes.font.size = 12,
      table.font.size = 16,
      heading.align = "left",
      grand_summary_row.background.color = "#990000",
      
      ...
    ) 
}


