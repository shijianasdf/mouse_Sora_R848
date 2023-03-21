#' 单细胞分析(老鼠)
#' @author shijian
# Vehicle
{
  N.sce <- FastSeurat(count=NULL,
                      dir.name="/Volumes/HDD4/tianjin_supple/tianjin_singleCell/mRNA/N/mRNA/outs/filtered_feature_bc_matrix/",
                      obj=NULL,#min.cells=3,
                      species=c("human","mouse")[2],
                      min.features=0,
                      max.features=300000,
                      percent.mt.num=100,
                      plot=F,
                      pcSelect=30,
                      project="Vehicle",
                      nfeatures=2000,
                      vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                      all.scale = F,
                      npcs = 50,
                      resolution = 0.5,
                      harmony=F,
                      doublet=F,
                      isMarkers=F,
                      perplexity = 30,
                      cellCycle=F,
                      features=NULL,
                      filepath=NULL,
                      outdir="Results",
                      names="love")
  p1 <- N.sce$v
  p1
  p2 <- DimPlot(N.sce$sce, reduction = "umap",group.by="seurat_clusters",label.size = 4,repel=T,label=T)
  p3 <- DimPlot(N.sce$sce, reduction = "umap",group.by="orig.ident",label.size = 4,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(N.sce$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(N.sce$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(N.sce$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  N.sce <- FastSeurat(count=NULL,
                      dir.name="/Volumes/HDD4/tianjin_supple/tianjin_singleCell/mRNA/N/mRNA/outs/filtered_feature_bc_matrix/",
                      obj=NULL,#min.cells=3,
                      species=c("human","mouse")[2],
                      min.features=200,
                      max.features=6000,
                      percent.mt.num=50,
                      plot=F,
                      pcSelect=30,
                      project="Vehicle",
                      nfeatures=2000,
                      vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                      all.scale = F,
                      npcs = 50,
                      resolution = 0.5,
                      harmony=F,
                      doublet=T,
                      isMarkers=T,
                      perplexity = 30,
                      cellCycle=T,
                      features=NULL,
                      filepath=NULL,
                      outdir="Results",
                      names="love")
  saveRDS(N.sce,file = "./Results/data/Vehicle.rds")
  p1 <- N.sce$v
  p1
  p2 <- DimPlot(N.sce$sce, reduction = "umap",group.by="seurat_clusters",label.size = 4,repel=T,label=T)
  p3 <- DimPlot(N.sce$sce, reduction = "umap",group.by="orig.ident",label.size = 4,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(N.sce$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(N.sce$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(N.sce$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  cl_markers <- N.sce$sce.markers
  top5_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  dh <- DoHeatmap(N.sce$sce, features = top5_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  
  N.sce <- readRDS("./Results/data/Vehicle.rds")
  All.Immune <- "Ptprc" #免疫细胞
  B.cell <- c("Bank1","Cd19","Cd79a","Cd79b","Ms4a1")
  T.cell <- c("Cd3d","Cd3e","Cd8a","Cd4")
  Dendritic.cell <- c("Siglech","Bst2","Itgax","Ccr7","Cst3")
  Macrophage <- c("Cd86","Cd80","Nos2","Cd163","Apoe","Spp1","Ccl4","Lyz2")
  Monocyte <- c("S100a8","S100a9","Cd14")
  NK <- c("Klrk1","Klrb1c","Ncr1")
  plasma <- c("Tnfrsf17","Sdc1")
  EPCAM <- "Cryab" #上皮细胞
  p.immune <- FeaturePlot(N.sce$sce, features = c(All.Immune))
  p.B.cell <- FeaturePlot(N.sce$sce, features = B.cell)
  p.T.cell <- FeaturePlot(N.sce$sce, features = T.cell)
  p.dendritic.cell <- FeaturePlot(N.sce$sce, features = Dendritic.cell)
  p.macrophage <- FeaturePlot(N.sce$sce, features = Macrophage)
  p.monocyte <- FeaturePlot(N.sce$sce, features = Monocyte)
  p.nk <- FeaturePlot(N.sce$sce, features = NK)
  p.plasma <- FeaturePlot(N.sce$sce, features = plasma)
  p.epcam <- FeaturePlot(N.sce$sce, features = EPCAM)
  
  #mouse celltype annotation from guoguoji<Mapping the Mouse Cell Atlas by Microwell-Seq>
  Endothelial_cell <- c("Kdr","Nrp1","Eng","Gpihbp1","Ptprb","Esam","Egfl7","Pecam1")
  Kuppfer_cell <- c("Clec4f","Cd68")
  Erythroblast_Hbb_bs <- c("Hbb-bs","Hba-a2m","Hba-a1")
  T_cell_Trbc2_high <- c("Trbc2","Cd3d","Trbc1","Thy1","Cd3g","Lat","Ccl5")
  Dendritic.cell <- c("Cd74","Irf5","Irf8")
  Dendritic_cell_Siglech <- c("Siglech","Ly6d","Irf8","Bst2")
  Granulocyte.cell <- c("S100a9","S100a8")
  Macrophage_Chil3_high <- c("Ccl9","Cd14","Mpeg1","Cd74","Chil3")
  Hepatocyte_Fabp1_high <- c("Fabp1")
  Pericentral_hepatocytes <- c("Cyp2e1")
  T_cell_Gzma <- c("Gzma","Ccl5","Cd7","Cd3g","Trbc1","Trbc2")
  Epithelia_Spp1_cell <- c("Spp1","Krt8")
  Erythroblast_Hbb_bt_high <- c("Hbb-bt","Hba-a2","Hba-a1","Hbb-bs")
  Periportal_hepatocyte <- c("Apoa1","Apoa2","Apoa5","Apob")
  Epithelial_cell <- c("Krt19","Krt8","Mmp7","Krt18","Muc1","Krt23")
  B_cell_Fcmr <- c("Fcmr","Cd79a","Cd79b","Cd19")
  B_cell_Jchain <- c("Jchain","Mzb1","Cd79b")
  Stromal_cell <- c("Timp2","Col1a2","Col3a1","Col14a1","Fbln2","Vim")
  Hepatocyte_mt_Nd4_high <- c("Echs1","Krt19","Krt8","Mmp7","Krt18","Muc1","Krt23","mt-Nd4")
  Neutrophil_Ngp_high <- c("Ngp","S100a9","Lcn2","S100a8","Cd177")
  
  FeaturePlot(N.sce$sce, features = Endothelial_cell)
  FeaturePlot(N.sce$sce, features = Kuppfer_cell)
  FeaturePlot(N.sce$sce, features = Erythroblast_Hbb_bs)
  FeaturePlot(N.sce$sce, features = T_cell_Trbc2_high)
  FeaturePlot(N.sce$sce, features = Dendritic.cell)
  FeaturePlot(N.sce$sce, features = Dendritic_cell_Siglech)
  FeaturePlot(N.sce$sce, features = T.cell)
  FeaturePlot(N.sce$sce, features = Granulocyte.cell)
  FeaturePlot(N.sce$sce, features = Macrophage_Chil3_high)
  FeaturePlot(N.sce$sce, features = Hepatocyte_Fabp1_high)
  FeaturePlot(N.sce$sce, features = Pericentral_hepatocytes)
  FeaturePlot(N.sce$sce, features = T_cell_Gzma)
  FeaturePlot(N.sce$sce, features = Epithelia_Spp1_cell)
  FeaturePlot(N.sce$sce, features = Erythroblast_Hbb_bt_high)
  FeaturePlot(N.sce$sce, features = Periportal_hepatocyte)
  FeaturePlot(N.sce$sce, features = Epithelial_cell)
  FeaturePlot(N.sce$sce, features = B_cell_Fcmr)
  FeaturePlot(N.sce$sce, features = B_cell_Jchain)
  FeaturePlot(N.sce$sce, features = Stromal_cell)
  FeaturePlot(N.sce$sce, features = Hepatocyte_mt_Nd4_high)
  FeaturePlot(N.sce$sce, features = Neutrophil_Ngp_high)
}
# R848
{
  R.sce <- FastSeurat(count=NULL,
                      dir.name="/Volumes/HDD4/tianjin_supple/tianjin_singleCell/mRNA/R/mRNA/outs/filtered_feature_bc_matrix/",
                      obj=NULL,#min.cells=3,
                      species=c("human","mouse")[2],
                      min.features=0,
                      max.features=300000,
                      percent.mt.num=100,
                      plot=F,
                      pcSelect=30,
                      project="R848",
                      nfeatures=2000,
                      vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                      all.scale = F,
                      npcs = 50,
                      resolution = 0.5,
                      harmony=F,
                      doublet=F,
                      isMarkers=F,
                      perplexity = 30,
                      cellCycle=F,
                      features=NULL,
                      filepath=NULL,
                      outdir="Results",
                      names="love")
  p1 <- R.sce$v
  p1
  p2 <- DimPlot(R.sce$sce, reduction = "umap",group.by="seurat_clusters",label.size = 4,repel=T,label=T)
  p3 <- DimPlot(R.sce$sce, reduction = "umap",group.by="orig.ident",label.size = 4,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(R.sce$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(R.sce$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(R.sce$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  R.sce <- FastSeurat(count=NULL,
                      dir.name="/Volumes/HDD4/tianjin_supple/tianjin_singleCell/mRNA/R/mRNA/outs/filtered_feature_bc_matrix/",
                      obj=NULL,#min.cells=3,
                      species=c("human","mouse")[2],
                      min.features=200,
                      max.features=6000,
                      percent.mt.num=50,
                      plot=F,
                      pcSelect=30,
                      project="R848",
                      nfeatures=2000,
                      vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                      all.scale = F,
                      npcs = 50,
                      resolution = 0.5,
                      harmony=F,
                      doublet=T,
                      isMarkers=T,
                      perplexity = 30,
                      cellCycle=T,
                      features=NULL,
                      filepath=NULL,
                      outdir="Results",
                      names="love")
  saveRDS(R.sce,file = "./Results/data/R848.rds")
  p1 <- R.sce$v
  p1
  p2 <- DimPlot(R.sce$sce, reduction = "umap",group.by="seurat_clusters",label.size = 4,repel=T,label=T)
  p3 <- DimPlot(R.sce$sce, reduction = "umap",group.by="orig.ident",label.size = 4,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(R.sce$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(R.sce$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(R.sce$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  cl_markers <- R.sce$sce.markers
  top5_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  dh <- DoHeatmap(R.sce$sce, features = top5_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  
  #细胞类型注释
  R.sce <- readRDS("./Results/data/R848.rds")
  All.Immune <- "Ptprc" #免疫细胞
  B.cell <- c("Bank1","Cd19","Cd79a","Cd79b","Ms4a1")
  T.cell <- c("Cd3d","Cd3e","Cd8a","Cd4")
  Dendritic.cell <- c("Siglech","Bst2","Itgax","Ccr7","Cst3")
  Macrophage <- c("Cd86","Cd80","Nos2","Cd163","Apoe","Spp1","Ccl4","Lyz2")
  Monocyte <- c("S100a8","S100a9","Cd14")
  NK <- c("Klrk1","Klrb1c","Ncr1")
  plasma <- c("Tnfrsf17","Sdc1")
  EPCAM <- "Cryab" #上皮细胞
  p.immune <- FeaturePlot(R.sce$sce, features = c(All.Immune))
  p.immune
  p.B.cell <- FeaturePlot(R.sce$sce, features = B.cell)
  p.B.cell
  p.T.cell <- FeaturePlot(R.sce$sce, features = T.cell)
  p.T.cell
  p.dendritic.cell <- FeaturePlot(R.sce$sce, features = Dendritic.cell)
  p.dendritic.cell
  p.macrophage <- FeaturePlot(R.sce$sce, features = Macrophage)
  p.macrophage
  p.monocyte <- FeaturePlot(R.sce$sce, features = Monocyte)
  p.monocyte
  p.nk <- FeaturePlot(R.sce$sce, features = NK)
  p.nk
  p.plasma <- FeaturePlot(R.sce$sce, features = plasma)
  p.plasma
  p.epcam <- FeaturePlot(R.sce$sce, features = EPCAM)
  p.epcam
  #mouse celltype annotation from guoguoji<Mapping the Mouse Cell Atlas by Microwell-Seq>
  Endothelial_cell <- c("Kdr","Nrp1","Eng","Gpihbp1","Ptprb","Esam","Egfl7","Pecam1")
  Kuppfer_cell <- c("Clec4f","Cd68")
  Erythroblast_Hbb_bs <- c("Hbb-bs","Hba-a2m","Hba-a1")
  T_cell_Trbc2_high <- c("Trbc2","Cd3d","Trbc1","Thy1","Cd3g","Lat","Ccl5")
  Dendritic.cell <- c("Cd74","Irf5","Irf8")
  Dendritic_cell_Siglech <- c("Siglech","Ly6d","Irf8","Bst2")
  Granulocyte.cell <- c("S100a9","S100a8")
  Macrophage_Chil3_high <- c("Ccl9","Cd14","Mpeg1","Cd74","Chil3")
  Hepatocyte_Fabp1_high <- c("Fabp1")
  Pericentral_hepatocytes <- c("Cyp2e1")
  T_cell_Gzma <- c("Gzma","Ccl5","Cd7","Cd3g","Trbc1","Trbc2")
  Epithelia_Spp1_cell <- c("Spp1","Krt8")
  Erythroblast_Hbb_bt_high <- c("Hbb-bt","Hba-a2","Hba-a1","Hbb-bs")
  Periportal_hepatocyte <- c("Apoa1","Apoa2","Apoa5","Apob")
  Epithelial_cell <- c("Krt19","Krt8","Mmp7","Krt18","Muc1","Krt23")
  B_cell_Fcmr <- c("Fcmr","Cd79a","Cd79b","Cd19")
  B_cell_Jchain <- c("Jchain","Mzb1","Cd79b")
  Stromal_cell <- c("Timp2","Col1a2","Col3a1","Col14a1","Fbln2","Vim")
  Hepatocyte_mt_Nd4_high <- c("Echs1","Krt19","Krt8","Mmp7","Krt18","Muc1","Krt23","mt-Nd4")
  Neutrophil_Ngp_high <- c("Ngp","S100a9","Lcn2","S100a8","Cd177")
  FeaturePlot(R.sce$sce, features = Endothelial_cell)
  FeaturePlot(R.sce$sce, features = Kuppfer_cell)
  FeaturePlot(R.sce$sce, features = Erythroblast_Hbb_bs)
  FeaturePlot(R.sce$sce, features = T_cell_Trbc2_high)
  FeaturePlot(R.sce$sce, features = Dendritic.cell)
  FeaturePlot(R.sce$sce, features = Dendritic_cell_Siglech)
  FeaturePlot(R.sce$sce, features = T.cell)
  FeaturePlot(R.sce$sce, features = Granulocyte.cell)
  FeaturePlot(R.sce$sce, features = Macrophage_Chil3_high)
  FeaturePlot(R.sce$sce, features = Hepatocyte_Fabp1_high)
  FeaturePlot(R.sce$sce, features = Pericentral_hepatocytes)
  FeaturePlot(R.sce$sce, features = T_cell_Gzma)
  FeaturePlot(R.sce$sce, features = Epithelia_Spp1_cell)
  FeaturePlot(R.sce$sce, features = Erythroblast_Hbb_bt_high)
  FeaturePlot(R.sce$sce, features = Periportal_hepatocyte)
  FeaturePlot(R.sce$sce, features = Epithelial_cell)
  FeaturePlot(R.sce$sce, features = B_cell_Fcmr)
  FeaturePlot(R.sce$sce, features = B_cell_Jchain)
  FeaturePlot(R.sce$sce, features = Stromal_cell)
  FeaturePlot(R.sce$sce, features = Hepatocyte_mt_Nd4_high)
  FeaturePlot(R.sce$sce, features = Neutrophil_Ngp_high)
}
# Sofa
{
  S.sce <- FastSeurat(count=NULL,
                      dir.name="/Volumes/HDD4/tianjin_supple/tianjin_singleCell/mRNA/S/mRNA/outs/filtered_feature_bc_matrix/",
                      obj=NULL,#min.cells=3,
                      species=c("human","mouse")[2],
                      min.features=0,
                      max.features=300000,
                      percent.mt.num=100,
                      plot=F,
                      pcSelect=30,
                      project="Sofa",
                      nfeatures=2000,
                      vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                      all.scale = F,
                      npcs = 50,
                      resolution = 0.5,
                      harmony=F,
                      doublet=F,
                      isMarkers=F,
                      perplexity = 30,
                      cellCycle=F,
                      features=NULL,
                      filepath=NULL,
                      outdir="Results",
                      names="love")
  p1 <- S.sce$v
  p1
  p2 <- DimPlot(S.sce$sce, reduction = "umap",group.by="seurat_clusters",label.size = 4,repel=T,label=T)
  p3 <- DimPlot(S.sce$sce, reduction = "umap",group.by="orig.ident",label.size = 4,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(S.sce$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(S.sce$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(S.sce$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  S.sce <- FastSeurat(count=NULL,
                      dir.name="/Volumes/HDD4/tianjin_supple/tianjin_singleCell/mRNA/S/mRNA/outs/filtered_feature_bc_matrix/",
                      obj=NULL,#min.cells=3,
                      species=c("human","mouse")[2],
                      min.features=200,
                      max.features=6000,
                      percent.mt.num=50,
                      plot=F,
                      pcSelect=30,
                      project="Sofa",
                      nfeatures=2000,
                      vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                      all.scale = F,
                      npcs = 50,
                      resolution = 0.5,
                      harmony=F,
                      doublet=T,
                      isMarkers=T,
                      perplexity = 30,
                      cellCycle=T,
                      features=NULL,
                      filepath=NULL,
                      outdir="Results",
                      names="love")
  saveRDS(S.sce,file = "./Results/data/Sofa.rds")
  p1 <- S.sce$v
  p1
  p2 <- DimPlot(S.sce$sce, reduction = "umap",group.by="seurat_clusters",label.size = 4,repel=T,label=T)
  p3 <- DimPlot(S.sce$sce, reduction = "umap",group.by="orig.ident",label.size = 4,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(S.sce$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(S.sce$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(S.sce$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  cl_markers <- S.sce$sce.markers
  top5_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  dh <- DoHeatmap(S.sce$sce, features = top5_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  dh
  #细胞类型注释
  S.sce <- readRDS("./Results/data/Sofa.rds")
  All.Immune <- "Ptprc" #免疫细胞
  B.cell <- c("Bank1","Cd19","Cd79a","Cd79b","Ms4a1")
  T.cell <- c("Cd3d","Cd3e","Cd8a","Cd4")
  Dendritic.cell <- c("Siglech","Bst2","Itgax","Ccr7","Cst3")
  Macrophage <- c("Cd86","Cd80","Nos2","Cd163","Apoe","Spp1","Ccl4","Lyz2")
  Monocyte <- c("S100a8","S100a9","Cd14")
  NK <- c("Klrk1","Klrb1c","Ncr1")
  plasma <- c("Tnfrsf17","Sdc1")
  EPCAM <- "Cryab" #上皮细胞
  p.immune <- FeaturePlot(S.sce$sce, features = c(All.Immune))
  p.immune
  p.B.cell <- FeaturePlot(S.sce$sce, features = B.cell)
  p.B.cell
  p.T.cell <- FeaturePlot(S.sce$sce, features = T.cell)
  p.T.cell
  p.dendritic.cell <- FeaturePlot(S.sce$sce, features = Dendritic.cell)
  p.dendritic.cell
  p.macrophage <- FeaturePlot(S.sce$sce, features = Macrophage)
  p.macrophage
  p.monocyte <- FeaturePlot(S.sce$sce, features = Monocyte)
  p.monocyte
  p.nk <- FeaturePlot(S.sce$sce, features = NK)
  p.nk
  p.plasma <- FeaturePlot(S.sce$sce, features = plasma)
  p.plasma
  p.epcam <- FeaturePlot(S.sce$sce, features = EPCAM)
  p.epcam
  #mouse celltype annotation from guoguoji<Mapping the Mouse Cell Atlas by Microwell-Seq>
  Endothelial_cell <- c("Kdr","Nrp1","Eng","Gpihbp1","Ptprb","Esam","Egfl7","Pecam1")
  Kuppfer_cell <- c("Clec4f","Cd68")
  Erythroblast_Hbb_bs <- c("Hbb-bs","Hba-a2m","Hba-a1")
  T_cell_Trbc2_high <- c("Trbc2","Cd3d","Trbc1","Thy1","Cd3g","Lat","Ccl5")
  Dendritic.cell <- c("Cd74","Irf5","Irf8")
  Dendritic_cell_Siglech <- c("Siglech","Ly6d","Irf8","Bst2")
  Granulocyte.cell <- c("S100a9","S100a8")
  Macrophage_Chil3_high <- c("Ccl9","Cd14","Mpeg1","Cd74","Chil3")
  Hepatocyte_Fabp1_high <- c("Fabp1")
  Pericentral_hepatocytes <- c("Cyp2e1")
  T_cell_Gzma <- c("Gzma","Ccl5","Cd7","Cd3g","Trbc1","Trbc2")
  Epithelia_Spp1_cell <- c("Spp1","Krt8")
  Erythroblast_Hbb_bt_high <- c("Hbb-bt","Hba-a2","Hba-a1","Hbb-bs")
  Periportal_hepatocyte <- c("Apoa1","Apoa2","Apoa5","Apob")
  Epithelial_cell <- c("Krt19","Krt8","Mmp7","Krt18","Muc1","Krt23")
  B_cell_Fcmr <- c("Fcmr","Cd79a","Cd79b","Cd19")
  B_cell_Jchain <- c("Jchain","Mzb1","Cd79b")
  Stromal_cell <- c("Timp2","Col1a2","Col3a1","Col14a1","Fbln2","Vim")
  Hepatocyte_mt_Nd4_high <- c("Echs1","Krt19","Krt8","Mmp7","Krt18","Muc1","Krt23","mt-Nd4")
  Neutrophil_Ngp_high <- c("Ngp","S100a9","Lcn2","S100a8","Cd177")
  FeaturePlot(S.sce$sce, features = Endothelial_cell)
  FeaturePlot(S.sce$sce, features = Kuppfer_cell)
  FeaturePlot(S.sce$sce, features = Erythroblast_Hbb_bs)
  FeaturePlot(S.sce$sce, features = T_cell_Trbc2_high)
  FeaturePlot(S.sce$sce, features = Dendritic.cell)
  FeaturePlot(S.sce$sce, features = Dendritic_cell_Siglech)
  FeaturePlot(S.sce$sce, features = T.cell)
  FeaturePlot(S.sce$sce, features = Granulocyte.cell)
  FeaturePlot(S.sce$sce, features = Macrophage_Chil3_high)
  FeaturePlot(S.sce$sce, features = Hepatocyte_Fabp1_high)
  FeaturePlot(S.sce$sce, features = Pericentral_hepatocytes)
  FeaturePlot(S.sce$sce, features = T_cell_Gzma)
  FeaturePlot(S.sce$sce, features = Epithelia_Spp1_cell)
  FeaturePlot(S.sce$sce, features = Erythroblast_Hbb_bt_high)
  FeaturePlot(S.sce$sce, features = Periportal_hepatocyte)
  FeaturePlot(S.sce$sce, features = Epithelial_cell)
  FeaturePlot(S.sce$sce, features = B_cell_Fcmr)
  FeaturePlot(S.sce$sce, features = B_cell_Jchain)
  FeaturePlot(S.sce$sce, features = Stromal_cell)
  FeaturePlot(S.sce$sce, features = Hepatocyte_mt_Nd4_high)
  FeaturePlot(S.sce$sce, features = Neutrophil_Ngp_high)
}
# SR3
{
  SR3.sce <- FastSeurat(count=NULL,
                      dir.name="/Volumes/HDD4/tianjin_supple/tianjin_singleCell/mRNA/SR3/mRNA/outs/filtered_feature_bc_matrix/",
                      obj=NULL,#min.cells=3,
                      species=c("human","mouse")[2],
                      min.features=0,
                      max.features=300000,
                      percent.mt.num=100,
                      plot=F,
                      pcSelect=30,
                      project="SR3",
                      nfeatures=2000,
                      vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                      all.scale = F,
                      npcs = 50,
                      resolution = 0.5,
                      harmony=F,
                      doublet=F,
                      isMarkers=F,
                      perplexity = 30,
                      cellCycle=F,
                      features=NULL,
                      filepath=NULL,
                      outdir="Results",
                      names="love")
  p1 <- SR3.sce$v
  p1
  p2 <- DimPlot(SR3.sce$sce, reduction = "umap",group.by="seurat_clusters",label.size = 4,repel=T,label=T)
  p3 <- DimPlot(SR3.sce$sce, reduction = "umap",group.by="orig.ident",label.size = 4,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(SR3.sce$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(SR3.sce$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(SR3.sce$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  SR3.sce <- FastSeurat(count=NULL,
                      dir.name="/Volumes/HDD4/tianjin_supple/tianjin_singleCell/mRNA/SR3/mRNA/outs/filtered_feature_bc_matrix/",
                      obj=NULL,#min.cells=3,
                      species=c("human","mouse")[2],
                      min.features=200,
                      max.features=5000,
                      percent.mt.num=50,
                      plot=F,
                      pcSelect=30,
                      project="SR3",
                      nfeatures=2000,
                      vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                      all.scale = F,
                      npcs = 50,
                      resolution = 0.5,
                      harmony=F,
                      doublet=T,
                      isMarkers=T,
                      perplexity = 30,
                      cellCycle=T,
                      features=NULL,
                      filepath=NULL,
                      outdir="Results",
                      names="love")
  saveRDS(SR3.sce,file = "./Results/data/SR3.rds")
  p1 <- SR3.sce$v
  p1
  p2 <- DimPlot(SR3.sce$sce, reduction = "umap",group.by="seurat_clusters",label.size = 4,repel=T,label=T)
  p3 <- DimPlot(SR3.sce$sce, reduction = "umap",group.by="orig.ident",label.size = 4,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(SR3.sce$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(SR3.sce$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(SR3.sce$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  cl_markers <- SR3.sce$sce.markers
  top5_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  dh <- DoHeatmap(SR3.sce$sce, features = top5_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  #细胞类型注释
  SR3.sce <- readRDS("./Results/data/SR3.rds")
  All.Immune <- "Ptprc" #免疫细胞
  B.cell <- c("Bank1","Cd19","Cd79a","Cd79b","Ms4a1")
  T.cell <- c("Cd3d","Cd3e","Cd8a","Cd4")
  Dendritic.cell <- c("Siglech","Bst2","Itgax","Ccr7","Cst3")
  Macrophage <- c("Cd86","Cd80","Nos2","Cd163","Apoe","Spp1","Ccl4","Lyz2")
  Monocyte <- c("S100a8","S100a9","Cd14")
  NK <- c("Klrk1","Klrb1c","Ncr1")
  plasma <- c("Tnfrsf17","Sdc1")
  EPCAM <- "Cryab" #上皮细胞
  p.immune <- FeaturePlot(SR3.sce$sce, features = c(All.Immune))
  p.immune
  p.B.cell <- FeaturePlot(SR3.sce$sce, features = B.cell)
  p.B.cell
  p.T.cell <- FeaturePlot(SR3.sce$sce, features = T.cell)
  p.T.cell
  p.dendritic.cell <- FeaturePlot(SR3.sce$sce, features = Dendritic.cell)
  p.dendritic.cell
  p.macrophage <- FeaturePlot(SR3.sce$sce, features = Macrophage)
  p.macrophage
  p.monocyte <- FeaturePlot(SR3.sce$sce, features = Monocyte)
  p.monocyte
  p.nk <- FeaturePlot(SR3.sce$sce, features = NK)
  p.nk
  p.plasma <- FeaturePlot(SR3.sce$sce, features = plasma)
  p.plasma
  p.epcam <- FeaturePlot(SR3.sce$sce, features = EPCAM)
  p.epcam
  #mouse celltype annotation from guoguoji<Mapping the Mouse Cell Atlas by Microwell-Seq>
  Endothelial_cell <- c("Kdr","Nrp1","Eng","Gpihbp1","Ptprb","Esam","Egfl7","Pecam1")
  Kuppfer_cell <- c("Clec4f","Cd68")
  Erythroblast_Hbb_bs <- c("Hbb-bs","Hba-a2m","Hba-a1")
  T_cell_Trbc2_high <- c("Trbc2","Cd3d","Trbc1","Thy1","Cd3g","Lat","Ccl5")
  Dendritic.cell1 <- c("Cd74","Irf5","Irf8")
  Dendritic_cell_Siglech <- c("Siglech","Ly6d","Irf8","Bst2")
  Granulocyte.cell <- c("S100a9","S100a8")
  Macrophage_Chil3_high <- c("Ccl9","Cd14","Mpeg1","Cd74","Chil3")
  Hepatocyte_Fabp1_high <- c("Fabp1")
  Pericentral_hepatocytes <- c("Cyp2e1")
  T_cell_Gzma <- c("Gzma","Ccl5","Cd7","Cd3g","Trbc1","Trbc2")
  Epithelia_Spp1_cell <- c("Spp1","Krt8")
  Erythroblast_Hbb_bt_high <- c("Hbb-bt","Hba-a2","Hba-a1","Hbb-bs")
  Periportal_hepatocyte <- c("Apoa1","Apoa2","Apoa5","Apob")
  Epithelial_cell <- c("Krt19","Krt8","Mmp7","Krt18","Muc1","Krt23")
  B_cell_Fcmr <- c("Fcmr","Cd79a","Cd79b","Cd19")
  B_cell_Jchain <- c("Jchain","Mzb1","Cd79b")
  Stromal_cell <- c("Timp2","Col1a2","Col3a1","Col14a1","Fbln2","Vim")
  Hepatocyte_mt_Nd4_high <- c("Echs1","Krt19","Krt8","Mmp7","Krt18","Muc1","Krt23","mt-Nd4")
  Neutrophil_Ngp_high <- c("Ngp","S100a9","Lcn2","S100a8","Cd177")
  FeaturePlot(SR3.sce$sce, features = Endothelial_cell)
  FeaturePlot(SR3.sce$sce, features = Kuppfer_cell)
  FeaturePlot(SR3.sce$sce, features = Erythroblast_Hbb_bs)
  FeaturePlot(SR3.sce$sce, features = T_cell_Trbc2_high)
  FeaturePlot(SR3.sce$sce, features = Dendritic.cell1)
  FeaturePlot(SR3.sce$sce, features = Dendritic_cell_Siglech)
  FeaturePlot(SR3.sce$sce, features = T.cell)
  FeaturePlot(SR3.sce$sce, features = Granulocyte.cell)
  FeaturePlot(SR3.sce$sce, features = Macrophage_Chil3_high)
  FeaturePlot(SR3.sce$sce, features = Hepatocyte_Fabp1_high)
  FeaturePlot(SR3.sce$sce, features = Pericentral_hepatocytes)
  FeaturePlot(SR3.sce$sce, features = T_cell_Gzma)
  FeaturePlot(SR3.sce$sce, features = Epithelia_Spp1_cell)
  FeaturePlot(SR3.sce$sce, features = Erythroblast_Hbb_bt_high)
  FeaturePlot(SR3.sce$sce, features = Periportal_hepatocyte)
  FeaturePlot(SR3.sce$sce, features = Epithelial_cell)
  FeaturePlot(SR3.sce$sce, features = B_cell_Fcmr)
  FeaturePlot(SR3.sce$sce, features = B_cell_Jchain)
  FeaturePlot(SR3.sce$sce, features = Stromal_cell)
  FeaturePlot(SR3.sce$sce, features = Hepatocyte_mt_Nd4_high)
  FeaturePlot(SR3.sce$sce, features = Neutrophil_Ngp_high)
}
# merge 4 seurat object
{
  N.sce <- readRDS(file = "./Results/data/Vehicle.rds")
  R848.sce <- readRDS(file = "./Results/data/R848.rds")
  Sofa.sce <- readRDS(file = "./Results/data/Sofa.rds")
  SR3.sce <- readRDS(file = "./Results/data/SR3.rds")
  sce_merge <- FastMergeSeurat(objList=list(N.sce$sce,R848.sce$sce,Sofa.sce$sce,SR3.sce$sce),
                               object.names=c("Vehicle","R848","Sofa","SR3"),
                               species="mouse",
                               project.name.default="ICC_mouse") 
  sce_merge@meta.data <- sce_merge@meta.data[,1:4]
  names <- names(sce_merge@active.ident)
  sce_merge@active.ident <- factor(sce_merge@meta.data[,1],levels = c("Vehicle","R848","Sofa","SR3"))
  names(sce_merge@active.ident) <- names
  sce_merge_res <- FastSeurat(count=NULL,
                              dir.name=NULL,
                              obj=sce_merge,#min.cells=3,
                              species=c("human","mouse")[2],
                              min.features=0,
                              max.features=300000,
                              percent.mt.num=100,
                              plot=F,
                              pcSelect=30,
                              project="ICC_mouse",
                              nfeatures=2000,
                              vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                              all.scale = F,
                              npcs = 50,
                              resolution = 0.5,
                              harmony=F,
                              doublet=F,
                              perplexity = 30,
                              isMarkers=F,
                              cellCycle=F,
                              features=NULL,
                              filepath=NULL,
                              outdir="Results",
                              names="love")
  
  p1 <- sce_merge_res$v
  p1
  
  dittoDimPlot(sce_merge_res$sce, var = "seurat_clusters",reduction.use="umap",do.ellipse=T,do.label = T)
  dittoDimPlot(sce_merge_res$sce, var = "orig.ident", reduction.use="umap",do.ellipse=T,do.label = T)
  p2 <- DimPlot(sce_merge_res$sce, reduction = "umap",group.by="seurat_clusters",label.size = 4,repel=T,label=T)
  p3 <- DimPlot(sce_merge_res$sce, reduction = "umap",group.by="orig.ident",label.size = 4,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(sce_merge_res$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_merge_res$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_merge_res$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  sce_merge_res <- FastSeurat(count=NULL,
                        dir.name=NULL,
                        obj=sce_merge,#min.cells=3,
                        species=c("human","mouse")[2],
                        min.features=200,
                        max.features=6000,
                        percent.mt.num=50,
                        plot=F,
                        pcSelect=30,
                        project="ICC_mouse",
                        nfeatures=2000,
                        vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                        all.scale = F,
                        npcs = 50,
                        resolution = 0.5,
                        harmony=F,
                        doublet=F,
                        isMarkers=T,
                        perplexity = 30,
                        cellCycle=T,
                        features=NULL,
                        filepath=NULL,
                        outdir="Results",
                        names="love")
  saveRDS(sce_merge_res,file = "./Results/data/sce_merge_res.rds")
  sce_merge_res <- readRDS(file = "./Results/data/sce_merge_res.rds")
  bardata <- sce_merge_res$sce@meta.data %>% 
    dplyr::group_by(orig.ident) %>% 
    dplyr::summarise(num=dplyr::n())
  p0 <- ggplot(data=bardata,mapping=aes(x=orig.ident,y=num,fill=orig.ident)) +
          geom_bar(stat = "identity") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  p1 <- sce_merge_res$v
  p1
  p2 <- DimPlot(sce_merge_res$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
  p3 <- DimPlot(sce_merge_res$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(sce_merge_res$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_merge_res$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_merge_res$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  cl_markers <- sce_merge_res$sce.markers
  top20_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  top20_cl_markers %>% filter(.,cluster==10) -> top20_cl10_markers
  top20_cl_markers %>% filter(.,cluster==5) -> top20_cl5_markers
  top20_cl_markers %>% filter(.,cluster==15) -> top20_cl15_markers
  dh <- DoHeatmap(sce_merge_res$sce, features = top20_cl10_markers$gene,slot = "data") + NoLegend()
  dh <- DoHeatmap(sce_merge_res$sce, features = top20_cl5_markers$gene,slot = "data") + NoLegend()
  dh <- DoHeatmap(sce_merge_res$sce, features = top20_cl15_markers$gene,slot = "data") + NoLegend()
  print(dh)
  
  #细胞类型注释
  sce_merge_res <- readRDS("./Results/data/sce_merge_res1.rds")
  All.Immune <- "Ptprc" #免疫细胞
  B.cell <- c("Bank1","Cd19","Cd79a","Cd79b","Ms4a1")
  T.cell <- c("Cd3d","Cd3e","Cd8a","Cd4")
  Dendritic.cell1 <- c("Siglech","Bst2","Itgax","Ccr7","Cst3")
  Macrophage <- c("Cd86","Cd80","Nos2","Cd163","Apoe","Spp1","Ccl4","Lyz2")
  Monocyte <- c("S100a8","S100a9","Cd14")
  NK <- c("Klrk1","Klrb1c","Ncr1")
  plasma <- c("Tnfrsf17","Sdc1","Iglc2","Ighg1","Igkc")
  EPCAM <- "Cryab" #上皮细胞
  p.immune <- FeaturePlot(sce_merge_res$sce, features = c(All.Immune))
  p.immune
  p.B.cell <- FeaturePlot(sce_merge_res$sce, features = B.cell)
  p.B.cell
  p.T.cell <- FeaturePlot(sce_merge_res$sce, features = T.cell)
  p.T.cell
  p.dendritic.cell <- FeaturePlot(sce_merge_res$sce, features = Dendritic.cell1)
  p.dendritic.cell
  p.macrophage <- FeaturePlot(sce_merge_res$sce, features = Macrophage)
  p.macrophage
  p.monocyte <- FeaturePlot(sce_merge_res$sce, features = Monocyte)
  p.monocyte
  p.nk <- FeaturePlot(sce_merge_res$sce, features = NK)
  p.nk
  p.plasma <- FeaturePlot(sce_merge_res$sce, features = plasma)
  p.plasma
  p.epcam <- FeaturePlot(sce_merge_res$sce, features = EPCAM)
  p.epcam
  #mouse celltype annotation from guoguoji<Mapping the Mouse Cell Atlas by Microwell-Seq>
  Endothelial_cell <- c("Kdr","Nrp1","Eng","Gpihbp1","Ptprb","Esam","Egfl7","Pecam1","Vwf")
  Kuppfer_cell <- c("Clec4f","Cd68")
  Erythroblast_Hbb_bs <- c("Hbb-bs","Hba-a2m","Hba-a1")
  T_cell_Trbc2_high <- c("Trbc2","Cd3d","Trbc1","Thy1","Cd3g","Lat","Ccl5")
  Dendritic.cell <- c("Cd74","Irf5","Irf8")
  Dendritic_cell_Siglech <- c("Siglech","Ly6d","Irf8","Bst2")
  Granulocyte.cell <- c("S100a9","S100a8")
  Macrophage_Chil3_high <- c("Ccl9","Cd14","Mpeg1","Cd74","Chil3")
  Hepatocyte_Fabp1_high <- c("Fabp1")
  Pericentral_hepatocytes <- c("Cyp2e1")
  T_cell_Gzma <- c("Gzma","Ccl5","Cd7","Cd3g","Trbc1","Trbc2")
  Epithelia_Spp1_cell <- c("Spp1","Krt8")
  Erythroblast_Hbb_bt_high <- c("Hbb-bt","Hba-a2","Hba-a1","Hbb-bs")
  Periportal_hepatocyte <- c("Apoa1","Apoa2","Apoa5","Apob")
  Epithelial_cell <- c("Krt19","Krt8","Mmp7","Krt18","Muc1","Krt23")
  B_cell_Fcmr <- c("Fcmr","Cd79a","Cd79b","Cd19")
  B_cell_Jchain <- c("Jchain","Mzb1","Cd79b")
  Stromal_cell <- c("Timp2","Col1a2","Col3a1","Col14a1","Fbln2","Vim")
  Hepatocyte_mt_Nd4_high <- c("Echs1","Krt19","Krt8","Mmp7","Krt18","Muc1","Krt23","mt-Nd4")
  Neutrophil_Ngp_high <- c("Ngp","S100a9","Lcn2","S100a8","Cd177")
  FeaturePlot(sce_merge_res$sce, features = Endothelial_cell)
  FeaturePlot(sce_merge_res$sce, features = Kuppfer_cell)
  FeaturePlot(sce_merge_res$sce, features = Erythroblast_Hbb_bs)
  FeaturePlot(sce_merge_res$sce, features = T_cell_Trbc2_high)
  FeaturePlot(sce_merge_res$sce, features = Dendritic.cell)
  FeaturePlot(sce_merge_res$sce, features = Dendritic_cell_Siglech)
  FeaturePlot(sce_merge_res$sce, features = T.cell)
  FeaturePlot(sce_merge_res$sce, features = Granulocyte.cell)
  FeaturePlot(sce_merge_res$sce, features = Macrophage_Chil3_high)
  FeaturePlot(sce_merge_res$sce, features = Hepatocyte_Fabp1_high)
  FeaturePlot(sce_merge_res$sce, features = Pericentral_hepatocytes)
  FeaturePlot(sce_merge_res$sce, features = T_cell_Gzma)
  FeaturePlot(sce_merge_res$sce, features = Epithelia_Spp1_cell)
  FeaturePlot(sce_merge_res$sce, features = Erythroblast_Hbb_bt_high)
  FeaturePlot(sce_merge_res$sce, features = Periportal_hepatocyte)
  FeaturePlot(sce_merge_res$sce, features = Epithelial_cell)
  FeaturePlot(sce_merge_res$sce, features = B_cell_Fcmr)
  FeaturePlot(sce_merge_res$sce, features = B_cell_Jchain)
  FeaturePlot(sce_merge_res$sce, features = Stromal_cell)
  FeaturePlot(sce_merge_res$sce, features = Hepatocyte_mt_Nd4_high)
  FeaturePlot(sce_merge_res$sce, features = Neutrophil_Ngp_high)
  
  FeaturePlot(sce_merge_res$sce, features = c("Sdc1","Clu","Ehd4","Cypla1","Plac8","Col8a1","Lrg1","Mgp","Mecom","Cd200"))
  gzh_Granulocyte <- c("Ly6c1","Ly6g","Ngp","Camp")
  gzh_DC <- c("Cd1c","Cd11c", "Itgax","Cd80","Cd83","Siglech","Flt3","Cd209a","Nrp1","Il3ra")
  gzh_Monocytes <- c("F13a1","Cd14","Fcgr3","Ly6c1","Itgam","Ly6g","F10","Lilrb4a")
  gzh_Mφ <- c("Cd14","Cd68","Itgam","Adgre1","Fcgr1","Mpeg1","Grm1","Marco","Fcer1a","Ear2","Cxcl2","Pf4","Csf1r","Trem2","Cd9")
  gzh_NK <- c("Ncam1","Cd7","Il2rb","Cd56","Cd49b","Nkp46","Nk1.4","Nk1.2","Zbtb32","Jun","Ccl3","Fos","Ccl4","Nlfkbia","Fosb","Klf2","Pp1r15a","Pim1","Ler5","Ccl5","Nkg7","Ncr1","Gzma","Klrb1c")
  gzh_Mast <- c("Kit","Pf4","Mcpt4","Cd32","Cd33","Cd117","Fcer1a","Integrinb7","Furin","Il1rl1")
  FeaturePlot(sce_merge_res$sce, features = gzh_Granulocyte)
  FeaturePlot(sce_merge_res1$sce, features = gzh_DC)
  FeaturePlot(sce_merge_res1$sce, features = gzh_Monocytes)
  FeaturePlot(sce_merge_res1$sce, features = gzh_Mφ)
  FeaturePlot(sce_merge_res$sce, features = gzh_NK)
  FeaturePlot(sce_merge_res$sce, features = gzh_Mast)
  FeaturePlot(sce_merge_res$sce, features = genes_hep)
  FeaturePlot(sce_merge_res$sce, features = genes_endo)
  FeaturePlot(sce_merge_res$sce, features = genes_kuppfer)
  FeaturePlot(sce_merge_res$sce, features = genes_nk)
  FeaturePlot(sce_merge_res$sce, features = c("Plet1","Wdfy4"))
  
  
  #过滤双细胞 cluster13
  sce_merge_res$sce <- subset(sce_merge_res$sce,seurat_clusters != 13)
  sce_merge_res1 <- FastSeurat(count=NULL,
                              dir.name=NULL,
                              obj=sce_merge_res$sce,#min.cells=3,
                              species=c("human","mouse")[2],
                              min.features=200,
                              max.features=6000,
                              percent.mt.num=50,
                              plot=F,
                              pcSelect=30,
                              project="ICC_mouse",
                              nfeatures=2000,
                              vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                              all.scale = F,
                              npcs = 50,
                              resolution = 0.5,
                              harmony=F,
                              doublet=F,
                              isMarkers=T,
                              perplexity = 30,
                              cellCycle=T,
                              features=NULL,
                              filepath=NULL,
                              outdir="Results",
                              names="love")
  saveRDS(sce_merge_res1,"./Results/data/sce_merge_res1.rds")
  sce_merge_res1 <- readRDS("./Results/data/sce_merge_res1.rds")
  p1 <- sce_merge_res1$v
  p1
  bardata <- sce_merge_res1$sce@meta.data %>% 
    group_by(orig.ident) %>% 
    summarise(num=n())
  p0 <- ggplot(data=bardata,mapping=aes(x=orig.ident,y=num,fill=orig.ident)) +
    geom_bar(stat = "identity") +
    #geom_text_repel() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  
  p2 <- DimPlot(sce_merge_res1$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
  p3 <- DimPlot(sce_merge_res1$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(sce_merge_res1$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_merge_res1$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_merge_res1$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  head(sce_merge_res1$sce@meta.data)
  
  
  Epithelial_cell <- c("Krt8","Mmp7","Krt18")
  Stromal_cell <- c("Col1a2","Col3a1","Col14a1","Fbln2")
  All.Immune <- "Ptprc" #免疫细胞
  B.cell <- c("Cd19","Cd79a","Cd79b","Ms4a1")
  T.cell <- c("Cd3d","Cd3e","Cd8a","Cd4")
  NK <- c("Gzma","Klrb1c","Ncr1")
  Mast <- c("Kit","Mcpt4","Fcer1a")
  Neutrophil <- c("S100a9","Lcn2","S100a8")
  Dendritic_cell <- c("Ccr7","Plet1","Wdfy4")
  Marcrophage <- c("Apoe","Lyz2")
  tail(sce_merge_res1$sce@meta.data)
  FeaturePlot(sce_merge_res1$sce, features = Epithelial_cell)
  FeaturePlot(sce_merge_res1$sce, features = Stromal_cell)
  FeaturePlot(sce_merge_res1$sce, features = All.Immune)
  FeaturePlot(sce_merge_res1$sce, features = B.cell)
  FeaturePlot(sce_merge_res1$sce, features = T.cell)
  FeaturePlot(sce_merge_res1$sce, features = NK)
  FeaturePlot(sce_merge_res1$sce, features = Mast)
  FeaturePlot(sce_merge_res1$sce, features = Neutrophil)
  FeaturePlot(sce_merge_res1$sce, features = Dendritic_cell)
  FeaturePlot(sce_merge_res1$sce, features = Marcrophage)
  cl_markers <- sce_merge_res1$sce.markers
  top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  dh <- DoHeatmap(sce_merge_res1$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  top20_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  top20_cl_markers %>% filter(.,cluster==10) -> top20_cl10_markers
  top20_cl_markers %>% filter(.,cluster==4) -> top20_cl4_markers
  top20_cl_markers %>% filter(.,cluster==14) -> top20_cl14_markers
  dh <- DoHeatmap(sce_merge_res1$sce, features = top20_cl10_markers$gene,slot = "data") + NoLegend()
  dh <- DoHeatmap(sce_merge_res1$sce, features = top20_cl5_markers$gene,slot = "data") + NoLegend()
  dh <- DoHeatmap(sce_merge_res1$sce, features = top20_cl15_markers$gene,slot = "data") + NoLegend()
  print(dh)
  head(sce_merge_res1$sce@meta.data)
  CD8_T <- subset(sce_merge_res1$sce,seurat_clusters %in% c(1,6,12))
  sce_merge_res1$sce@meta.data$celltype <- mapvalues(sce_merge_res1$sce@meta.data$seurat_clusters,from=0:16,to=c("Epithelial","CD8_T","Neutrophil","Neutrophil","Marcrophage","Epithelial","CD8_T","NK","CD4_T","Stromal","Dendritic_cell","Epithelial","CD8_T","B_cell","Dendritic_cell","Plasma","Mast"))
  sample_cell_proportion <- sce_merge_res1$sce@meta.data %>% dplyr::group_by(orig.ident)   %>% dplyr::summarise(proportion=as.numeric(table(celltype)/length(celltype)),names=str_replace_all(names(table(celltype))," ","_")) 
  sample_cell_proportion <- sce_merge_res1$sce@meta.data %>% filter(.,!celltype %in% c("Epithelial","Stromal")) %>% dplyr::group_by(orig.ident)   %>% dplyr::summarise(proportion=as.numeric(table(celltype)/length(celltype)),names=str_replace_all(names(table(celltype))," ","_")) 
  sample_cell_proportion <- sample_cell_proportion %>% filter(.,!names %in% c("Epithelial","Stromal"))
  library(RColorBrewer)
  mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
  ##柱状图整体比例情况
  ggplot(sample_cell_proportion,aes(x=orig.ident,y=proportion,fill=names))+
    geom_bar(stat="identity",position="fill")+
    ggtitle("Proportion of immune cell subtypes")+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size=20)) +
    guides(fill=guide_legend(title="sample cell subtypes"))+
    scale_fill_manual(values = mypalette(11))
  sample_cell_proportion %>% dcast(.,names~orig.ident,value.var="proportion") %>% select(.,2:5) 
  chisq.test(table(sce_merge_res1$sce@meta.data$orig.ident,sce_merge_res1$sce@meta.data$celltype))
  pos <- which(!sce_merge_res1$sce@meta.data$celltype %in% c("Epithelial","Stromal"))
  chisq.test(table(sce_merge_res1$sce@meta.data$orig.ident[pos],sce_merge_res1$sce@meta.data$celltype[pos])[,-c(1,7)])
  
  #提取巨噬细胞,进行精细划分
  Marcrophage_sce <- subset(sce_merge_res1$sce,seurat_clusters %in% 4)
  Marcrophage1 <- FastSeurat(count=NULL,
                            dir.name=NULL,
                            obj=Marcrophage_sce,#min.cells=3,
                            species=c("human","mouse")[2],
                            min.features=200,
                            max.features=6000,
                            percent.mt.num=50,
                            plot=F,
                            pcSelect=30,
                            project="ICC_mouse",
                            nfeatures=2000,
                            vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                            all.scale = F,
                            npcs = 50,
                            resolution = 0.5,
                            harmony=F,
                            doublet=F,
                            isMarkers=T,
                            perplexity = 30,
                            cellCycle=F,
                            features=NULL,
                            filepath=NULL,
                            outdir="Results",
                            names="love")
  p2 <- DimPlot(Marcrophage1$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
  p3 <- DimPlot(Marcrophage1$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T)
  p2 | p3
  cl_markers <- Marcrophage1$sce.markers
  top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  dh <- DoHeatmap(Marcrophage1$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  write.table(cl_markers,file="./Results/markers.txt",sep="\t",quote=F,row.names=F,col.names=T)

  ## 计算巨噬细胞marker得分，查看我们的巨噬细胞倾向于什么
  eset <- Marcrophage_sce@assays$RNA@data
  M1_down <- format_msigdb("/Users/biofly/project/shijian/data/GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN.gmt", ont = "term", gene = "gene")
  M1_up <- format_msigdb("/Users/biofly/project/shijian/data/GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP.gmt", ont = "term", gene = "gene")
  library(stringr)
  M1_down$GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN <- str_to_title(unlist(M1_down))
  M1_up$GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP <- str_to_title(unlist(M1_up))
  signature <- c(M1_down,M1_up)
  Marcrophage_zscore <- calculate_sig_score_zscore(pdata = NULL,
                                eset=eset,
                                signature=signature,
                                mini_gene_count=2,
                                column_of_sample="ID",
                                adjust_eset = FALSE)
  Marcrophage_ssgsea <- calculate_sig_score_ssgsea(pdata = NULL,
                             eset=eset,
                             signature=signature,
                             mini_gene_count=2,
                             column_of_sample="ID",
                             adjust_eset = FALSE)
  Marcrophage1$sce@meta.data %>%tibble::rownames_to_column("cell") %>% left_join(.,Marcrophage_zscore,by=c("cell"="ID")) -> meta.data
  rownames(meta.data) <- meta.data$cell
  Marcrophage1$sce@meta.data <- meta.data
  saveRDS(Marcrophage1,file="./Results/data/Marcrophage1.rds")
  Marcrophage1 <- readRDS("./Results/data/Marcrophage1.rds")
  melt(Marcrophage_zscore,id.vars=c("Index","ID")) %>% mutate(sample=str_split(ID,"_",simplify = T)[,1]) %>% 
    ggplot(data=.,mapping = aes(x=sample,y=value,fill=sample)) +
    geom_boxplot(outlier.shape = NA) +  # 不显示离群点
    geom_jitter(size=.5) + 
    facet_wrap(vars(variable),nrow = 3, scales = "free_y")
  melt(Marcrophage_ssgsea,id.vars=c("Index","ID")) %>% 
    ggplot(data=.,mapping = aes(x=variable,y=value,fill=variable)) +
    geom_boxplot(outlier.shape = NA) +  # 不显示离群点
    geom_jitter(size=.5)
  Marcrophage1$sce@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% left_join(.,Marcrophage_zscore,by=c("cell"="ID")) %>% mutate(type=ifelse(GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP > GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN,"M1","M2")) %>% 
    ggplot(data=.,mapping = aes(x=UMAP_1,y=UMAP_2,color=GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(type = "viridis")
  Marcrophage1$sce@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% left_join(.,Marcrophage_zscore,by=c("cell"="ID")) %>% mutate(type=ifelse(GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP > GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN,"M1","M2")) %>% 
    ggplot(data=.,mapping = aes(x=UMAP_1,y=UMAP_2,color=GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP)) +
    geom_point() +
    theme_bw() +
    scale_color_continuous(type = "viridis")
  ggplot(data=Marcrophage1$sce@meta.data,mapping = aes(x=orig.ident,y=GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN,color=orig.ident)) +
    geom_boxplot() +
    theme_minimal_grid()
  
  
}

{
  #四组小鼠所有细胞
  {
    dir.create("./code/tianjin_Supple/singlecell_mouse/pictures")
    setwd("./code/tianjin_Supple/singlecell_mouse/pictures")
    sce_merge_res1 <- readRDS("./Results/data/mouse_sce_merge_res1.rds")
    #将cluster5去掉，因为是低活性的肿瘤
    sce_merge_res1$sce <- subset(sce_merge_res1$sce, seurat_clusters!=5)
    #依据样本着色
    p1 <- dittoDimPlot(sce_merge_res1$sce, reduction.use  = "umap", var="orig.ident", size = 1, do.label = TRUE, legend.show = T, main = "4-samples",
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 12))
    p1
    p2 <- dittoDimPlot(sce_merge_res1$sce, reduction.use  = "umap", var="seurat_clusters", size = 1, do.label = TRUE, legend.show = T, main = "4-samples",
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 12))
    p2
    p1 | p2
    head(sce_merge_res1$sce@meta.data)
    p4 <- VlnPlot(sce_merge_res1$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
      geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
    p5 <- VlnPlot(sce_merge_res1$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
      geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
    p6 <- VlnPlot(sce_merge_res1$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
      geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
      theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
    CombinePlots(plots = list(p4,p5,p6),nrow=3)
    Epithelial_cell <- c("Krt8","Mmp7","Krt18")
    Stromal_cell <- c("Col1a2","Col3a1","Col14a1","Fbln2")
    All.Immune <- "Ptprc" #免疫细胞
    B.cell <- c("Cd19","Cd79a","Cd79b","Ms4a1")
    T.cell <- c("Cd3d","Cd3e","Cd8a","Cd4")
    NK <- c("Gzma","Klrb1c","Ncr1")
    Mast <- c("Kit","Mcpt4","Fcer1a")
    Neutrophil <- c("S100a9","Lcn2","S100a8")
    Dendritic_cell <- c("Ccr7","Plet1","Wdfy4")
    Marcrophage <- c("Apoe","Lyz2")
    FeaturePlot_scCustom(sce_merge_res1$sce, features = Epithelial_cell, order = F)
    FeaturePlot_scCustom(sce_merge_res1$sce, features = Stromal_cell, order = F)
    FeaturePlot_scCustom(sce_merge_res1$sce, features = All.Immune, order = F)
    FeaturePlot_scCustom(sce_merge_res1$sce, features = B.cell, order = F)
    FeaturePlot_scCustom(sce_merge_res1$sce, features = T.cell, order = F)
    FeaturePlot_scCustom(sce_merge_res1$sce, features = NK, order = F)
    FeaturePlot_scCustom(sce_merge_res1$sce, features = Mast, order = F)
    FeaturePlot_scCustom(sce_merge_res1$sce, features = Neutrophil, order = F)
    FeaturePlot_scCustom(sce_merge_res1$sce, features = Dendritic_cell, order = F)
    FeaturePlot_scCustom(sce_merge_res1$sce, features = Marcrophage, order = F)
    #sce_merge_res1$sce <- subset(sce_merge_res1$sce,seurat_clusters != 15) #去除cluster15
    #sce_merge_res1$sce@meta.data$seurat_clusters <- factor(sce_merge_res1$sce@meta.data$seurat_clusters,levels = 0:14)
    # 0,11  Epithelial_cell
    # 9 Stromal_cell
    # 13 B_cell
    # 1,12,6 CD8+T_cell
    # 8 CD4+T_cell
    # 7 NK
    # 16 Mast_cell
    # 2,3 Neutrophil_cell
    # 10,14  Dendritic_cell
    # 4 Marcrophage
    # 15 plasma
    #细胞打标签
    sce_merge_res1$sce@meta.data$celltype <- as.character(plyr::mapvalues(sce_merge_res1$sce@meta.data$seurat_clusters,from = c(0:4,6:16),
                                                                          to=c("Epithelial_cell","CD8+T_cell","Neutrophil_cell","Neutrophil_cell","Marcrophage",
                                                                               "CD8+T_cell","NK","CD4+T_cell","Stromal_cell","Dendritic_cell","Epithelial_cell",
                                                                               "CD8+T_cell","B_cell","Dendritic_cell","Plasma","Mast_cell")))
    sce_merge_res1$sce@meta.data$celltype1 <- as.character(plyr::mapvalues(sce_merge_res1$sce@meta.data$seurat_clusters,from = c(0:4,6:16),
                                                                           to=c("Epithelial_cell","T_cell","Neutrophil_cell","Neutrophil_cell","Marcrophage",
                                                                                "T_cell","NK","T_cell","Stromal_cell","Dendritic_cell","Epithelial_cell",
                                                                                "T_cell","B_cell","Dendritic_cell","Plasma","Mast_cell")))
    saveRDS(sce_merge_res1,file = "./Results/data/mouse_sce_merge_res1.rds")
    feature_genes <- c("Krt8","Col3a1","Cd79a","Cd3d","Cd8a","Cd4","Gzma","Mcpt4","S100a8","Ccr7","Plet1","Lyz2")
    p <- FeaturePlot_scCustom(sce_merge_res1$sce, features = feature_genes, order = F)
    ggsave(p,filename = "./code/tianjin_Supple/singlecell_mouse/pictures/all_features_plot.pdf",width=17,height = 14)
    p2 <- dittoDimPlot(sce_merge_res1$sce, reduction.use  = "umap", var="celltype", size = 1, do.label = TRUE, 
                       do.ellipse = T, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500, legend.show = T,main = "4-samples") + theme(legend.text = element_text(face="bold",size = 12))
    p2
    ggsave(p2,filename = "./code/tianjin_Supple/singlecell_mouse/pictures/4_samples_cell_annotation.pdf",width=8,height=8)
    p3 <- dittoDimPlot(sce_merge_res1$sce, reduction.use  = "umap", var="celltype", size = 1, do.label = TRUE, cells.use = grepl("Vehicle",sce_merge_res1$sce@meta.data$orig.ident),legend.show = F, main = "Vehicle",
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 12))
    p3
    p4 <- dittoDimPlot(sce_merge_res1$sce, reduction.use  = "umap", var="celltype", size = 1, do.label = TRUE, cells.use = grepl("Sofa",sce_merge_res1$sce@meta.data$orig.ident),legend.show = F, main = "Sofa",
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 12))
    p4
    p5 <- dittoDimPlot(sce_merge_res1$sce, reduction.use  = "umap", var="celltype", size = 1, do.label = TRUE, cells.use = grepl("R848",sce_merge_res1$sce@meta.data$orig.ident),legend.show = F, main = "R848",
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 12))
    p5
    p6 <- dittoDimPlot(sce_merge_res1$sce, reduction.use  = "umap", var="celltype", size = 1, do.label = TRUE, cells.use = grepl("SR3",sce_merge_res1$sce@meta.data$orig.ident),legend.show = T, main = "SR3",
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 12))
    p6
    p <- p3 | p4 | p5 | p6
    ggsave(p,filename = "4_sample_umap.pdf",width = 16,height=5)
    head(sce_merge_res1$sce@meta.data)
    sce_merge_res1$sce@meta.data %>% dplyr::group_by(orig.ident) %>% 
                                     dplyr::summarise(proportion=as.numeric(table(celltype)/length(celltype)),names=str_replace_all(names(table(celltype))," ","_")) %>% 
                                     ungroup() %>% as.data.frame() -> plot.data
    tt.df <- sce_merge_res1$sce@meta.data %>% 
                    dplyr::group_by(orig.ident) %>% 
                    dplyr::summarise(proportion=as.numeric(table(celltype)/length(celltype)),names=str_replace_all(names(table(celltype))," ","_"))  
    tt.df <- tt.df %>% 
                 filter(.,!grepl("Stromal|Epithelial",names)) %>% 
                 group_by(orig.ident) %>% 
                 dplyr::summarise(proportion=sum(proportion)) %>% 
                 bind_cols(names=rep("immune_cell",4)) %>% 
                 bind_rows(plot.data) 
    # sce_merge_res1$sce@meta.data %>% dplyr::group_by(orig.ident) %>% dplyr::summarise(proportion=as.numeric(table(celltype1)/length(celltype1)),names=str_replace_all(names(table(celltype1))," ","_")) %>% 
    #   filter(.,names == "T_cell" | names == "Neutrophil_cell") %>% ungroup() %>% as.data.frame() -> plot.data1
    tt.df$orig.ident <- factor(tt.df$orig.ident,levels = c("Vehicle","Sofa","R848","SR3"))
    tt.df$names <- factor(tt.df$names,levels = c("immune_cell","CD4+T_cell","CD8+T_cell","Neutrophil_cell"))
    p <- ggplot(data = tt.df,mapping = aes(x=orig.ident,y=proportion,fill=orig.ident)) + 
      geom_bar(stat="identity") +
      facet_wrap(vars(names), scales = "free_y") +
      scale_fill_manual(values = dittoColors()) +
      theme_bw()
    ggsave(p,filename = "cellproportion.pdf",width = 16,height=12)
    #
    #barplot的检验不是多组进行，而是相比于对照组，是否显著高的单侧检验；
    #
    chisq.test(table(sce_merge_res1$sce@meta.data$orig.ident,sce_merge_res1$sce@meta.data$celltype))
    ##免疫细胞比例的差异
    head(sce_merge_res1$sce@meta.data)
    table(sce_merge_res1$sce@meta.data$celltype)
    meta.data <- sce_merge_res1$sce@meta.data %>% 
                     mutate(celltype3=if_else(celltype != "Epithelial_cell" & celltype != "Stromal_cell","immune_cell","others"))
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="Sofa")
    table(aa$orig.ident,aa$celltype3)
    fisher.test(table(aa$orig.ident,aa$celltype3),alternative = "greater")
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="R848")
    table(aa$orig.ident,aa$celltype3)
    fisher.test(table(aa$orig.ident,aa$celltype3),alternative = "greater")
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype3)
    fisher.test(table(aa$orig.ident,aa$celltype3),alternative = "greater")
    aa <- meta.data %>% filter(., orig.ident=="Sofa" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype3)
    fisher.test(table(aa$orig.ident,aa$celltype3),alternative = "less")
    aa <- meta.data %>% filter(., orig.ident=="R848" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype3)
    fisher.test(table(aa$orig.ident,aa$celltype3),alternative = "greater")
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype3)
    fisher.test(table(aa$orig.ident,aa$celltype3),alternative = "greater")
    
    meta.data <- sce_merge_res1$sce@meta.data %>% mutate(celltype4=if_else(celltype == "CD4+T_cell","CD4+T_cell","others"))
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="Sofa")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "less")
    
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="R848")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "less")
    
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "less")
    
    meta.data <- sce_merge_res1$sce@meta.data %>% mutate(celltype4=if_else(celltype == "CD8+T_cell","CD8+T_cell","others"))
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="Sofa")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "greater")
    
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="R848")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "greater")
    
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "greater")
    
    meta.data <- sce_merge_res1$sce@meta.data %>% mutate(celltype4=if_else(celltype == "Neutrophil_cell","Neutrophil_cell","others"))
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="Sofa")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "greater")
    
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="R848")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "greater")
    
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4))
    
    aa <- meta.data %>% filter(., orig.ident=="Sofa" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4))
    
    meta.data <- sce_merge_res1$sce@meta.data %>% mutate(celltype4=if_else(celltype == "Marcrophage","Marcrophage","others"))
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="Sofa")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "greater")
    
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="R848")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "greater")
    
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "greater")
    
    aa <- meta.data %>% filter(., orig.ident=="Sofa" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4))
    
    aa <- meta.data %>% filter(., orig.ident=="R848" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4))
    
    meta.data <- sce_merge_res1$sce@meta.data %>% mutate(celltype4=if_else(celltype == "Epithelial_cell","Epithelial_cell","others"))
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="Sofa")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "greater")
    
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="R848")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "greater")
    
    aa <- meta.data %>% filter(., orig.ident=="Vehicle" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4))
    
    aa <- meta.data %>% filter(., orig.ident=="Sofa" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4),alternative = "greater")
    
    aa <- meta.data %>% filter(., orig.ident=="R848" | orig.ident=="SR3")
    table(aa$orig.ident,aa$celltype4)
    fisher.test(table(aa$orig.ident,aa$celltype4))
    # ggplot(data = plot.data1,mapping = aes(x=orig.ident,y=proportion,fill=orig.ident)) + 
    #   geom_bar(stat="identity") +
    #   #geom_errorbar() +
    #   facet_wrap(vars(names), scales = "free_y") +
    #   scale_fill_manual(values = dittoColors()) +
    #   theme_bw()
    #TNFa，IFNr，Gzmb，Gzma, Tim3
    FeaturePlot_scCustom(sce_merge_res1$sce, features = c("Tnf","Ifng","Gzmb","Havcr2"), order = F)
    FeaturePlot_scCustom(sce_merge_res1$sce, features = "nFeature_RNA", order = F)
    features <- c("Tnf","Ifng","Gzmb","Gzma","Havcr2")
    VlnPlot(sce_merge_res1$sce, features = features, group.by="celltype",split.by="orig.ident",stack=T,raster = T,flip=T,sort = "decreasing")
    Stacked_VlnPlot(sce_merge_res1$sce, features = features, group.by="celltype",split.by="orig.ident",x_lab_rotate=T)
    sce_merge_res1$sce@assays$RNA@data[1:4,1:4]
    t(sce_merge_res1$sce@assays$RNA@data[c("Tnf","Ifng","Gzmb","Gzma","Havcr2"),])
    metadata <- sce_merge_res1$sce@meta.data %>% tibble::rownames_to_column("cellID") %>% 
                  dplyr::select(.,c(1,2,11,12)) %>% 
                  bind_cols(as.data.frame(t(as.matrix(sce_merge_res1$sce@assays$RNA@data[c("Tnf","Ifng","Gzmb","Gzma","Havcr2"),])))) %>% 
                  reshape2::melt(.,id.vars=c("cellID","orig.ident","celltype","celltype1" )) %>% 
                  data.table::setnames(.,old = c("variable","value"),new = c("gene","exp")) %>% 
                  mutate(orig.ident=factor(orig.ident,levels = c("Vehicle","Sofa","R848","SR3")))
   
    p <- metadata %>% 
               group_by(orig.ident, celltype, gene) %>% 
               dplyr::summarise(mean=mean(exp)) %>% 
               right_join(.,metadata,by=c("orig.ident","gene","celltype")) %>% 
               filter(.,celltype %in% c("CD4+T_cell","CD8+T_cell","Neutrophil_cell","Marcrophage")) %>% 
                ggplot(mapping=aes(x=orig.ident,y=exp,fill=mean)) +
                geom_violin(scale = "width", adjust = 2, trim = TRUE) +
                scale_fill_gradient(low = "#F2E80D", high = "#E65224") +
                ggsignif::geom_signif(comparisons = GetComb(c("Vehicle","R848","Sofa","SR3")),test = "t.test",step_increase = 0.1) +
                facet_grid(gene~celltype,scales = "free") +
                cowplot::theme_cowplot(font_size = 12) +
                theme(panel.spacing = unit(0, "lines"),
                      panel.background = element_rect(fill = NA, color = "black"),
                      strip.background = element_blank(),
                      strip.text = element_text(face = "bold"),
                      strip.text.y.left = element_text(angle = 0),
                      axis.text.x = element_text(angle=60, vjust=1,hjust=1)) +
                coord_flip()    
    ggsave(p,filename = "figure7_violin_plot.pdf",width = 12,height=14)
    
    metadata %>% 
      group_by(orig.ident, celltype, gene) %>% 
      dplyr::summarise(mean=mean(exp)) %>% 
      right_join(.,metadata,by=c("orig.ident","gene","celltype")) %>% 
      filter(.,celltype=="Neutrophil_cell",gene=="Ifng") %>% 
      filter(.,orig.ident %in% c("R848","Vehicle")) %>% 
      t.test(exp~orig.ident,data=.)
    
    
  }

  #T细胞精确注释
  {
    T_sce <- subset(sce_merge_res1$sce,celltype1=="T_cell")
    T_sce <- FastSeurat(count=NULL,
                        dir.name=NULL,
                        obj=T_sce,#min.cells=3,
                        species=c("human","mouse")[2],
                        min.features=0,
                        max.features=300000,
                        percent.mt.num=100,
                        plot=F,
                        pcSelect=20,
                        project="mouse",
                        nfeatures=2000,
                        vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                        all.scale = F,
                        npcs = 50,
                        resolution = 0.5,
                        harmony=F,
                        doublet=F,
                        isMarkers=T,
                        cellCycle=T,
                        rmOtherGene=T,
                        perplexity = 30,
                        features=NULL,
                        filepath=NULL,
                        outdir="Results",
                        names="love")
    p1 <- dittoDimPlot(T_sce$sce, reduction.use  = "umap", var="seurat_clusters", size = 1, do.label = T, 
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "T cells") + theme(legend.text = element_text(face="bold",size = 12))
    p2 <- dittoDimPlot(T_sce$sce, reduction.use  = "umap", var="orig.ident", size = 1, do.label = TRUE, 
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "T cells") + theme(legend.text = element_text(face="bold",size = 12))
    p1 | p2
    cl_markers <- T_sce$sce.markers
    top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    dh <- DoHeatmap(T_sce$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
    dh <- DoHeatmap(T_sce$sce, features = top10_cl_markers$gene) + NoLegend()
    print(dh)
    Epithelial_cell <- c("Krt8","Mmp7","Krt18")
    Stromal_cell <- c("Col1a2","Col3a1","Col14a1","Fbln2")
    All.Immune <- "Ptprc" #免疫细胞
    B.cell <- c("Cd19","Cd79a","Cd79b","Ms4a1")
    T.cell <- c("Cd3d","Cd3e","Cd8a","Cd4")
    Treg <- c("Foxp3","Il2ra","Ikzf2","Stat5","Il7r","Il10","Tgfb1","Trac","Trbc2","Foxp3","Il32","Ctla4","Tnfrsf18")
    NKT <- c("Cd56","Cd49b","Zbtb16","Tbx21","Klrc3","Klrc2","Klrc1","Klrb1c","Nkg7","Cd3","Klrb1","Fcgr3a","Ncam","Zbtb16","Tyrobp")
    NK <- c("Gzma","Klrb1c","Ncr1")
    Mast <- c("Kit","Mcpt4","Fcer1a")
    Neutrophil <- c("S100a9","Lcn2","S100a8")
    Dendritic_cell <- c("Ccr7","Plet1","Wdfy4")
    Marcrophage <- c("Apoe","Lyz2")                                                                                                                                                                                                                                                                                                                                                                
    FeaturePlot_scCustom(T_sce$sce, features = Epithelial_cell, order = F)
    FeaturePlot_scCustom(T_sce$sce, features = Stromal_cell, order = F)
    FeaturePlot_scCustom(T_sce$sce, features = All.Immune, order = F)
    FeaturePlot_scCustom(T_sce$sce, features = B.cell, order = F)
    FeaturePlot_scCustom(T_sce$sce, features = T.cell, order = F)
    FeaturePlot_scCustom(T_sce$sce, features = Treg, order = F)
    FeaturePlot_scCustom(T_sce$sce, features = NKT, order = F)
    FeaturePlot_scCustom(T_sce$sce, features = NK, order = F)
    FeaturePlot_scCustom(T_sce$sce, features = Mast, order = F)
    FeaturePlot_scCustom(T_sce$sce, features = Neutrophil, order = F)
    FeaturePlot_scCustom(T_sce$sce, features = Dendritic_cell, order = F)
    FeaturePlot_scCustom(T_sce$sce, features = Marcrophage, order = F)
    #去除双细胞9
    T_sce1 <- subset(T_sce$sce,seurat_clusters != 9)
    
    cl_markers <- FindAllMarkers(T_sce1)
    top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    dh <- DoHeatmap(T_sce1, features = top10_cl_markers$gene) + NoLegend()
    print(dh)
    p1 <- dittoDimPlot(T_sce1, reduction.use  = "umap", var="seurat_clusters", size = 1, do.label = T, 
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "T cells") + theme(legend.text = element_text(face="bold",size = 12))
    p2 <- dittoDimPlot(T_sce1, reduction.use  = "umap", var="orig.ident", size = 1, do.label = TRUE, 
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "T cells") + theme(legend.text = element_text(face="bold",size = 12))
    p1 | p2
    proliferative <- c("Mki67","Pcna","Tk1","Top2a","Tubb","Tyms") #增殖，细胞复制S期，处于增殖阶段
    CD8.T <- c("Cd8a","Gzmk","Gzmb","Pik3r1","Tuba4a","Tnfrsf9","Top2a","Mki67")
    CD4.T <- c("Cd4","Il7r")
    Tn <- c("Tcf7","Ccr7","Sell","IL7R","IL2RG","LEF1","S1PR1","IL2RG") #naive
    Tcm <- c("Gpr183","Ccl5","Gzma","Prf1","Ccr7", "Sell","Il7r","S1pr1","Ptger2","Anxa1") #memory
    Trm <- c("ITGAE","CD69","ZNF683","ITGA1","IFNG","GZMH","PRF1","KLRG1","CD44","PRDM1","XCL1","XCL2","ITGA1","CXCR6","ENTPD1","RUNX3")
    Treg <- c("FOXP3","BATF","CTLA4", "TIGIT", "IL2RA", "CCR8", "IL10RA", "IKZF2", "RTKN2","CDC25B", "TNFRSF4" ,"TNFRSF18", "SAT1", "NRP1" ,"LAG3", "ENTPD1", "CD73")
    Teff <- c("GZMA","GZMK","NKG7","IFNG","PRF1","GZMB","CX3CR1","FCGR3A","FGFBP2","KLRG1","GZMH","TBX21","EOMES","GNLY") #毒性
    T_early_after_activation <- c("JUNB","FOS", "ATF3", "HSPA1A" ,"DNAJB1")
    Gama_delt_T<- c("Trdc","Trgv9" ,"Klrb1")
    Exhausted.cell <- c("LAG3","TIGIT","TIM3","HAVCR1","PDCD1","HAVCR2","TOX","NR4A1","CTLA4","TIGIT","LAYN","TNFRSF9","CXCL13","TNFRSF18")
    Th1 <- c("IFNG","TBX21","STAT4","STAT1","CD119","CXCR3","CXCR6","CCR1","CCR5")
    MAIT <- c("KLRB1","SLC4A10","ZBTB16","NCR3","TRAV1-2")
    Th1 <- c("IFNG","TBX21","STAT4","STAT1","CD119","CXCR3","CXCR6","CCR1","CCR5","RANKL","TIM3")
    Th2 <- c("STAT6","GATA3","CCR4","PTGDR2","IL4","IL5","IL13")
    Tfh <- c("IL6ST","STAT3","BCL6","CXCR5","IL21","CXCL13","CCR7","ICOS","PDCD1","CD200","ICA1","TOX","TOX2","MAGEH1","BTLA")
    Th17 <- c("IL17A","IL23R","RORC","KLRB1","CCR6","CAPG","ITGAE","FURIN","CTSH")
    NK_like <- c("KLRB1","FCGR3A","GZMH","GNLY")
    NKT <- c("CD3","KLRB1","FCGR3A","NCAM","ZBTB16","PLZF","TYROBP")  
    FeaturePlot_scCustom(T_sce1, features = Epithelial_cell, order = F)
    FeaturePlot_scCustom(T_sce1, features = Stromal_cell, order = F)
    FeaturePlot_scCustom(T_sce1, features = All.Immune, order = F)
    FeaturePlot_scCustom(T_sce1, features = B.cell, order = F)
    FeaturePlot_scCustom(T_sce1, features = T.cell, order = F)
    FeaturePlot_scCustom(T_sce1, features = Gama_delt_T, order = F)
    FeaturePlot_scCustom(T_sce1, features = Treg, order = F)
    FeaturePlot_scCustom(T_sce1, features = NKT, order = F)
    FeaturePlot_scCustom(T_sce1, features = NK, order = F)
    FeaturePlot_scCustom(T_sce1, features = Mast, order = F)
    FeaturePlot_scCustom(T_sce1, features = Neutrophil, order = F)
    FeaturePlot_scCustom(T_sce1, features = Dendritic_cell, order = F)
    FeaturePlot_scCustom(T_sce1, features = Marcrophage, order = F)
    T_sce1@meta.data$celltype <- as.character(plyr::mapvalues(T_sce1@meta.data$seurat_clusters,from = c(0:8,10:11),
                                                                  to=c("CD8+T","CD4+T","CD8+T","CD8+T","CD8+T","CD8+T","Treg","CD4+T","gama_T","CD8+T","NKT")))
    saveRDS(T_sce1,file="./Results/mouse/T_cell.rds")
    T_sce1 <- readRDS(file = "./Results/mouse/T_cell.rds")
    p3 <- dittoDimPlot(T_sce1, reduction.use  = "umap", var="celltype", size = 1, do.label = TRUE, 
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "T cells") + theme(legend.text = element_text(face="bold",size = 12))
    p2 <- dittoDimPlot(T_sce1, reduction.use  = "umap", var="orig.ident", size = 1, do.label = TRUE, 
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "T cells") + theme(legend.text = element_text(face="bold",size = 12))
    p1 <- dittoDimPlot(T_sce1, reduction.use  = "umap", var="seurat_clusters", size = 1, do.label = TRUE, 
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "T cells") + theme(legend.text = element_text(face="bold",size = 12))
    p1 | p2 | p3
    #特征基因展示
    feature_genes <- c("Cd8a","Cd4","Foxp3","Trdc","Gzma","Klrb1c")
    FeaturePlot_scCustom(T_sce1, features =  feature_genes, order = F)
    #统计四组Treg细胞比例
    T_sce1@meta.data %>% 
                   dplyr::group_by(orig.ident) %>% 
                   dplyr::summarise(proportion=as.numeric(table(celltype)/length(celltype)),names=str_replace_all(names(table(celltype))," ","_")) %>% 
                   ungroup() %>% as.data.frame() %>% 
                   dplyr::filter(.,names=="Treg") %>% 
                   dplyr::mutate(orig.ident = factor(orig.ident,levels=c("Vehicle","Sofa","R848","SR3"))) -> treg.plot
    p <- ggplot(data = treg.plot,mapping = aes(x=orig.ident,y=proportion,fill=orig.ident)) + 
      geom_bar(stat="identity") +
      #facet_wrap(vars(names), nrow=1,scales = "free_y") +
      scale_fill_manual(values = dittoColors()) +
      theme_bw()
    ggsave(p,filename = "Treg_proportion.pdf")
    T_sce1@meta.data$celltype[which(T_sce1@meta.data$celltype != "Treg")] <- "non_treg"
    chisq.test(table(T_sce1@meta.data$orig.ident,T_sce1@meta.data$celltype))
    
    #Vehicle和Sofa比较
    pos <- Fastgrep2(c("Vehicle","Sofa"),T_sce1@meta.data$orig.ident) 
    aa <- T_sce1@meta.data[unlist(pos),]
    table(aa$orig.ident)
    table(aa$celltype)
    aa$celltype <- ifelse(aa$celltype=="Treg","Treg","others")
    fisher.test(table(aa$orig.ident,aa$celltype),alternative = "greater")  #"greater" or "less"
    
    #Vehicle和R848比较
    pos <- Fastgrep2(c("Vehicle","R848"),T_sce1@meta.data$orig.ident) 
    aa <- T_sce1@meta.data[unlist(pos),]
    table(aa$orig.ident)
    table(aa$celltype)
    aa$celltype <- ifelse(aa$celltype=="Treg","Treg","others")
    fisher.test(table(aa$orig.ident,aa$celltype),alternative = "greater")
    
    #Vehicle和SR3比较
    pos <- Fastgrep2(c("Vehicle","SR3"),T_sce1@meta.data$orig.ident) 
    aa <- T_sce1@meta.data[unlist(pos),]
    table(aa$orig.ident)
    table(aa$celltype)
    aa$celltype <- ifelse(aa$celltype=="Treg","Treg","others")
    fisher.test(table(aa$orig.ident,aa$celltype))
    
    #Sofa和SR3比较
    pos <- Fastgrep2(c("Sofa","SR3"),T_sce1@meta.data$orig.ident) 
    aa <- T_sce1@meta.data[unlist(pos),]
    table(aa$orig.ident)
    table(aa$celltype)
    aa$celltype <- ifelse(aa$celltype=="Treg","Treg","others")
    fisher.test(table(aa$orig.ident,aa$celltype))
    
    #R848和SR3比较
    pos <- Fastgrep2(c("R848","SR3"),T_sce1@meta.data$orig.ident) 
    aa <- T_sce1@meta.data[unlist(pos),]
    table(aa$orig.ident)
    table(aa$celltype)
    aa$celltype <- ifelse(aa$celltype=="Treg","Treg","others")
    fisher.test(table(aa$orig.ident,aa$celltype))
  }
  
  #对髓样细胞精细分群注释
  {
    myeloid <- subset(sce_merge_res1$sce,celltype1 == "Dendritic_cell" | celltype1 == "Marcrophage" | celltype1 == "Mast_cell" | celltype1 == "Neutrophil_cell")
    myeloid <- FastSeurat(count=NULL,
                          dir.name=NULL,
                          obj=myeloid,#min.cells=3,
                          species=c("human","mouse")[2],
                          min.features=0,
                          max.features=300000,
                          percent.mt.num=100,
                          plot=F,
                          pcSelect=30,
                          project="mouse",
                          nfeatures=2000,
                          vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                          all.scale = F,
                          npcs = 50,
                          resolution = 0.2,
                          harmony=F,
                          doublet=F,
                          isMarkers=T,
                          cellCycle=T,
                          rmOtherGene=T,
                          perplexity = 30,
                          features=NULL,
                          filepath=NULL,
                          outdir="Results",
                          names="love")
    p1 <- dittoDimPlot(myeloid$sce, reduction.use  = "umap", var="seurat_clusters", size = 1, do.label = TRUE, 
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "myeloid cells") + theme(legend.text = element_text(face="bold",size = 12))
    p2 <- dittoDimPlot(myeloid$sce, reduction.use  = "umap", var="orig.ident", size = 1, do.label = TRUE, 
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "myeloid cells") + theme(legend.text = element_text(face="bold",size = 12))
    p1 | p2
    cl_markers <- myeloid$sce.markers
    top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    dh <- DoHeatmap(myeloid$sce, features = top10_cl_markers$gene) + NoLegend()
    print(dh)
    Epithelial_cell <- c("Krt8","Mmp7","Krt18")
    Stromal_cell <- c("Col1a2","Col3a1","Col14a1","Fbln2")
    All.Immune <- "Ptprc" #免疫细胞
    B.cell <- c("Cd19","Cd79a","Cd79b","Ms4a1")
    T.cell <- c("Cd3d","Cd3e","Cd8a","Cd4")
    Treg <- c("Foxp3","Il2ra","Ikzf2","Stat5","Il7r","Il10","Tgfb1","Trac","Trbc2","Foxp3","Il32","Ctla4","Tnfrsf18")
    NKT <- c("Cd56","Cd49b","Zbtb16","Tbx21","Klrc3","Klrc2","Klrc1","Klrb1c","Nkg7","Cd3","Klrb1","Fcgr3a","Ncam","Zbtb16","Tyrobp")
    NK <- c("Gzma","Klrb1c","Ncr1")
    Mast <- c("Kit","Mcpt4","Fcer1a")
    Neutrophil <- c("S100a9","Lcn2","S100a8")
    Dendritic_cell <- c("Ccr7","Plet1","Wdfy4")
    Marcrophage <- c("Apoe","Lyz2")
    FeaturePlot_scCustom(myeloid$sce, features = Epithelial_cell, order = F)
    FeaturePlot_scCustom(myeloid$sce, features = Stromal_cell, order = F)
    FeaturePlot_scCustom(myeloid$sce, features = All.Immune, order = F)
    FeaturePlot_scCustom(myeloid$sce, features = B.cell, order = F)
    FeaturePlot_scCustom(myeloid$sce, features = T.cell, order = F)
    FeaturePlot_scCustom(myeloid$sce, features = Treg, order = F)
    FeaturePlot_scCustom(myeloid$sce, features = NKT, order = F)
    FeaturePlot_scCustom(myeloid$sce, features = NK, order = F)
    FeaturePlot_scCustom(myeloid$sce, features = Mast, order = F)
    FeaturePlot_scCustom(myeloid$sce, features = Neutrophil, order = F)
    FeaturePlot_scCustom(myeloid$sce, features = Dendritic_cell, order = F)
    FeaturePlot_scCustom(myeloid$sce, features = Marcrophage, order = F)
    #去除cluster3双细胞重新分群 
    myeloid_sce <- subset(myeloid$sce,seurat_clusters != 3)
    myeloid_sce <- FastSeurat(count=NULL,
                              dir.name=NULL,
                              obj=myeloid_sce,#min.cells=3,
                              species=c("human","mouse")[2],
                              min.features=0,
                              max.features=300000,
                              percent.mt.num=100,
                              plot=F,
                              pcSelect=30,
                              project="mouse",
                              nfeatures=2000,
                              vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                              all.scale = F,
                              npcs = 50,
                              resolution = 0.2,
                              harmony=F,
                              doublet=F,
                              isMarkers=T,
                              cellCycle=T,
                              rmOtherGene=T,
                              perplexity = 30,
                              features=NULL,
                              filepath=NULL,
                              outdir="Results",
                              names="love")
    saveRDS(myeloid_sce,file = "./Results/mouse/Myeloid.rds")
    myeloid_sce <- readRDS(file = "./Results/mouse/Myeloid.rds")
    Epithelial_cell <- c("Krt8","Mmp7","Krt18")
    Stromal_cell <- c("Col1a2","Col3a1","Col14a1","Fbln2")
    All.Immune <- "Ptprc" #免疫细胞
    B.cell <- c("Cd19","Cd79a","Cd79b","Ms4a1")
    T.cell <- c("Cd3d","Cd3e","Cd8a","Cd4")
    Treg <- c("Foxp3","Il2ra","Ikzf2","Stat5","Il7r","Il10","Tgfb1","Trac","Trbc2","Foxp3","Il32","Ctla4","Tnfrsf18")
    NKT <- c("Cd56","Cd49b","Zbtb16","Tbx21","Klrc3","Klrc2","Klrc1","Klrb1c","Nkg7","Cd3","Klrb1","Fcgr3a","Ncam","Zbtb16","Tyrobp")
    NK <- c("Gzma","Klrb1c","Ncr1")
    Mast <- c("Kit","Mcpt4","Fcer1a")
    Neutrophil <- c("S100a9","Lcn2","S100a8")
    Dendritic_cell <- c("Ccr7","Plet1","Wdfy4")
    Marcrophage <- c("Apoe","Lyz2")
    FeaturePlot_scCustom(myeloid_sce$sce, features = Epithelial_cell, order = F)
    FeaturePlot_scCustom(myeloid_sce$sce, features = Stromal_cell, order = F)
    FeaturePlot_scCustom(myeloid_sce$sce, features = All.Immune, order = F)
    FeaturePlot_scCustom(myeloid_sce$sce, features = B.cell, order = F)
    FeaturePlot_scCustom(myeloid_sce$sce, features = T.cell, order = F)
    FeaturePlot_scCustom(myeloid_sce$sce, features = Treg, order = F)
    FeaturePlot_scCustom(myeloid_sce$sce, features = NKT, order = F)
    FeaturePlot_scCustom(myeloid_sce$sce, features = NK, order = F)
    FeaturePlot_scCustom(myeloid_sce$sce, features = Mast, order = F)
    FeaturePlot_scCustom(myeloid_sce$sce, features = Neutrophil, order = F)
    FeaturePlot_scCustom(myeloid_sce$sce, features = Dendritic_cell, order = F)
    FeaturePlot_scCustom(myeloid_sce$sce, features = Marcrophage, order = F)
    myeloid_features <- c("Lyz2","Mcpt4","S100a9","S100a8","ccr7","Plet1")
    FeaturePlot_scCustom(myeloid_sce$sce, features = myeloid_features, order = F)
    myeloid_sce$sce@meta.data$celltype <- as.character(plyr::mapvalues(myeloid_sce$sce@meta.data$seurat_clusters,from = c(0:7),
                                                              to=c("Neutrophil","Neutrophil","Marcrophage","Dendritic","Dendritic","Mast","6","Neutrophil")))
    p1 <- dittoDimPlot(myeloid_sce$sce, reduction.use  = "umap", var="seurat_clusters", size = 1, do.label = TRUE, 
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "myeloid cells") + theme(legend.text = element_text(face="bold",size = 12))
    p2 <- dittoDimPlot(myeloid_sce$sce, reduction.use  = "umap", var="orig.ident", size = 1, do.label = TRUE, 
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "myeloid cells") + theme(legend.text = element_text(face="bold",size = 12))
    p3 <- dittoDimPlot(myeloid_sce$sce, reduction.use  = "umap", var="celltype", size = 1, do.label = TRUE, 
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "myeloid cells") + theme(legend.text = element_text(face="bold",size = 12))
    p1 | p2 | p3
    head(myeloid_sce$sce@meta.data)
    p1 <- dittoDimPlot(myeloid_sce$sce, reduction.use  = "umap", var="celltype", size = 1, do.label = TRUE, cells.use = grepl("Vehicle",myeloid_sce$sce@meta.data$orig.ident) ,
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "Vehicle") + theme(legend.text = element_text(face="bold",size = 12)) & NoLegend()
    p2 <- dittoDimPlot(myeloid_sce$sce, reduction.use  = "umap", var="celltype", size = 1, do.label = TRUE, cells.use = grepl("Sofa",myeloid_sce$sce@meta.data$orig.ident) ,
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "Sofa") + theme(legend.text = element_text(face="bold",size = 12)) & NoLegend()
    p3 <- dittoDimPlot(myeloid_sce$sce, reduction.use  = "umap", var="celltype", size = 1, do.label = TRUE, cells.use = grepl("R848",myeloid_sce$sce@meta.data$orig.ident) ,
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "R848") + theme(legend.text = element_text(face="bold",size = 12)) & NoLegend()
    p4 <- dittoDimPlot(myeloid_sce$sce, reduction.use  = "umap", var="celltype", size = 1, do.label = TRUE, cells.use = grepl("SR3",myeloid_sce$sce@meta.data$orig.ident) ,
                       do.ellipse = F, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "SR3") + theme(legend.text = element_text(face="bold",size = 12))
    p <- p1 | p2 | p3 | p4
    ggsave(p,filename = "figure5.myeloid_umap.pdf",width=16,height=5)
    cl_markers <- myeloid_sce$sce.markers
    top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    dh <- DoHeatmap(myeloid_sce$sce, features = top10_cl_markers$gene) + NoLegend()
    print(dh)
    
    ##############   提取巨噬细胞，计算巨噬细胞marker得分，查看我们的巨噬细胞倾向于什么  #############
    Marcrophage_sce <- subset(myeloid_sce$sce,celltype == "Marcrophage")
    Marcrophage_sce <- FastSeurat(count=NULL,
                              dir.name=NULL,
                              obj=Marcrophage_sce,#min.cells=3,
                              species=c("human","mouse")[2],
                              min.features=0,
                              max.features=300000,
                              percent.mt.num=100,
                              plot=F,
                              pcSelect=30,
                              project="mouse",
                              nfeatures=2000,
                              vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                              all.scale = F,
                              npcs = 50,
                              resolution = 0.1,
                              harmony=F,
                              doublet=F,
                              isMarkers=T,
                              cellCycle=T,
                              rmOtherGene=T,
                              perplexity = 30,
                              features=NULL,
                              filepath=NULL,
                              outdir="Results",
                              names="love")
    Epithelial_cell <- c("Krt8","Mmp7","Krt18")
    Stromal_cell <- c("Col1a2","Col3a1","Col14a1","Fbln2")
    All.Immune <- "Ptprc" #免疫细胞
    B.cell <- c("Cd19","Cd79a","Cd79b","Ms4a1")
    T.cell <- c("Cd3d","Cd3e","Cd8a","Cd4")
    Treg <- c("Foxp3","Il2ra","Ikzf2","Stat5","Il7r","Il10","Tgfb1","Trac","Trbc2","Foxp3","Il32","Ctla4","Tnfrsf18")
    NKT <- c("Cd56","Cd49b","Zbtb16","Tbx21","Klrc3","Klrc2","Klrc1","Klrb1c","Nkg7","Cd3","Klrb1","Fcgr3a","Ncam","Zbtb16","Tyrobp")
    NK <- c("Gzma","Klrb1c","Ncr1")
    Mast <- c("Kit","Mcpt4","Fcer1a")
    Neutrophil <- c("S100a9","Lcn2","S100a8")
    Dendritic_cell <- c("Ccr7","Plet1","Wdfy4")
    Marcrophage <- c("Apoe","Lyz2")
    FeaturePlot_scCustom(Marcrophage_sce$sce, features = Epithelial_cell, order = F)
    FeaturePlot_scCustom(Marcrophage_sce$sce, features = Stromal_cell, order = F)
    FeaturePlot_scCustom(Marcrophage_sce$sce, features = All.Immune, order = F)
    FeaturePlot_scCustom(Marcrophage_sce$sce, features = B.cell, order = F)
    FeaturePlot_scCustom(Marcrophage_sce$sce, features = T.cell, order = F)
    FeaturePlot_scCustom(Marcrophage_sce$sce, features = Treg, order = F)
    FeaturePlot_scCustom(Marcrophage_sce$sce, features = NKT, order = F)
    FeaturePlot_scCustom(Marcrophage_sce$sce, features = NK, order = F)
    FeaturePlot_scCustom(Marcrophage_sce$sce, features = Mast, order = F)
    FeaturePlot_scCustom(Marcrophage_sce$sce, features = Neutrophil, order = F)
    FeaturePlot_scCustom(Marcrophage_sce$sce, features = Dendritic_cell, order = F)
    FeaturePlot_scCustom(Marcrophage_sce$sce, features = Marcrophage, order = F)
    p1 <- dittoDimPlot(Marcrophage_sce$sce, reduction.use  = "umap", var="seurat_clusters", size = 1, do.label = TRUE, 
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "myeloid cells") + theme(legend.text = element_text(face="bold",size = 12))
    p2 <- dittoDimPlot(Marcrophage_sce$sce, reduction.use  = "umap", var="orig.ident", size = 1, do.label = TRUE, 
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "myeloid cells") + theme(legend.text = element_text(face="bold",size = 12))
    p3 <- dittoDimPlot(Marcrophage_sce$sce, reduction.use  = "umap", var="celltype", size = 1, do.label = TRUE, 
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "myeloid cells") + theme(legend.text = element_text(face="bold",size = 12))
    p1 | p2 | p3
    
    p1 <- dittoDimPlot(Marcrophage_sce$sce, reduction.use  = "umap", var="celltype2", size = 1, do.label = TRUE, cells.use = grepl("Vehicle",myeloid_sce$sce@meta.data$orig.ident) ,
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "Vehicle") + theme(legend.text = element_text(face="bold",size = 12))
    p2 <- dittoDimPlot(Marcrophage_sce$sce, reduction.use  = "umap", var="celltype2", size = 1, do.label = TRUE, cells.use = grepl("Sofa",myeloid_sce$sce@meta.data$orig.ident) ,
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "Sofa") + theme(legend.text = element_text(face="bold",size = 12))
    p3 <- dittoDimPlot(Marcrophage_sce$sce, reduction.use  = "umap", var="celltype2", size = 1, do.label = TRUE, cells.use = grepl("R848",myeloid_sce$sce@meta.data$orig.ident) ,
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "R848") + theme(legend.text = element_text(face="bold",size = 12))
    p4 <- dittoDimPlot(Marcrophage_sce$sce, reduction.use  = "umap", var="celltype2", size = 1, do.label = TRUE, cells.use = grepl("SR3",myeloid_sce$sce@meta.data$orig.ident) ,
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "SR3") + theme(legend.text = element_text(face="bold",size = 12))
    p1 + p2 + p3 + p4
    cl_markers <- Marcrophage_sce$sce.markers
    top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    dh <- DoHeatmap(myeloid_sce$sce, features = top10_cl_markers$gene) + NoLegend()
    print(dh)
    #去除巨噬细胞里的双细胞cluster2
    Marcrophage_sce <- subset(Marcrophage_sce$sce,seurat_clusters != 2)
    saveRDS(Marcrophage_sce,file="./Results/mouse/Marcrophage.rds")
    Marcrophage_sce <- readRDS(file="./Results/mouse/Marcrophage.rds")
    cl_markers <- FindAllMarkers(Marcrophage_sce)
    top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
    dh <- DoHeatmap(Marcrophage_sce, features = top10_cl_markers$gene, slot="data")
    print(dh)
    ggsave(dh,filename = "./code/tianjin_Supple/singlecell_mouse/pictures/marker_heatmap.pdf",width=10,height=6)
    p1 <- dittoDimPlot(Marcrophage_sce, reduction.use  = "umap", var="seurat_clusters", size = 1, do.label = TRUE, 
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "Marcrophage") + theme(legend.text = element_text(face="bold",size = 12))
    p2 <- dittoDimPlot(Marcrophage_sce, reduction.use  = "umap", var="orig.ident", size = 1, do.label = TRUE, 
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "Marcrophage") + scale_color_manual(values = c(Vehicle = "#3A8757", Sofa = "#ADDEAD", R848 = "#FFD0A0",SR3="#FF8000")) + theme(legend.text = element_text(face="bold",size = 12))
    p <- p1 | p2
    ggsave(p2,filename = "figure5_marcrofhage_umap1.pdf")
    eset <- Marcrophage_sce@assays$RNA@data
    M1_down <- format_msigdb("/Users/biofly/project/shijian/data/GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN.gmt", ont = "term", gene = "gene")
    M1_up <- format_msigdb("/Users/biofly/project/shijian/data/GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP.gmt", ont = "term", gene = "gene")
    M1_down$GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN <- str_to_title(unlist(M1_down))
    M1_up$GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP <- str_to_title(unlist(M1_up))
    signature <- c(M1_down,M1_up)
    Marcrophage_zscore <- calculate_sig_score_zscore(pdata = NULL,
                                                     eset=as.matrix(eset),
                                                     signature=signature,
                                                     mini_gene_count=2,
                                                     column_of_sample="ID",
                                                     adjust_eset = FALSE)
    Marcrophage_ssgsea <- calculate_sig_score_ssgsea(pdata = NULL,
                                                     eset=as.matrix(eset),
                                                     signature=signature,
                                                     mini_gene_count=2,
                                                     column_of_sample="ID",
                                                     adjust_eset = FALSE)
    Marcrophage_sce@meta.data %>%
                              tibble::rownames_to_column("cell") %>% 
                              left_join(.,Marcrophage_zscore,by=c("cell"="ID")) -> meta.data
    rownames(meta.data) <- meta.data$cell
    plot.data <- melt(Marcrophage_zscore,id.vars=c("Index","ID")) %>% 
          mutate(sample=str_split(ID,"_",simplify = T)[,1]) %>% 
          left_join(.,meta.data[,c(1,2,7,17,18)],by=c("ID"="cell")) 
    p <- ggplot(data=plot.data,mapping = aes(x=seurat_clusters,y=value,fill=seurat_clusters)) +
      geom_boxplot(outlier.shape = NA) +  # 不显示离群点
      ggsignif::geom_signif(comparisons = list(unique(as.character(plot.data$seurat_clusters))),test = t.test)+
      geom_jitter(size=.5) + 
      facet_wrap(vars(variable), scales = "free_y") + 
      theme_bw() +
      scale_fill_manual(values = dittoColors())
    ggsave(p,filename = "figure5_M1_M2_zscore.pdf")
      #cowplot::theme_cowplot(font_size = 14)
    plot.data1 <- melt(Marcrophage_ssgsea,id.vars=c("Index","ID")) %>% 
      mutate(sample=str_split(ID,"_",simplify = T)[,1]) %>% 
      left_join(.,meta.data[,c(1,2,7,17,18)],by=c("ID"="cell")) 
    ggplot(data=plot.data1,mapping = aes(x=seurat_clusters,y=value,fill=seurat_clusters)) +
      geom_boxplot(outlier.shape = NA) +  # 不显示离群点
      geom_signif(comparisons = list(unique(as.character(plot.data1$seurat_clusters))),test = t.test)+
      geom_jitter(size=.5) + 
      facet_wrap(vars(variable), scales = "free_y")+
      theme_cowplot(font_size = 14)
   Marcrophage_sce@reductions$umap@cell.embeddings %>% 
         as.data.frame() %>% 
         tibble::rownames_to_column("cell") %>% 
         left_join(.,Marcrophage_zscore,by=c("cell"="ID")) %>% 
         mutate(type=ifelse(GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP > GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN,"M1","M2")) %>% 
          ggplot(data=.,mapping = aes(x=UMAP_1,y=UMAP_2,color=GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN)) +
          geom_point() +
          theme_bw() +
          scale_color_gradientn(colours = rainbow(2))
   Marcrophage_sce@reductions$umap@cell.embeddings %>% 
                                       as.data.frame() %>% 
                                       tibble::rownames_to_column("cell") %>% left_join(.,Marcrophage_zscore,by=c("cell"="ID")) %>% 
                                       mutate(type=ifelse(GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP > GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN,"M1","M2")) %>% 
                                        ggplot(data=.,mapping = aes(x=UMAP_1,y=UMAP_2,color=GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP)) +
                                        geom_point() +
                                        theme_bw() +
                                        scale_color_gradientn(colours = rainbow(2))
                                        #scale_color_continuous(type = "viridis")
    ggplot(data=meta.data,mapping = aes(x=factor(orig.ident,levels=c("Vehicle","Sofa","R848","SR3")),y=GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN,color=orig.ident)) +
      geom_boxplot() +
      labs(x="orig.ident") +
      theme_minimal_grid()
    ggplot(data=meta.data,mapping = aes(x=orig.ident,y=GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP,color=orig.ident)) +
      geom_boxplot() +
      labs(x="orig.ident") +
      theme_minimal_grid()
  
    #M2巨噬细胞比例
    meta.data %>% dplyr::group_by(orig.ident) %>% 
      dplyr::summarise(proportion=as.numeric(table(seurat_clusters)/length(seurat_clusters)),names=str_replace_all(names(table(seurat_clusters))," ","_")) %>% 
      filter(.,names==0) %>% 
      mutate(orig.ident=factor(orig.ident,levels = c("Vehicle","Sofa","R848","SR3"))) %>% 
      ggplot(mapping = aes(x=orig.ident,y=proportion,fill=orig.ident)) + 
      geom_bar(stat="identity") +
      #facet_wrap(vars(names), nrow=1,scales = "free_y") +
      scale_fill_manual(values = dittoColors()) +
      theme_bw() 
    
    #Vehicle vs Sofa巨噬细胞比例差异
    head(meta.data)
    meta.data %>% mutate(celltype2=if_else(seurat_clusters==1,"M1","M2")) -> meta.data
    aa <- meta.data %>% filter(.,orig.ident=="Vehicle" | orig.ident=="Sofa") 
    fisher.test(table(aa$orig.ident,aa$celltype2),alternative = "less")
    aa <- meta.data %>% filter(.,orig.ident=="Vehicle" | orig.ident=="R848") 
    fisher.test(table(aa$orig.ident,aa$celltype2),alternative = "less")
    aa <- meta.data %>% filter(.,orig.ident=="Vehicle" | orig.ident=="SR3") 
    fisher.test(table(aa$orig.ident,aa$celltype2),alternative = "less")
    
    #这个图对于cluster0和1，分别以每个cluster为单位，四组所占的比例做一个barplot或者piechart。
    tt <- Marcrophage_sce@meta.data %>% 
                mutate(seurat_clusters=paste0("cluster",as.character(seurat_clusters))) %>% 
                dplyr::group_by(seurat_clusters) %>% 
                dplyr::summarise(proportion=table(orig.ident)/length(orig.ident),names=names(table(orig.ident))) %>% 
                mutate(proportion=as.numeric(proportion)) %>% 
                mutate(labels=scales::percent(proportion)) 
      
   
    ggplot(data= tt,mapping=aes(x="",y=proportion,fill=names)) +
                    geom_col(color="black") +
                    geom_text(aes(label = labels),
                            position = position_stack(vjust = 0.5)) +
                    scale_fill_viridis_d() +
                    coord_polar(theta = "y") +
                    facet_wrap(vars(seurat_clusters)) +
                    #cowplot::theme_cowplot(font_size = 14) +
                    theme_void()
                
    p <- ggplot(data= tt,mapping=aes(x="",y=proportion,fill=names)) +
      geom_bar(stat="identity", width=1, color="white") +
      geom_text(aes(label = labels),
                position = position_stack(vjust = 0.5)) +
      scale_fill_viridis_d() +
      coord_polar(theta = "y",start = 0) +
      facet_wrap(vars(seurat_clusters)) +
      cowplot::theme_cowplot(font_size = 14)
    ggsave(p,filename="figure5_pie_plot.pdf")
  }
  
  #提取巨噬细胞cluster2进行GSVA分析
  {
    Marcrophage_sce <- subset(myeloid_sce$sce,seurat_clusters == 2)
    mar_metadata <- Marcrophage_sce@meta.data %>% tibble::rownames_to_column("cellID") %>% 
      dplyr::select(.,c(1,2,7,11,12)) %>% 
      bind_cols(as.data.frame(t(as.matrix(Marcrophage_sce@assays$RNA@data[c("Tnf","Ifng","Gzmb","Gzma","Havcr2"),])))) %>% 
      reshape2::melt(.,id.vars=c("cellID","orig.ident","celltype","celltype1","seurat_clusters" )) %>% 
      data.table::setnames(.,old = c("variable","value"),new = c("gene","exp")) %>% 
      mutate(orig.ident=factor(orig.ident,levels = c("Vehicle","Sofa","R848","SR3"))) %>% 
      mutate(celltype=plyr::mapvalues(seurat_clusters,from=0:1,to=c("M2","M1")))
    #1倾向于M1，0倾向于M2
    p <- mar_metadata %>% 
      group_by(orig.ident, celltype, gene) %>% 
      dplyr::summarise(mean=mean(exp)) %>% 
      right_join(.,mar_metadata,by=c("orig.ident","gene","celltype")) %>% 
      ggplot(mapping=aes(x=orig.ident,y=exp,fill=mean)) +
      geom_violin(scale = "width", adjust = 2, trim = TRUE) +
      scale_fill_gradient(low = "#F2E80D", high = "#E65224") +
      ggsignif::geom_signif(comparisons = GetComb(c("Vehicle","Sofa","R848","SR3")), step_increase = 0.1, method=t.test) +
      facet_grid(gene~celltype,scales = "free") +
      theme_cowplot(font_size = 12) +
      theme(panel.spacing = unit(0, "lines"),
            panel.background = element_rect(fill = NA, color = "black"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            strip.text.y.left = element_text(angle = 0),
            axis.text.x = element_text(angle=60, vjust=1,hjust=1)) +
      coord_flip()   
    ggsave(p, filename = "./code/tianjin_Supple/singlecell_mouse/pictures/marcrophage.violin.plot.pdf", width = 12, height=14)
    
    #计算ifng基因在M1和M2中任意两组的差异
    mar_metadata %>% 
      group_by(orig.ident, celltype, gene) %>% 
      dplyr::summarise(mean=mean(exp)) %>% 
      right_join(.,mar_metadata,by=c("orig.ident","gene","celltype")) -> diff.data
    head(diff.data)
    
    diff.data %>% 
      filter(.,gene=="Ifng",orig.ident %in% c("R848","SR3"),celltype=="M2") %>% 
      as.data.frame() %>% 
      t.test(exp~orig.ident,data=.)            
      #rstatix::t_test(exp~orig.ident)
    
    diff.data %>% 
      filter(.,gene=="Ifng",orig.ident %in% c("Vehicle","Sofa"),celltype=="M1") %>% 
      as.data.frame() %>% 
      t.test(exp~orig.ident,data=.) 
    
    
    mar_metadata  
    signature <- format_msigdb(gmt="/Users/biofly/project/shijian/data/tianjian_ICC/h.all.v7.5.symbols.gmt", ont = "term", gene = "gene")
    kegg.signature <- format_msigdb(gmt="/Users/biofly/project/shijian/data/tianjian_ICC/c2.all.v7.5.1.symbols.gmt", ont = "term", gene = "gene")
    go.signature <- format_msigdb(gmt="/Users/biofly/project/shijian/data/tianjian_ICC/c5.all.v7.5.1.symbols.gmt", ont = "term", gene = "gene")
    eset <- Marcrophage_sce@assays$RNA@data
    signature <- lapply(signature,function(x){stringr::str_to_title(x)})
    kegg.signature <- lapply(kegg.signature,function(x){stringr::str_to_title(x)})
    go.signature <- lapply(go.signature,function(x){stringr::str_to_title(x)})
    eset[1:4,1:4]
    signature.score <- calculate_sig_score(pdata = NULL,
                                           eset=as.matrix(eset),
                                           signature =signature,
                                           method = "zscore",
                                           mini_gene_count = 3,
                                           column_of_sample = "ID",
                                           print_gene_propotion = FALSE,
                                           adjust_eset = FALSE,
                                           print_filtered_signatures = FALSE)
    kegg.score <- calculate_sig_score(pdata = NULL,
                                      eset=as.matrix(eset),
                                      signature = kegg.signature,
                                      method = "zscore",
                                      mini_gene_count = 3,
                                      column_of_sample = "ID",
                                      print_gene_propotion = FALSE,
                                      adjust_eset = FALSE,
                                      print_filtered_signatures = FALSE)
    go.score <- calculate_sig_score(pdata = NULL,
                                    eset=as.matrix(eset),
                                    signature = go.signature,
                                    method = "zscore",
                                    mini_gene_count = 3,
                                    column_of_sample = "ID",
                                    print_gene_propotion = FALSE,
                                    adjust_eset = FALSE,
                                    print_filtered_signatures = FALSE)
    hallmark.ssgsea.matrix <- signature.score %>% dplyr::select(.,-1) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    annotation_col <- Marcrophage_sce@meta.data 
    table(annotation_col$orig.ident)
    vehicle.R848.sig.res <- ScreenGenes(hallmark.ssgsea.matrix, design = annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Vehicle","R848"),
                                        contrast.control = c("Vehicle"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    vehicle.Sofa.sig.res <- ScreenGenes(hallmark.ssgsea.matrix, design = annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Vehicle","Sofa"),
                                        contrast.control = c("Vehicle"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    vehicle.SR3.sig.res <- ScreenGenes(hallmark.ssgsea.matrix, design = annotation_col, 
                                       contrast.col = "orig.ident", 
                                       contrast.level = c("Vehicle","SR3"),
                                       contrast.control = c("Vehicle"),
                                       method = c("t.test","wilcox.test")[1],
                                       cutoff.p=0.05,
                                       cutoff.q=0.05,
                                       cutoff.logFC=1,
                                       pmethod = "BH",
                                       verbose = F)
    Sofa.SR3.sig.res <- ScreenGenes(hallmark.ssgsea.matrix, design = annotation_col, 
                                    contrast.col = "orig.ident", 
                                    contrast.level = c("Sofa","SR3"),
                                    contrast.control = c("Sofa"),
                                    method = c("t.test","wilcox.test")[1],
                                    cutoff.p=0.05,
                                    cutoff.q=0.05,
                                    cutoff.logFC=1,
                                    pmethod = "BH",
                                    verbose = F)
    R848.SR3.sig.res <- ScreenGenes(hallmark.ssgsea.matrix, design = annotation_col, 
                                    contrast.col = "orig.ident", 
                                    contrast.level = c("R848","SR3"),
                                    contrast.control = c("R848"),
                                    method = c("t.test","wilcox.test")[1],
                                    cutoff.p=0.05,
                                    cutoff.q=0.05,
                                    cutoff.logFC=1,
                                    pmethod = "BH",
                                    verbose = F)
    Sofa.R848.sig.res <- ScreenGenes(hallmark.ssgsea.matrix, design = annotation_col, 
                                     contrast.col = "orig.ident", 
                                     contrast.level = c("Sofa","R848"),
                                     contrast.control = c("Sofa"),
                                     method = c("t.test","wilcox.test")[1],
                                     cutoff.p=0.05,
                                     cutoff.q=0.05,
                                     cutoff.logFC=1,
                                     pmethod = "BH",
                                     verbose = F)
    vehicle.R848.sig.res <- vehicle.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    vehicle.Sofa.sig.res <- vehicle.Sofa.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    vehicle.SR3.sig.res <- vehicle.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    Sofa.SR3.sig.res <- Sofa.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    R848.SR3.sig.res <- R848.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    Sofa.R848.sig.res <- Sofa.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    sig.hallmark <- unique(c(rownames(vehicle.R848.sig.res[vehicle.R848.sig.res$lab != "NO",]),rownames(vehicle.Sofa.sig.res[vehicle.Sofa.sig.res$lab != "NO",]),rownames(vehicle.SR3.sig.res[vehicle.SR3.sig.res$lab != "NO",])
                             ,rownames(Sofa.SR3.sig.res[Sofa.SR3.sig.res$lab != "NO",]),rownames(R848.SR3.sig.res[R848.SR3.sig.res$lab != "NO",]),rownames(Sofa.R848.sig.res[Sofa.R848.sig.res$lab != "NO",])))
    plotHeat(hallmark.ssgsea.matrix[sig.hallmark,], 
             annotation_col=annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    # col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    # col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols <- colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col <- circlize::colorRamp2(c(-5,0,5), cols)
    sam_anno <- HeatmapAnnotation(orig.ident = annotation_col$orig.ident)
    p.up <- Heatmap(t(scale(t(hallmark.ssgsea.matrix[sig.hallmark,]))),show_column_names = F, column_split = annotation_col$orig.ident,col=col,top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),row_names_gp = gpar(fontsize = 6,fontface="bold")) 
    draw(p.up,annotation_legend_side = "left")
    
    pathway.ssgsea.matrix <- kegg.score %>% dplyr::select(.,-1) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    annotation_col <- Marcrophage_sce@meta.data 
    table(annotation_col$orig.ident)
    vehicle.R848.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","R848"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    vehicle.Sofa.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","Sofa"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    vehicle.SR3.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                           contrast.col = "orig.ident", 
                                           contrast.level = c("Vehicle","SR3"),
                                           contrast.control = c("Vehicle"),
                                           method = c("t.test","wilcox.test")[1],
                                           cutoff.p=0.05,
                                           cutoff.q=0.05,
                                           cutoff.logFC=1,
                                           pmethod = "BH",
                                           verbose = F)
    Sofa.SR3.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Sofa","SR3"),
                                        contrast.control = c("Sofa"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    R848.SR3.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("R848","SR3"),
                                        contrast.control = c("R848"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    Sofa.R848.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                         contrast.col = "orig.ident", 
                                         contrast.level = c("Sofa","R848"),
                                         contrast.control = c("Sofa"),
                                         method = c("t.test","wilcox.test")[1],
                                         cutoff.p=0.05,
                                         cutoff.q=0.05,
                                         cutoff.logFC=1,
                                         pmethod = "BH",
                                         verbose = F)
    vehicle.R848.pathway.res <- vehicle.R848.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    vehicle.Sofa.pathway.res <- vehicle.Sofa.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    vehicle.SR3.pathway.res <- vehicle.SR3.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    Sofa.SR3.pathway.res <- Sofa.SR3.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    R848.SR3.pathway.res <- R848.SR3.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    Sofa.R848.pathway.res <- Sofa.R848.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    up.sig.pathway <- unique(c(rownames(vehicle.R848.pathway.res[vehicle.R848.pathway.res$lab == "UP",]),rownames(vehicle.Sofa.pathway.res[vehicle.Sofa.pathway.res$lab == "UP",]),rownames(vehicle.SR3.pathway.res[vehicle.SR3.pathway.res$lab == "UP",])
                               ,rownames(Sofa.SR3.pathway.res[Sofa.SR3.pathway.res$lab == "UP",]),rownames(R848.SR3.pathway.res[R848.SR3.pathway.res$lab == "UP",]),rownames(Sofa.R848.pathway.res[Sofa.R848.pathway.res$lab == "UP",])))
    down.sig.pathway <- unique(c(rownames(vehicle.R848.pathway.res[vehicle.R848.pathway.res$lab == "DOWN",]),rownames(vehicle.Sofa.pathway.res[vehicle.Sofa.pathway.res$lab == "DOWN",]),rownames(vehicle.SR3.pathway.res[vehicle.SR3.pathway.res$lab == "DOWN",])
                                 ,rownames(Sofa.SR3.pathway.res[Sofa.SR3.pathway.res$lab == "DOWN",]),rownames(R848.SR3.pathway.res[R848.SR3.pathway.res$lab == "DOWN",]),rownames(Sofa.R848.pathway.res[Sofa.R848.pathway.res$lab == "DOWN",])))
    sig.pathway <- unique(c(rownames(vehicle.R848.pathway.res[vehicle.R848.pathway.res$lab != "NO",]),rownames(vehicle.Sofa.pathway.res[vehicle.Sofa.pathway.res$lab != "NO",]),rownames(vehicle.SR3.pathway.res[vehicle.SR3.pathway.res$lab != "NO",])
                            ,rownames(Sofa.SR3.pathway.res[Sofa.SR3.pathway.res$lab != "NO",]),rownames(R848.SR3.pathway.res[R848.SR3.pathway.res$lab != "NO",]),rownames(Sofa.R848.pathway.res[Sofa.R848.pathway.res$lab != "NO",])))
    plotHeat(pathway.ssgsea.matrix[sig.pathway,], 
             annotation_col=annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    #col_scale <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
    #col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols <- colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col <- circlize::colorRamp2(c(-5,0,5), cols)
    sam_anno <- HeatmapAnnotation(orig.ident = annotation_col$orig.ident)
    p <- Heatmap(t(scale(t(pathway.ssgsea.matrix[sig.pathway[1:183],]))),show_column_names = F, column_split = annotation_col$orig.ident,col=col,top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),row_names_gp = gpar(fontsize = 5,fontface="bold")) 
    draw(p,annotation_legend_side = "left")
    
    
    go.ssgsea.matrix <- go.score %>% dplyr::select(.,-1) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    annotation_col <- Marcrophage_sce@meta.data 
    table(annotation_col$orig.ident)
    vehicle.R848.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                       contrast.col = "orig.ident", 
                                       contrast.level = c("Vehicle","R848"),
                                       contrast.control = c("Vehicle"),
                                       method = c("t.test","wilcox.test")[1],
                                       cutoff.p=0.05,
                                       cutoff.q=0.05,
                                       cutoff.logFC=1,
                                       pmethod = "BH",
                                       verbose = F)
    vehicle.Sofa.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                       contrast.col = "orig.ident", 
                                       contrast.level = c("Vehicle","Sofa"),
                                       contrast.control = c("Vehicle"),
                                       method = c("t.test","wilcox.test")[1],
                                       cutoff.p=0.05,
                                       cutoff.q=0.05,
                                       cutoff.logFC=1,
                                       pmethod = "BH",
                                       verbose = F)
    vehicle.SR3.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                      contrast.col = "orig.ident", 
                                      contrast.level = c("Vehicle","SR3"),
                                      contrast.control = c("Vehicle"),
                                      method = c("t.test","wilcox.test")[1],
                                      cutoff.p=0.05,
                                      cutoff.q=0.05,
                                      cutoff.logFC=1,
                                      pmethod = "BH",
                                      verbose = F)
    Sofa.SR3.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                   contrast.col = "orig.ident", 
                                   contrast.level = c("Sofa","SR3"),
                                   contrast.control = c("Sofa"),
                                   method = c("t.test","wilcox.test")[1],
                                   cutoff.p=0.05,
                                   cutoff.q=0.05,
                                   cutoff.logFC=1,
                                   pmethod = "BH",
                                   verbose = F)
    R848.SR3.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                   contrast.col = "orig.ident", 
                                   contrast.level = c("R848","SR3"),
                                   contrast.control = c("R848"),
                                   method = c("t.test","wilcox.test")[1],
                                   cutoff.p=0.05,
                                   cutoff.q=0.05,
                                   cutoff.logFC=1,
                                   pmethod = "BH",
                                   verbose = F)
    Sofa.R848.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                    contrast.col = "orig.ident", 
                                    contrast.level = c("Sofa","R848"),
                                    contrast.control = c("Sofa"),
                                    method = c("t.test","wilcox.test")[1],
                                    cutoff.p=0.05,
                                    cutoff.q=0.05,
                                    cutoff.logFC=1,
                                    pmethod = "BH",
                                    verbose = F)
    vehicle.R848.go.res <- vehicle.R848.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    vehicle.Sofa.go.res <- vehicle.Sofa.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    vehicle.SR3.go.res <- vehicle.SR3.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    Sofa.SR3.go.res <- Sofa.SR3.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    R848.SR3.go.res <- R848.SR3.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    Sofa.R848.go.res <- Sofa.R848.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    sig.go <- unique(c(rownames(vehicle.R848.go.res[vehicle.R848.go.res$lab != "NO",]),rownames(vehicle.Sofa.go.res[vehicle.Sofa.go.res$lab != "NO",]),rownames(vehicle.SR3.go.res[vehicle.SR3.go.res$lab != "NO",])
                       ,rownames(Sofa.SR3.go.res[Sofa.SR3.go.res$lab != "NO",]),rownames(R848.SR3.go.res[R848.SR3.go.res$lab != "NO",]),rownames(Sofa.R848.go.res[Sofa.R848.go.res$lab != "NO",])))
    plotHeat(go.ssgsea.matrix[sig.go,], 
             annotation_col=annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols<-colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col=circlize::colorRamp2(c(-5,0,5), cols)
    sam_anno <- HeatmapAnnotation(orig.ident = annotation_col$orig.ident)
    p <- Heatmap(t(scale(t(go.ssgsea.matrix[sig.pathway,]))),column_split = annotation_col$orig.ident,col=col,top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),row_names_gp = gpar(fontsize = 6,fontface="bold")) 
    draw(p,annotation_legend_side = "left")
  }
  
  #T细胞GSVA分析
  {
    signature <- format_msigdb(gmt="/Users/biofly/project/shijian/data/tianjian_ICC/h.all.v7.5.symbols.gmt", ont = "term", gene = "gene")
    kegg.signature <- format_msigdb(gmt="/Users/biofly/project/shijian/data/tianjian_ICC/c2.all.v7.5.1.symbols.gmt", ont = "term", gene = "gene")
    go.signature <- format_msigdb(gmt="/Users/biofly/project/shijian/data/tianjian_ICC/c5.all.v7.5.1.symbols.gmt", ont = "term", gene = "gene")
    eset <- T_sce1@assays$RNA@data
    CD8_sce <- subset(T_sce1,celltype=="CD8+T")
    CD4_sce <- subset(T_sce1,celltype=="CD4+T")
    CD8_eset <- CD8_sce@assays$RNA@data
    CD4_eset <- CD4_sce@assays$RNA@data
    signature <- lapply(signature,function(x){stringr::str_to_title(x)})
    kegg.signature <- lapply(kegg.signature,function(x){stringr::str_to_title(x)})
    go.signature <- lapply(go.signature,function(x){stringr::str_to_title(x)})
    eset[1:4,1:4]
    signature.score <- calculate_sig_score(pdata = NULL,
                                           eset=eset,
                                           signature =signature,
                                           method = "zscore",
                                           mini_gene_count = 3,
                                           column_of_sample = "ID",
                                           print_gene_propotion = FALSE,
                                           adjust_eset = FALSE,
                                           print_filtered_signatures = FALSE)
    kegg.score <- calculate_sig_score(pdata = NULL,
                                      eset=eset,
                                      signature = kegg.signature,
                                      method = "zscore",
                                      mini_gene_count = 3,
                                      column_of_sample = "ID",
                                      print_gene_propotion = FALSE,
                                      adjust_eset = FALSE,
                                      print_filtered_signatures = FALSE)
    go.score <- calculate_sig_score(pdata = NULL,
                                    eset=eset,
                                    signature = go.signature,
                                    method = "zscore",
                                    mini_gene_count = 3,
                                    column_of_sample = "ID",
                                    print_gene_propotion = FALSE,
                                    adjust_eset = FALSE,
                                    print_filtered_signatures = FALSE)
    hallmark.ssgsea.matrix <- signature.score %>% dplyr::select(.,-1) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    match(colnames(hallmark.ssgsea.matrix),rownames(T_sce1@meta.data))
    annotation_col <- T_sce1@meta.data 
    table(annotation_col$orig.ident)
    pos <- Fastgrep2(c("Vehicle","R848"),T_sce1@meta.data$orig.ident)  
    vehicle.R848.sig.res <- ScreenGenes(hallmark.ssgsea.matrix[,unlist(pos)], design = annotation_col[unlist(pos),], 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Vehicle","R848"),
                                        contrast.control = c("Vehicle"),
                                        method = c("t.test","wilcox.test")[2],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = T)
    pos <- Fastgrep2(c("Vehicle","Sofa"),T_sce1@meta.data$orig.ident)  
    vehicle.Sofa.sig.res <- ScreenGenes(hallmark.ssgsea.matrix[,unlist(pos)], design = annotation_col[unlist(pos),], 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Vehicle","Sofa"),
                                        contrast.control = c("Vehicle"),
                                        method = c("t.test","wilcox.test")[2],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = T)
    pos <- Fastgrep2(c("Vehicle","SR3"),T_sce1@meta.data$orig.ident)  
    vehicle.SR3.sig.res <- ScreenGenes(hallmark.ssgsea.matrix[,unlist(pos)], design = annotation_col[unlist(pos),], 
                                       contrast.col = "orig.ident", 
                                       contrast.level = c("Vehicle","SR3"),
                                       contrast.control = c("Vehicle"),
                                       method = c("t.test","wilcox.test")[2],
                                       cutoff.p=0.05,
                                       cutoff.q=0.05,
                                       cutoff.logFC=1,
                                       pmethod = "BH",
                                       verbose = T)
    pos <- Fastgrep2(c("Sofa","SR3"),T_sce1@meta.data$orig.ident)
    Sofa.SR3.sig.res <- ScreenGenes(hallmark.ssgsea.matrix[,unlist(pos)], design = annotation_col[unlist(pos),], 
                                    contrast.col = "orig.ident", 
                                    contrast.level = c("Sofa","SR3"),
                                    contrast.control = c("Sofa"),
                                    method = c("t.test","wilcox.test")[2],
                                    cutoff.p=0.05,
                                    cutoff.q=0.05,
                                    cutoff.logFC=1,
                                    pmethod = "BH",
                                    verbose = T)
    pos <- Fastgrep2(c("R848","SR3"),T_sce1@meta.data$orig.ident)
    R848.SR3.sig.res <- ScreenGenes(hallmark.ssgsea.matrix[,unlist(pos)], design = annotation_col[unlist(pos),], 
                                    contrast.col = "orig.ident", 
                                    contrast.level = c("R848","SR3"),
                                    contrast.control = c("R848"),
                                    method = c("t.test","wilcox.test")[2],
                                    cutoff.p=0.05,
                                    cutoff.q=0.05,
                                    cutoff.logFC=1,
                                    pmethod = "BH",
                                    verbose = T)
    pos <- Fastgrep2(c("Sofa","R848"),T_sce1@meta.data$orig.ident)
    Sofa.R848.sig.res <- ScreenGenes(hallmark.ssgsea.matrix[,unlist(pos)], design = annotation_col[unlist(pos),], 
                                     contrast.col = "orig.ident", 
                                     contrast.level = c("Sofa","R848"),
                                     contrast.control = c("Sofa"),
                                     method = c("t.test","wilcox.test")[2],
                                     cutoff.p=0.05,
                                     cutoff.q=0.05,
                                     cutoff.logFC=1,
                                     pmethod = "BH",
                                     verbose = T)
    vehicle.R848.sig.res <- vehicle.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    vehicle.Sofa.sig.res <- vehicle.Sofa.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    vehicle.SR3.sig.res <- vehicle.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    Sofa.SR3.sig.res <- Sofa.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    R848.SR3.sig.res <- R848.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    Sofa.R848.sig.res <- Sofa.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    sig.hallmark <- unique(c(rownames(vehicle.R848.sig.res[vehicle.R848.sig.res$lab != "NO",]),rownames(vehicle.Sofa.sig.res[vehicle.Sofa.sig.res$lab != "NO",]),rownames(vehicle.SR3.sig.res[vehicle.SR3.sig.res$lab != "NO",])
                             ,rownames(Sofa.SR3.sig.res[Sofa.SR3.sig.res$lab != "NO",]),rownames(R848.SR3.sig.res[R848.SR3.sig.res$lab != "NO",]),rownames(Sofa.R848.sig.res[Sofa.R848.sig.res$lab != "NO",])))
    plotHeat(hallmark.ssgsea.matrix[sig.hallmark,], 
             annotation_col=annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    # col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    # col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols <- colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col <- circlize::colorRamp2(c(-5,0,5), cols)
    sam_anno <- HeatmapAnnotation(orig.ident = annotation_col$orig.ident)
    p <- Heatmap(t(scale(t(hallmark.ssgsea.matrix[sig.hallmark,]))),
                 show_column_names = F, 
                 column_split = annotation_col$orig.ident,
                 col=col,
                 top_annotation=sam_anno,
                 row_names_max_width = unit(25, "cm"),
                 row_names_gp = gpar(fontsize = 6,fontface="bold")) 
    draw(p, annotation_legend_side = "left")
    
    pathway.ssgsea.matrix <- kegg.score %>% dplyr::select(.,-1) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    annotation_col <- T_sce1@meta.data 
    table(annotation_col$orig.ident)
    vehicle.R848.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","R848"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    vehicle.Sofa.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","Sofa"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    vehicle.SR3.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                           contrast.col = "orig.ident", 
                                           contrast.level = c("Vehicle","SR3"),
                                           contrast.control = c("Vehicle"),
                                           method = c("t.test","wilcox.test")[1],
                                           cutoff.p=0.05,
                                           cutoff.q=0.05,
                                           cutoff.logFC=1,
                                           pmethod = "BH",
                                           verbose = F)
    Sofa.SR3.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Sofa","SR3"),
                                        contrast.control = c("Sofa"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    R848.SR3.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("R848","SR3"),
                                        contrast.control = c("R848"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    Sofa.R848.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                         contrast.col = "orig.ident", 
                                         contrast.level = c("Sofa","R848"),
                                         contrast.control = c("Sofa"),
                                         method = c("t.test","wilcox.test")[1],
                                         cutoff.p=0.05,
                                         cutoff.q=0.05,
                                         cutoff.logFC=1,
                                         pmethod = "BH",
                                         verbose = F)
    vehicle.R848.pathway.res <- vehicle.R848.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    vehicle.Sofa.pathway.res <- vehicle.Sofa.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    vehicle.SR3.pathway.res <- vehicle.SR3.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    Sofa.SR3.pathway.res <- Sofa.SR3.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    R848.SR3.pathway.res <- R848.SR3.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    Sofa.R848.pathway.res <- Sofa.R848.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    up.sig.pathway <- unique(c(rownames(vehicle.R848.pathway.res[vehicle.R848.pathway.res$lab == "UP",]),rownames(vehicle.Sofa.pathway.res[vehicle.Sofa.pathway.res$lab == "UP",]),rownames(vehicle.SR3.pathway.res[vehicle.SR3.pathway.res$lab == "UP",])
                               ,rownames(Sofa.SR3.pathway.res[Sofa.SR3.pathway.res$lab == "UP",]),rownames(R848.SR3.pathway.res[R848.SR3.pathway.res$lab == "UP",]),rownames(Sofa.R848.pathway.res[Sofa.R848.pathway.res$lab == "UP",])))
    down.sig.pathway <- unique(c(rownames(vehicle.R848.pathway.res[vehicle.R848.pathway.res$lab == "DOWN",]),rownames(vehicle.Sofa.pathway.res[vehicle.Sofa.pathway.res$lab == "DOWN",]),rownames(vehicle.SR3.pathway.res[vehicle.SR3.pathway.res$lab == "DOWN",])
                                 ,rownames(Sofa.SR3.pathway.res[Sofa.SR3.pathway.res$lab == "DOWN",]),rownames(R848.SR3.pathway.res[R848.SR3.pathway.res$lab == "DOWN",]),rownames(Sofa.R848.pathway.res[Sofa.R848.pathway.res$lab == "DOWN",])))
    sig.pathway <- unique(c(rownames(vehicle.R848.pathway.res[vehicle.R848.pathway.res$lab != "NO",]),rownames(vehicle.Sofa.pathway.res[vehicle.Sofa.pathway.res$lab != "NO",]),rownames(vehicle.SR3.pathway.res[vehicle.SR3.pathway.res$lab != "NO",])
                            ,rownames(Sofa.SR3.pathway.res[Sofa.SR3.pathway.res$lab != "NO",]),rownames(R848.SR3.pathway.res[R848.SR3.pathway.res$lab != "NO",]),rownames(Sofa.R848.pathway.res[Sofa.R848.pathway.res$lab != "NO",])))
    plotHeat(pathway.ssgsea.matrix[sig.pathway,], 
             annotation_col=annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    # col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    # col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols<-colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col_palette = circlize::colorRamp2(c(-5,0,5), cols)
    sam_anno <- HeatmapAnnotation(orig.ident = annotation_col$orig.ident)
    p <- Heatmap(t(scale(t(pathway.ssgsea.matrix[sig.pathway,]))),show_column_names = F, column_split = annotation_col$orig.ident,col=col_palette,top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),row_names_gp = gpar(fontsize = 5,fontface="bold")) 
    draw(p,annotation_legend_side = "left")
    
    
    go.ssgsea.matrix <- go.score %>% dplyr::select(.,-2) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    annotation_col <- T_sce@meta.data 
    table(annotation_col$orig.ident)
    vehicle.R848.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                       contrast.col = "orig.ident", 
                                       contrast.level = c("Vehicle","R848"),
                                       contrast.control = c("Vehicle"),
                                       method = c("t.test","wilcox.test")[1],
                                       cutoff.p=0.05,
                                       cutoff.q=0.05,
                                       cutoff.logFC=1,
                                       pmethod = "BH",
                                       verbose = F)
    vehicle.Sofa.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                       contrast.col = "orig.ident", 
                                       contrast.level = c("Vehicle","Sofa"),
                                       contrast.control = c("Vehicle"),
                                       method = c("t.test","wilcox.test")[1],
                                       cutoff.p=0.05,
                                       cutoff.q=0.05,
                                       cutoff.logFC=1,
                                       pmethod = "BH",
                                       verbose = F)
    vehicle.SR3.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                      contrast.col = "orig.ident", 
                                      contrast.level = c("Vehicle","SR3"),
                                      contrast.control = c("Vehicle"),
                                      method = c("t.test","wilcox.test")[1],
                                      cutoff.p=0.05,
                                      cutoff.q=0.05,
                                      cutoff.logFC=1,
                                      pmethod = "BH",
                                      verbose = F)
    Sofa.SR3.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                   contrast.col = "orig.ident", 
                                   contrast.level = c("Sofa","SR3"),
                                   contrast.control = c("Sofa"),
                                   method = c("t.test","wilcox.test")[1],
                                   cutoff.p=0.05,
                                   cutoff.q=0.05,
                                   cutoff.logFC=1,
                                   pmethod = "BH",
                                   verbose = F)
    R848.SR3.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                   contrast.col = "orig.ident", 
                                   contrast.level = c("R848","SR3"),
                                   contrast.control = c("R848"),
                                   method = c("t.test","wilcox.test")[1],
                                   cutoff.p=0.05,
                                   cutoff.q=0.05,
                                   cutoff.logFC=1,
                                   pmethod = "BH",
                                   verbose = F)
    Sofa.R848.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                    contrast.col = "orig.ident", 
                                    contrast.level = c("Sofa","R848"),
                                    contrast.control = c("Sofa"),
                                    method = c("t.test","wilcox.test")[1],
                                    cutoff.p=0.05,
                                    cutoff.q=0.05,
                                    cutoff.logFC=1,
                                    pmethod = "BH",
                                    verbose = F)
    vehicle.R848.go.res <- vehicle.R848.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    vehicle.Sofa.go.res <- vehicle.Sofa.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    vehicle.SR3.go.res <- vehicle.SR3.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    Sofa.SR3.go.res <- Sofa.SR3.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    R848.SR3.go.res <- R848.SR3.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    Sofa.R848.go.res <- Sofa.R848.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    sig.go <- unique(c(rownames(vehicle.R848.go.res[vehicle.R848.go.res$lab != "NO",]),rownames(vehicle.Sofa.go.res[vehicle.Sofa.go.res$lab != "NO",]),rownames(vehicle.SR3.go.res[vehicle.SR3.go.res$lab != "NO",])
                       ,rownames(Sofa.SR3.go.res[Sofa.SR3.go.res$lab != "NO",]),rownames(R848.SR3.go.res[R848.SR3.go.res$lab != "NO",]),rownames(Sofa.R848.go.res[Sofa.R848.go.res$lab != "NO",])))
    plotHeat(go.ssgsea.matrix[sig.go,], 
             annotation_col=annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    sam_anno <- HeatmapAnnotation(orig.ident = annotation_col$orig.ident)
    p <- Heatmap(t(scale(t(go.ssgsea.matrix[sig.pathway,]))),column_split = annotation_col$orig.ident,col=col_palette,top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),row_names_gp = gpar(fontsize = 6,fontface="bold")) 
    draw(p,annotation_legend_side = "left")
  }
  
  #infercnv
  {
    class(sce_merge_res1$sce@meta.data$celltype)
    table(sce_merge_res1$sce@meta.data$celltype)
    sce_merge_res1$sce@meta.data$celltypes <- NA
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==0] ="EPCAM_cell_0"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==1] = "CD8+T_cell_1"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==2] = "Neutrophil_cell_2"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==3] = "Neutrophil_cell_3"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==4] = "Marcrophage_4"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==5] = "CD8+T_cell_5"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==6] = "NK_6"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==7] = "CD4+T_cell_7"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==8] = "Stromal_cell_8"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==9] = "Dendritic_cell_9"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==10] = "EPCAM_cell_10"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==11] = "CD8+T_cell_11"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==12] = "B_cell_12"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==13] = "Dendritic_cell_13"
    sce_merge_res1$sce@meta.data$celltypes[sce_merge_res1$sce@meta.data$seurat_clusters==14] = "Mast_cell_14"
    res_mouse <- FastInferCNV(obj=sce_merge_res1$sce,
                              counts=NULL,
                              annotations_file = NULL,
                              cellType = "celltypes",
                              gene_order_file = NULL,
                              ref_group_names = c("CD8+T_cell_1","Neutrophil_cell_2","Neutrophil_cell_3","Marcrophage_4","CD8+T_cell_5","NK_6","CD4+T_cell_7","Stromal_cell_8","Dendritic_cell_9","CD8+T_cell_11","B_cell_12","Dendritic_cell_13","Mast_cell_14"),  
                              cutoff = c(0.1,1)[1], 
                              analysis_mode = c("samples", "subclusters", "cells")[1],
                              out_dir = "mouse_inferCNV",
                              cluster = T,
                              HMM = F,
                              denoise = T,
                              anno.db = "/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/gencode.vM28.annotation.gtf")
  }
  
  #提取肿瘤细胞(0，5,11都是肿瘤细胞？还是0，11是肿瘤细胞？)做GSVA分析
  {
    Epithelial_sce <- subset(sce_merge_res1$sce, seurat_clusters==0 | seurat_clusters==11)
    signature <- format_msigdb(gmt="/Users/biofly/project/shijian/data/tianjian_ICC/h.all.v7.5.symbols.gmt", ont = "term", gene = "gene")
    kegg.signature <- format_msigdb(gmt="/Users/biofly/project/shijian/data/tianjian_ICC/c2.all.v7.5.1.symbols.gmt", ont = "term", gene = "gene")
    go.signature <- format_msigdb(gmt="/Users/biofly/project/shijian/data/tianjian_ICC/c5.all.v7.5.1.symbols.gmt", ont = "term", gene = "gene")
    eset <- Epithelial_sce@assays$RNA@data
    signature <- lapply(signature,function(x){stringr::str_to_title(x)})
    kegg.signature <- lapply(kegg.signature,function(x){stringr::str_to_title(x)})
    go.signature <- lapply(go.signature,function(x){stringr::str_to_title(x)})
    eset[1:4,1:4]
    signature.score <- calculate_sig_score(pdata = NULL,
                                           eset=eset,
                                           signature =signature,
                                           method = "zscore",
                                           mini_gene_count = 3,
                                           column_of_sample = "ID",
                                           print_gene_propotion = FALSE,
                                           adjust_eset = FALSE,
                                           print_filtered_signatures = FALSE)
    kegg.score <- calculate_sig_score(pdata = NULL,
                                      eset=eset,
                                      signature = kegg.signature,
                                      method = "zscore",
                                      mini_gene_count = 3,
                                      column_of_sample = "ID",
                                      print_gene_propotion = FALSE,
                                      adjust_eset = FALSE,
                                      print_filtered_signatures = FALSE)
    go.score <- calculate_sig_score(pdata = NULL,
                                    eset=eset,
                                    signature = go.signature,
                                    method = "zscore",
                                    mini_gene_count = 3,
                                    column_of_sample = "ID",
                                    print_gene_propotion = FALSE,
                                    adjust_eset = FALSE,
                                    print_filtered_signatures = FALSE)
    hallmark.ssgsea.matrix <- signature.score %>% dplyr::select(.,-1) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    annotation_col <- Epithelial_sce@meta.data 
    table(annotation_col$orig.ident)
    vehicle.R848.sig.res <- ScreenGenes(hallmark.ssgsea.matrix, design = annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Vehicle","R848"),
                                        contrast.control = c("Vehicle"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    vehicle.Sofa.sig.res <- ScreenGenes(hallmark.ssgsea.matrix, design = annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Vehicle","Sofa"),
                                        contrast.control = c("Vehicle"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    vehicle.SR3.sig.res <- ScreenGenes(hallmark.ssgsea.matrix, design = annotation_col, 
                                       contrast.col = "orig.ident", 
                                       contrast.level = c("Vehicle","SR3"),
                                       contrast.control = c("Vehicle"),
                                       method = c("t.test","wilcox.test")[1],
                                       cutoff.p=0.05,
                                       cutoff.q=0.05,
                                       cutoff.logFC=1,
                                       pmethod = "BH",
                                       verbose = F)
    Sofa.SR3.sig.res <- ScreenGenes(hallmark.ssgsea.matrix, design = annotation_col, 
                                    contrast.col = "orig.ident", 
                                    contrast.level = c("Sofa","SR3"),
                                    contrast.control = c("Sofa"),
                                    method = c("t.test","wilcox.test")[1],
                                    cutoff.p=0.05,
                                    cutoff.q=0.05,
                                    cutoff.logFC=1,
                                    pmethod = "BH",
                                    verbose = F)
    R848.SR3.sig.res <- ScreenGenes(hallmark.ssgsea.matrix, design = annotation_col, 
                                    contrast.col = "orig.ident", 
                                    contrast.level = c("R848","SR3"),
                                    contrast.control = c("R848"),
                                    method = c("t.test","wilcox.test")[1],
                                    cutoff.p=0.05,
                                    cutoff.q=0.05,
                                    cutoff.logFC=1,
                                    pmethod = "BH",
                                    verbose = F)
    Sofa.R848.sig.res <- ScreenGenes(hallmark.ssgsea.matrix, design = annotation_col, 
                                     contrast.col = "orig.ident", 
                                     contrast.level = c("Sofa","R848"),
                                     contrast.control = c("Sofa"),
                                     method = c("t.test","wilcox.test")[1],
                                     cutoff.p=0.05,
                                     cutoff.q=0.05,
                                     cutoff.logFC=1,
                                     pmethod = "BH",
                                     verbose = F)
    vehicle.R848.sig.res <- vehicle.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    vehicle.Sofa.sig.res <- vehicle.Sofa.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    vehicle.SR3.sig.res <- vehicle.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    Sofa.SR3.sig.res <- Sofa.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    R848.SR3.sig.res <- R848.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    Sofa.R848.sig.res <- Sofa.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    sig.hallmark <- unique(c(rownames(vehicle.R848.sig.res[vehicle.R848.sig.res$lab != "NO",]),rownames(vehicle.Sofa.sig.res[vehicle.Sofa.sig.res$lab != "NO",]),rownames(vehicle.SR3.sig.res[vehicle.SR3.sig.res$lab != "NO",])
                             ,rownames(Sofa.SR3.sig.res[Sofa.SR3.sig.res$lab != "NO",]),rownames(R848.SR3.sig.res[R848.SR3.sig.res$lab != "NO",]),rownames(Sofa.R848.sig.res[Sofa.R848.sig.res$lab != "NO",])))
    plotHeat(hallmark.ssgsea.matrix[sig.hallmark,], 
             annotation_col=annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols <- colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col <- circlize::colorRamp2(c(-5,0,5), cols)
    sam_anno <- HeatmapAnnotation(orig.ident = annotation_col$orig.ident)
    p <- Heatmap(t(scale(t(hallmark.ssgsea.matrix[sig.hallmark,]))),show_column_names = F, column_split = annotation_col$orig.ident,col=col,top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),row_names_gp = gpar(fontsize = 6,fontface="bold")) 
    draw(p, annotation_legend_side = "left")
    
    pathway.ssgsea.matrix <- kegg.score %>% dplyr::select(.,-1) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    annotation_col <- Epithelial_sce@meta.data 
    table(annotation_col$orig.ident)
    vehicle.R848.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","R848"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    vehicle.Sofa.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","Sofa"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    vehicle.SR3.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                           contrast.col = "orig.ident", 
                                           contrast.level = c("Vehicle","SR3"),
                                           contrast.control = c("Vehicle"),
                                           method = c("t.test","wilcox.test")[1],
                                           cutoff.p=0.05,
                                           cutoff.q=0.05,
                                           cutoff.logFC=1,
                                           pmethod = "BH",
                                           verbose = F)
    Sofa.SR3.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Sofa","SR3"),
                                        contrast.control = c("Sofa"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    R848.SR3.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("R848","SR3"),
                                        contrast.control = c("R848"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    Sofa.R848.pathway.res <- ScreenGenes(pathway.ssgsea.matrix, design = annotation_col, 
                                         contrast.col = "orig.ident", 
                                         contrast.level = c("Sofa","R848"),
                                         contrast.control = c("Sofa"),
                                         method = c("t.test","wilcox.test")[1],
                                         cutoff.p=0.05,
                                         cutoff.q=0.05,
                                         cutoff.logFC=1,
                                         pmethod = "BH",
                                         verbose = F)
    vehicle.R848.pathway.res <- vehicle.R848.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    vehicle.Sofa.pathway.res <- vehicle.Sofa.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    vehicle.SR3.pathway.res <- vehicle.SR3.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    Sofa.SR3.pathway.res <- Sofa.SR3.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    R848.SR3.pathway.res <- R848.SR3.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    Sofa.R848.pathway.res <- Sofa.R848.pathway.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    up.sig.pathway <- unique(c(rownames(vehicle.R848.pathway.res[vehicle.R848.pathway.res$lab == "UP",]),rownames(vehicle.Sofa.pathway.res[vehicle.Sofa.pathway.res$lab == "UP",]),rownames(vehicle.SR3.pathway.res[vehicle.SR3.pathway.res$lab == "UP",])
                               ,rownames(Sofa.SR3.pathway.res[Sofa.SR3.pathway.res$lab == "UP",]),rownames(R848.SR3.pathway.res[R848.SR3.pathway.res$lab == "UP",]),rownames(Sofa.R848.pathway.res[Sofa.R848.pathway.res$lab == "UP",])))
    down.sig.pathway <- unique(c(rownames(vehicle.R848.pathway.res[vehicle.R848.pathway.res$lab == "DOWN",]),rownames(vehicle.Sofa.pathway.res[vehicle.Sofa.pathway.res$lab == "DOWN",]),rownames(vehicle.SR3.pathway.res[vehicle.SR3.pathway.res$lab == "DOWN",])
                                 ,rownames(Sofa.SR3.pathway.res[Sofa.SR3.pathway.res$lab == "DOWN",]),rownames(R848.SR3.pathway.res[R848.SR3.pathway.res$lab == "DOWN",]),rownames(Sofa.R848.pathway.res[Sofa.R848.pathway.res$lab == "DOWN",])))
    sig.pathway <- unique(c(rownames(vehicle.R848.pathway.res[vehicle.R848.pathway.res$lab != "NO",]),rownames(vehicle.Sofa.pathway.res[vehicle.Sofa.pathway.res$lab != "NO",]),rownames(vehicle.SR3.pathway.res[vehicle.SR3.pathway.res$lab != "NO",])
                            ,rownames(Sofa.SR3.pathway.res[Sofa.SR3.pathway.res$lab != "NO",]),rownames(R848.SR3.pathway.res[R848.SR3.pathway.res$lab != "NO",]),rownames(Sofa.R848.pathway.res[Sofa.R848.pathway.res$lab != "NO",])))
    plotHeat(pathway.ssgsea.matrix[sig.pathway,], 
             annotation_col=annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    #col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    #col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols<-colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col_palette <- circlize::colorRamp2(c(-5,0,5), cols) #颜色插值
    sam_anno <- HeatmapAnnotation(orig.ident = annotation_col$orig.ident)
    p <- Heatmap(t(scale(t(pathway.ssgsea.matrix[sig.pathway,]))),
                 show_column_names = F, column_split = annotation_col$orig.ident,col=col_palette,
                 top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),
                 row_names_gp = gpar(fontsize = 5,fontface="bold"),use_raster = F) 
    draw(p,annotation_legend_side = "left")
    
    
    go.ssgsea.matrix <- go.score %>% dplyr::select(.,-2) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    annotation_col <- Epithelial_sce@meta.data 
    table(annotation_col$orig.ident)
    vehicle.R848.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                       contrast.col = "orig.ident", 
                                       contrast.level = c("Vehicle","R848"),
                                       contrast.control = c("Vehicle"),
                                       method = c("t.test","wilcox.test")[1],
                                       cutoff.p=0.05,
                                       cutoff.q=0.05,
                                       cutoff.logFC=1,
                                       pmethod = "BH",
                                       verbose = F)
    vehicle.Sofa.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                       contrast.col = "orig.ident", 
                                       contrast.level = c("Vehicle","Sofa"),
                                       contrast.control = c("Vehicle"),
                                       method = c("t.test","wilcox.test")[1],
                                       cutoff.p=0.05,
                                       cutoff.q=0.05,
                                       cutoff.logFC=1,
                                       pmethod = "BH",
                                       verbose = F)
    vehicle.SR3.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                      contrast.col = "orig.ident", 
                                      contrast.level = c("Vehicle","SR3"),
                                      contrast.control = c("Vehicle"),
                                      method = c("t.test","wilcox.test")[1],
                                      cutoff.p=0.05,
                                      cutoff.q=0.05,
                                      cutoff.logFC=1,
                                      pmethod = "BH",
                                      verbose = F)
    Sofa.SR3.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                   contrast.col = "orig.ident", 
                                   contrast.level = c("Sofa","SR3"),
                                   contrast.control = c("Sofa"),
                                   method = c("t.test","wilcox.test")[1],
                                   cutoff.p=0.05,
                                   cutoff.q=0.05,
                                   cutoff.logFC=1,
                                   pmethod = "BH",
                                   verbose = F)
    R848.SR3.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                   contrast.col = "orig.ident", 
                                   contrast.level = c("R848","SR3"),
                                   contrast.control = c("R848"),
                                   method = c("t.test","wilcox.test")[1],
                                   cutoff.p=0.05,
                                   cutoff.q=0.05,
                                   cutoff.logFC=1,
                                   pmethod = "BH",
                                   verbose = F)
    Sofa.R848.go.res <- ScreenGenes(go.ssgsea.matrix, design = annotation_col, 
                                    contrast.col = "orig.ident", 
                                    contrast.level = c("Sofa","R848"),
                                    contrast.control = c("Sofa"),
                                    method = c("t.test","wilcox.test")[1],
                                    cutoff.p=0.05,
                                    cutoff.q=0.05,
                                    cutoff.logFC=1,
                                    pmethod = "BH",
                                    verbose = F)
    vehicle.R848.go.res <- vehicle.R848.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    vehicle.Sofa.go.res <- vehicle.Sofa.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    vehicle.SR3.go.res <- vehicle.SR3.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    Sofa.SR3.go.res <- Sofa.SR3.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    R848.SR3.go.res <- R848.SR3.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    Sofa.R848.go.res <- Sofa.R848.go.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    sig.go <- unique(c(rownames(vehicle.R848.go.res[vehicle.R848.go.res$lab != "NO",]),rownames(vehicle.Sofa.go.res[vehicle.Sofa.go.res$lab != "NO",]),rownames(vehicle.SR3.go.res[vehicle.SR3.go.res$lab != "NO",])
                       ,rownames(Sofa.SR3.go.res[Sofa.SR3.go.res$lab != "NO",]),rownames(R848.SR3.go.res[R848.SR3.go.res$lab != "NO",]),rownames(Sofa.R848.go.res[Sofa.R848.go.res$lab != "NO",])))
    plotHeat(go.ssgsea.matrix[sig.go,], 
             annotation_col=annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    sam_anno <- HeatmapAnnotation(orig.ident = annotation_col$orig.ident)
    p <- Heatmap(t(scale(t(go.ssgsea.matrix[sig.pathway,]))),column_split = annotation_col$orig.ident,col=col_palette,top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),row_names_gp = gpar(fontsize = 6,fontface="bold")) 
    draw(p,annotation_legend_side = "left")
  }

  #CD8 CD4 hallmark
  {
    CD8_signature.score <- calculate_sig_score(pdata = NULL,
                                               eset=ass.matrix(CD8_eset),
                                               signature =signature,
                                               method = "zscore",
                                               mini_gene_count = 3,
                                               column_of_sample = "ID",
                                               print_gene_propotion = FALSE,
                                               adjust_eset = FALSE,
                                               print_filtered_signatures = FALSE)
    CD4_signature.score <- calculate_sig_score(pdata = NULL,
                                               eset=as.matrix(CD4_eset),
                                               signature =signature,
                                               method = "zscore",
                                               mini_gene_count = 3,
                                               column_of_sample = "ID",
                                               print_gene_propotion = FALSE,
                                               adjust_eset = FALSE,
                                               print_filtered_signatures = FALSE)
    CD8_hallmark.ssgsea.matrix <- CD8_signature.score %>% dplyr::select(.,-1) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    CD4_hallmark.ssgsea.matrix <- CD4_signature.score %>% dplyr::select(.,-1) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    CD8_annotation_col <- CD8_sce@meta.data 
    CD4_annotation_col <- CD4_sce@meta.data
    table(CD8_annotation_col$orig.ident)
    table(CD4_annotation_col$orig.ident)
    CD8_vehicle.R848.sig.res <- ScreenGenes(CD8_hallmark.ssgsea.matrix, design = CD8_annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","R848"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    CD8_vehicle.Sofa.sig.res <- ScreenGenes(CD8_hallmark.ssgsea.matrix, design = CD8_annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","Sofa"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    CD8_vehicle.SR3.sig.res <- ScreenGenes(CD8_hallmark.ssgsea.matrix, design = CD8_annotation_col, 
                                           contrast.col = "orig.ident", 
                                           contrast.level = c("Vehicle","SR3"),
                                           contrast.control = c("Vehicle"),
                                           method = c("t.test","wilcox.test")[1],
                                           cutoff.p=0.05,
                                           cutoff.q=0.05,
                                           cutoff.logFC=1,
                                           pmethod = "BH",
                                           verbose = F)
    CD8_Sofa.SR3.sig.res <- ScreenGenes(CD8_hallmark.ssgsea.matrix, design = CD8_annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Sofa","SR3"),
                                        contrast.control = c("Sofa"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    CD8_R848.SR3.sig.res <- ScreenGenes(CD8_hallmark.ssgsea.matrix, design = CD8_annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("R848","SR3"),
                                        contrast.control = c("R848"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    CD8_Sofa.R848.sig.res <- ScreenGenes(CD8_hallmark.ssgsea.matrix, design = CD8_annotation_col, 
                                         contrast.col = "orig.ident", 
                                         contrast.level = c("Sofa","R848"),
                                         contrast.control = c("Sofa"),
                                         method = c("t.test","wilcox.test")[1],
                                         cutoff.p=0.05,
                                         cutoff.q=0.05,
                                         cutoff.logFC=1,
                                         pmethod = "BH",
                                         verbose = F)
    CD8_vehicle.R848.sig.res <- CD8_vehicle.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    CD8_vehicle.Sofa.sig.res <- CD8_vehicle.Sofa.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    CD8_vehicle.SR3.sig.res <- CD8_vehicle.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    CD8_Sofa.SR3.sig.res <- CD8_Sofa.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    CD8_R848.SR3.sig.res <- CD8_R848.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    CD8_Sofa.R848.sig.res <- CD8_Sofa.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    CD8_sig.hallmark <- unique(c(rownames(CD8_vehicle.R848.sig.res[CD8_vehicle.R848.sig.res$lab != "NO",]),rownames(CD8_vehicle.Sofa.sig.res[CD8_vehicle.Sofa.sig.res$lab != "NO",]),rownames(CD8_vehicle.SR3.sig.res[CD8_vehicle.SR3.sig.res$lab != "NO",])
                                 ,rownames(CD8_Sofa.SR3.sig.res[CD8_Sofa.SR3.sig.res$lab != "NO",]),rownames(CD8_R848.SR3.sig.res[CD8_R848.SR3.sig.res$lab != "NO",]),rownames(CD8_Sofa.R848.sig.res[CD8_Sofa.R848.sig.res$lab != "NO",])))
    plotHeat(CD8_hallmark.ssgsea.matrix[CD8_sig.hallmark,], 
             annotation_col=annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols <- colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col <- circlize::colorRamp2(c(-5,0,5), cols)
    sam_anno <- HeatmapAnnotation(orig.ident = CD8_annotation_col$orig.ident)
    p <- Heatmap(t(scale(t(CD8_hallmark.ssgsea.matrix[CD8_sig.hallmark,]))),show_column_names = F, column_split = CD8_annotation_col$orig.ident,col=col,top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),row_names_gp = gpar(fontsize = 6,fontface="bold")) 
    draw(p, annotation_legend_side = "left")
    
    CD4_vehicle.R848.sig.res <- ScreenGenes(CD4_hallmark.ssgsea.matrix, design = CD4_annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","R848"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    CD4_vehicle.Sofa.sig.res <- ScreenGenes(CD4_hallmark.ssgsea.matrix, design = CD4_annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","Sofa"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    CD4_vehicle.SR3.sig.res <- ScreenGenes(CD4_hallmark.ssgsea.matrix, design = CD4_annotation_col, 
                                           contrast.col = "orig.ident", 
                                           contrast.level = c("Vehicle","SR3"),
                                           contrast.control = c("Vehicle"),
                                           method = c("t.test","wilcox.test")[1],
                                           cutoff.p=0.05,
                                           cutoff.q=0.05,
                                           cutoff.logFC=1,
                                           pmethod = "BH",
                                           verbose = F)
    CD4_Sofa.SR3.sig.res <- ScreenGenes(CD4_hallmark.ssgsea.matrix, design = CD4_annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Sofa","SR3"),
                                        contrast.control = c("Sofa"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    CD4_R848.SR3.sig.res <- ScreenGenes(CD4_hallmark.ssgsea.matrix, design = CD4_annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("R848","SR3"),
                                        contrast.control = c("R848"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    CD4_Sofa.R848.sig.res <- ScreenGenes(CD4_hallmark.ssgsea.matrix, design = CD4_annotation_col, 
                                         contrast.col = "orig.ident", 
                                         contrast.level = c("Sofa","R848"),
                                         contrast.control = c("Sofa"),
                                         method = c("t.test","wilcox.test")[1],
                                         cutoff.p=0.05,
                                         cutoff.q=0.05,
                                         cutoff.logFC=1,
                                         pmethod = "BH",
                                         verbose = F)
    CD4_vehicle.R848.sig.res <- CD4_vehicle.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    CD4_vehicle.Sofa.sig.res <- CD4_vehicle.Sofa.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    CD4_vehicle.SR3.sig.res <- CD4_vehicle.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    CD4_Sofa.SR3.sig.res <- CD4_Sofa.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    CD4_R848.SR3.sig.res <- CD4_R848.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    CD4_Sofa.R848.sig.res <- CD4_Sofa.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN")))
    CD4_sig.hallmark <- unique(c(rownames(CD4_vehicle.R848.sig.res[CD4_vehicle.R848.sig.res$lab != "NO",]),rownames(CD4_vehicle.Sofa.sig.res[CD4_vehicle.Sofa.sig.res$lab != "NO",]),rownames(CD4_vehicle.SR3.sig.res[CD4_vehicle.SR3.sig.res$lab != "NO",])
                                 ,rownames(CD4_Sofa.SR3.sig.res[CD4_Sofa.SR3.sig.res$lab != "NO",]),rownames(CD4_R848.SR3.sig.res[CD4_R848.SR3.sig.res$lab != "NO",]),rownames(CD4_Sofa.R848.sig.res[CD4_Sofa.R848.sig.res$lab != "NO",])))
    plotHeat(CD4_hallmark.ssgsea.matrix[sig.hallmark,], 
             annotation_col=CD4_annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols <- colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col <- circlize::colorRamp2(c(-5,0,5), cols)
    sam_anno <- HeatmapAnnotation(orig.ident = CD4_annotation_col$orig.ident)
    p <- Heatmap(t(scale(t(CD4_hallmark.ssgsea.matrix[CD4_sig.hallmark,]))),show_column_names = F, column_split = CD4_annotation_col$orig.ident,col=col,top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),row_names_gp = gpar(fontsize = 6,fontface="bold")) 
    draw(p, annotation_legend_side = "left")
  }
  
  #CD8 CD4 KEGG 
  {
    CD8_signature.score <- calculate_sig_score(pdata = NULL,
                                               eset=as.matrix(CD8_eset),
                                               signature =kegg.signature,
                                               method = "zscore",
                                               mini_gene_count = 3,
                                               column_of_sample = "ID",
                                               print_gene_propotion = FALSE,
                                               adjust_eset = FALSE,
                                               print_filtered_signatures = FALSE)
    CD4_signature.score <- calculate_sig_score(pdata = NULL,
                                               eset=as.matrix(CD4_eset),
                                               signature =kegg.signature,
                                               method = "zscore",
                                               mini_gene_count = 3,
                                               column_of_sample = "ID",
                                               print_gene_propotion = FALSE,
                                               adjust_eset = FALSE,
                                               print_filtered_signatures = FALSE)
    CD8_hallmark.ssgsea.matrix <- CD8_signature.score %>% dplyr::select(.,-1) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    CD4_hallmark.ssgsea.matrix <- CD4_signature.score %>% dplyr::select(.,-1) %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>% t() %>%  as.matrix()
    CD8_annotation_col <- CD8_sce@meta.data 
    CD4_annotation_col <- CD4_sce@meta.data
    table(CD8_annotation_col$orig.ident)
    table(CD4_annotation_col$orig.ident)
    CD8_vehicle.R848.sig.res <- ScreenGenes(CD8_hallmark.ssgsea.matrix, design = CD8_annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","R848"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    CD8_vehicle.Sofa.sig.res <- ScreenGenes(CD8_hallmark.ssgsea.matrix, design = CD8_annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","Sofa"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    CD8_vehicle.SR3.sig.res <- ScreenGenes(CD8_hallmark.ssgsea.matrix, design = CD8_annotation_col, 
                                           contrast.col = "orig.ident", 
                                           contrast.level = c("Vehicle","SR3"),
                                           contrast.control = c("Vehicle"),
                                           method = c("t.test","wilcox.test")[1],
                                           cutoff.p=0.05,
                                           cutoff.q=0.05,
                                           cutoff.logFC=1,
                                           pmethod = "BH",
                                           verbose = F)
    CD8_Sofa.SR3.sig.res <- ScreenGenes(CD8_hallmark.ssgsea.matrix, design = CD8_annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Sofa","SR3"),
                                        contrast.control = c("Sofa"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    CD8_R848.SR3.sig.res <- ScreenGenes(CD8_hallmark.ssgsea.matrix, design = CD8_annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("R848","SR3"),
                                        contrast.control = c("R848"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    CD8_Sofa.R848.sig.res <- ScreenGenes(CD8_hallmark.ssgsea.matrix, design = CD8_annotation_col, 
                                         contrast.col = "orig.ident", 
                                         contrast.level = c("Sofa","R848"),
                                         contrast.control = c("Sofa"),
                                         method = c("t.test","wilcox.test")[1],
                                         cutoff.p=0.05,
                                         cutoff.q=0.05,
                                         cutoff.logFC=1,
                                         pmethod = "BH",
                                         verbose = F)
    CD8_vehicle.R848.sig.res <- CD8_vehicle.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    CD8_vehicle.Sofa.sig.res <- CD8_vehicle.Sofa.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    CD8_vehicle.SR3.sig.res <- CD8_vehicle.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    CD8_Sofa.SR3.sig.res <- CD8_Sofa.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    CD8_R848.SR3.sig.res <- CD8_R848.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    CD8_Sofa.R848.sig.res <- CD8_Sofa.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    CD8_sig.hallmark <- unique(c(rownames(CD8_vehicle.R848.sig.res[CD8_vehicle.R848.sig.res$lab != "NO",]),rownames(CD8_vehicle.Sofa.sig.res[CD8_vehicle.Sofa.sig.res$lab != "NO",]),rownames(CD8_vehicle.SR3.sig.res[CD8_vehicle.SR3.sig.res$lab != "NO",])
                                 ,rownames(CD8_Sofa.SR3.sig.res[CD8_Sofa.SR3.sig.res$lab != "NO",]),rownames(CD8_R848.SR3.sig.res[CD8_R848.SR3.sig.res$lab != "NO",]),rownames(CD8_Sofa.R848.sig.res[CD8_Sofa.R848.sig.res$lab != "NO",])))
    plotHeat(CD8_hallmark.ssgsea.matrix[CD8_sig.hallmark,], 
             annotation_col=annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols <- colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col <- circlize::colorRamp2(c(-5,0,5), cols)
    sam_anno <- HeatmapAnnotation(orig.ident = CD8_annotation_col$orig.ident)
    p <- Heatmap(use_raster=F,t(scale(t(CD8_hallmark.ssgsea.matrix[CD8_sig.hallmark,]))),show_column_names = F, column_split = CD8_annotation_col$orig.ident,col=col,top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),row_names_gp = gpar(fontsize = 6,fontface="bold")) 
    draw(p, annotation_legend_side = "left")
    
    CD4_vehicle.R848.sig.res <- ScreenGenes(CD4_hallmark.ssgsea.matrix, design = CD4_annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","R848"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    CD4_vehicle.Sofa.sig.res <- ScreenGenes(CD4_hallmark.ssgsea.matrix, design = CD4_annotation_col, 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","Sofa"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[1],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = F)
    CD4_vehicle.SR3.sig.res <- ScreenGenes(CD4_hallmark.ssgsea.matrix, design = CD4_annotation_col, 
                                           contrast.col = "orig.ident", 
                                           contrast.level = c("Vehicle","SR3"),
                                           contrast.control = c("Vehicle"),
                                           method = c("t.test","wilcox.test")[1],
                                           cutoff.p=0.05,
                                           cutoff.q=0.05,
                                           cutoff.logFC=1,
                                           pmethod = "BH",
                                           verbose = F)
    CD4_Sofa.SR3.sig.res <- ScreenGenes(CD4_hallmark.ssgsea.matrix, design = CD4_annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Sofa","SR3"),
                                        contrast.control = c("Sofa"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    CD4_R848.SR3.sig.res <- ScreenGenes(CD4_hallmark.ssgsea.matrix, design = CD4_annotation_col, 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("R848","SR3"),
                                        contrast.control = c("R848"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
    CD4_Sofa.R848.sig.res <- ScreenGenes(CD4_hallmark.ssgsea.matrix, design = CD4_annotation_col, 
                                         contrast.col = "orig.ident", 
                                         contrast.level = c("Sofa","R848"),
                                         contrast.control = c("Sofa"),
                                         method = c("t.test","wilcox.test")[1],
                                         cutoff.p=0.05,
                                         cutoff.q=0.05,
                                         cutoff.logFC=1,
                                         pmethod = "BH",
                                         verbose = F)
    CD4_vehicle.R848.sig.res <- CD4_vehicle.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    CD4_vehicle.Sofa.sig.res <- CD4_vehicle.Sofa.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    CD4_vehicle.SR3.sig.res <- CD4_vehicle.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    CD4_Sofa.SR3.sig.res <- CD4_Sofa.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    CD4_R848.SR3.sig.res <- CD4_R848.SR3.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    CD4_Sofa.R848.sig.res <- CD4_Sofa.R848.sig.res %>% dplyr::mutate(lab=ifelse( qvalue > 0.05, "NO", ifelse(mean.treat > mean.control,"UP","DOWN"))) %>% dplyr::arrange(.,qvalue) %>%  dplyr::top_n(50,desc(qvalue))
    CD4_sig.hallmark <- unique(c(rownames(CD4_vehicle.R848.sig.res[CD4_vehicle.R848.sig.res$lab != "NO",]),rownames(CD4_vehicle.Sofa.sig.res[CD4_vehicle.Sofa.sig.res$lab != "NO",]),rownames(CD4_vehicle.SR3.sig.res[CD4_vehicle.SR3.sig.res$lab != "NO",])
                                 ,rownames(CD4_Sofa.SR3.sig.res[CD4_Sofa.SR3.sig.res$lab != "NO",]),rownames(CD4_R848.SR3.sig.res[CD4_R848.SR3.sig.res$lab != "NO",]),rownames(CD4_Sofa.R848.sig.res[CD4_Sofa.R848.sig.res$lab != "NO",])))
    plotHeat(CD4_hallmark.ssgsea.matrix[sig.hallmark,], 
             annotation_col=CD4_annotation_col[,1:2], 
             cluster_cols = T,
             cluster_rows = F, 
             show_colnames = F, 
             show_rownames= T)
    library(ComplexHeatmap)
    col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols <- colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col <- circlize::colorRamp2(c(-5,0,5), cols)
    sam_anno <- HeatmapAnnotation(orig.ident = CD4_annotation_col$orig.ident)
    p <- Heatmap(use_raster=F,t(scale(t(CD4_hallmark.ssgsea.matrix[CD4_sig.hallmark,]))),show_column_names = F, column_split = CD4_annotation_col$orig.ident,col=col,top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),row_names_gp = gpar(fontsize = 6,fontface="bold")) 
    draw(p, annotation_legend_side = "left")
  }
  
  #T细胞中筛选各组上调的基因，作功能注释
  {
    #sce_merge_res1 <- readRDS("./Results/data/sce_merge_res1.rds")
    T_eset <- T_sce1@assays$RNA@data
    {
      T_annotation_col <- T_sce1@meta.data
      pos <- Fastgrep2(c("Vehicle","R848"),T_sce1@meta.data$orig.ident)  
      vehicle.R848.sig.genes <- ScreenGenes(T_eset[,unlist(pos)], design = T_annotation_col[unlist(pos),], 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","R848"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[2],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = T)
      table(vehicle.R848.sig.genes$lab)
      table(na.omit(vehicle.R848.sig.genes)$lab)
      vol <- ggplot(na.omit(vehicle.R848.sig.genes), aes(x = logFC, y = -log10(qvalue), color = lab))+  
        ggtitle(label = "R848 vs vehicle", subtitle = "Colored by fold-change direction") +
        geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
        scale_color_manual(name = "Directionality",
                           values = c(UP = "#008B00", DOWN = "#CD4F39", NO = "darkgray")) +
        theme_bw(base_size = 14) + # change overall theme
        theme(legend.position = "right") + # change the legend
        xlab(expression(log[2]("FoldChange"))) + # Change X-Axis label
        ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
        geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
        scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
      SR3_pos <- Fastgrep2(c("Vehicle","SR3"),T_sce1@meta.data$orig.ident)  
      vehicle.SR3.sig.genes <- ScreenGenes(T_eset[,unlist(SR3_pos)], design = T_annotation_col[unlist(SR3_pos),], 
                                           contrast.col = "orig.ident", 
                                           contrast.level = c("Vehicle","SR3"),
                                           contrast.control = c("Vehicle"),
                                           method = c("t.test","wilcox.test")[2],
                                           cutoff.p=0.05,
                                           cutoff.q=0.05,
                                           cutoff.logFC=1,
                                           pmethod = "BH",
                                           verbose = T)
      vol <- ggplot(na.omit(vehicle.SR3.sig.genes), aes(x = logFC, y = -log10(qvalue), color = lab))+  
        ggtitle(label = "SR3 vs vehicle", subtitle = "Colored by fold-change direction") +
        geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
        scale_color_manual(name = "Directionality",
                           values = c(UP = "#008B00", DOWN = "#CD4F39", NO = "darkgray")) +
        theme_bw(base_size = 14) + # change overall theme
        theme(legend.position = "right") + # change the legend
        xlab(expression(log[2]("FoldChange"))) + # Change X-Axis label
        ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
        geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
        scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
      Sofa_pos <- Fastgrep2(c("Vehicle","Sofa"),T_sce1@meta.data$orig.ident)  
      vehicle.Sofa.sig.genes <- ScreenGenes(T_eset[,unlist(Sofa_pos)], design = T_annotation_col[unlist(Sofa_pos),], 
                                            contrast.col = "orig.ident", 
                                            contrast.level = c("Vehicle","Sofa"),
                                            contrast.control = c("Vehicle"),
                                            method = c("t.test","wilcox.test")[2],
                                            cutoff.p=0.05,
                                            cutoff.q=0.05,
                                            cutoff.logFC=1,
                                            pmethod = "BH",
                                            verbose = T)
      vol <- ggplot(na.omit(vehicle.Sofa.sig.genes), aes(x = logFC, y = -log10(qvalue), color = lab))+  
        ggtitle(label = "Sofa vs vehicle", subtitle = "Colored by fold-change direction") +
        geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
        scale_color_manual(name = "Directionality",
                           values = c(UP = "#008B00", DOWN = "#CD4F39", NO = "darkgray")) +
        theme_bw(base_size = 14) + # change overall theme
        theme(legend.position = "right") + # change the legend
        xlab(expression(log[2]("FoldChange"))) + # Change X-Axis label
        ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
        geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
        scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
      pos <- Fastgrep2(c("R848","SR3"),T_sce1@meta.data$orig.ident)  
      R848.SR3.sig.genes <- ScreenGenes(T_eset[,unlist(pos)], design = T_annotation_col[unlist(pos),], 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("R848","SR3"),
                                        contrast.control = c("R848"),
                                        method = c("t.test","wilcox.test")[2],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = T)
      table(R848.SR3.sig.genes$lab)
      table(na.omit(R848.SR3.sig.genes)$lab)
      vol <- ggplot(na.omit(R848.SR3.sig.genes), aes(x = logFC, y = -log10(qvalue), color = lab))+  
        ggtitle(label = "SR3 vs R848", subtitle = "Colored by fold-change direction") +
        geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
        scale_color_manual(name = "Directionality",
                           values = c(UP = "#008B00", DOWN = "#CD4F39", NO = "darkgray")) +
        theme_bw(base_size = 14) + # change overall theme
        theme(legend.position = "right") + # change the legend
        xlab(expression(log[2]("FoldChange"))) + # Change X-Axis label
        ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
        geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
        scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
      pos <- Fastgrep2(c("Sofa","SR3"),T_sce1@meta.data$orig.ident)  
      Sofa.SR3.sig.genes <- ScreenGenes(T_eset[,unlist(pos)], design = T_annotation_col[unlist(pos),], 
                                        contrast.col = "orig.ident", 
                                        contrast.level = c("Sofa","SR3"),
                                        contrast.control = c("Sofa"),
                                        method = c("t.test","wilcox.test")[2],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = T)
      table(Sofa.SR3.sig.genes$lab)
      table(na.omit(Sofa.SR3.sig.genes)$lab)
      vol <- ggplot(na.omit(Sofa.SR3.sig.genes), aes(x = logFC, y = -log10(qvalue), color = lab))+  
        ggtitle(label = "SR3 vs Sofa", subtitle = "Colored by fold-change direction") +
        geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
        scale_color_manual(name = "Directionality",
                           values = c(UP = "#008B00", DOWN = "#CD4F39", NO = "darkgray")) +
        theme_bw(base_size = 14) + # change overall theme
        theme(legend.position = "right") + # change the legend
        xlab(expression(log[2]("FoldChange"))) + # Change X-Axis label
        ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
        geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
        scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
      pos <- Fastgrep2(c("R848","Sofa"),T_sce1@meta.data$orig.ident)  
      Sofa.R848.sig.genes <- ScreenGenes(T_eset[,unlist(pos)], design = T_annotation_col[unlist(pos),], 
                                         contrast.col = "orig.ident", 
                                         contrast.level = c("Sofa","R848"),
                                         contrast.control = c("Sofa"),
                                         method = c("t.test","wilcox.test")[2],
                                         cutoff.p=0.05,
                                         cutoff.q=0.05,
                                         cutoff.logFC=1,
                                         pmethod = "BH",
                                         verbose = T)
      table(Sofa.R848.sig.genes$lab)
      table(na.omit(Sofa.R848.sig.genes)$lab)
      vol <- ggplot(na.omit(Sofa.R848.sig.genes), aes(x = logFC, y = -log10(qvalue), color = lab))+  
        ggtitle(label = "R848 vs Sofa", subtitle = "Colored by fold-change direction") +
        geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
        scale_color_manual(name = "Directionality",
                           values = c(UP = "#008B00", DOWN = "#CD4F39", NO = "darkgray")) +
        theme_bw(base_size = 14) + # change overall theme
        theme(legend.position = "right") + # change the legend
        xlab(expression(log[2]("FoldChange"))) + # Change X-Axis label
        ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
        geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
        scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
      save(vehicle.R848.sig.genes,vehicle.SR3.sig.genes,vehicle.Sofa.sig.genes,R848.SR3.sig.genes,Sofa.SR3.sig.genes,Sofa.R848.sig.genes,file="./Results/mouse/diffExp.RData")
    }
    {
      T_sce1 <- readRDS(file="./Results/mouse/T_cell.rds")
      T_eset <- T_sce1@assays$RNA@data
      head(T_sce1@meta.data)
      Vehicle.markers <- FindMarkers(T_sce1,ident.1 = "Vehicle",ident.2 = NULL, group.by="orig.ident", only.pos = T)
      Sofa.markers <- FindMarkers(T_sce1,ident.1 = "Sofa",ident.2 = NULL, group.by="orig.ident", only.pos = T)
      R848.markers <- FindMarkers(T_sce1,ident.1 = "R848",ident.2 = NULL, group.by="orig.ident", only.pos = T)
      SR3.markers <- FindMarkers(T_sce1,ident.1 = "SR3",ident.2 = NULL, group.by="orig.ident", only.pos = T)
    
      Vehicle.markers1 <- Vehicle.markers %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=200) %>% rownames()
      Sofa.markers1 <- Sofa.markers %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=200) %>% rownames()
      R848.markers1 <- R848.markers %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=200) %>% rownames()
      SR3.markers1 <- SR3.markers %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=200) %>% rownames()
      
      #plot
      plot.Vehicle.markers1 <- Vehicle.markers  %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=20) %>% rownames()
      plot.Sofa.markers1 <- Sofa.markers  %>% filter(.,p_val_adj < 0.05)%>% dplyr::slice_max(avg_log2FC,n=20) %>% rownames()
      plot.R848.markers1 <- R848.markers  %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=20) %>% rownames()
      plot.SR3.markers1 <- SR3.markers  %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=20) %>% rownames()

      plot.markers <- unique(c(plot.Vehicle.markers1,plot.Sofa.markers1,plot.R848.markers1,plot.SR3.markers1))
      dh <- DoHeatmap(T_sce1, features = plot.markers,group.by = "orig.ident") + NoLegend()
      print(dh)
    
      }
    
    
    ####功能注释
    {
      mmu.annot <- getAnnotation(gtf.path = "/Users/biofly/project/shijian/gencode.vM28.annotation.gtf.gz")
      mmu.annot <- getAnnotation2(mmu.annot,organism = "mouse")
      
      #up_genes <-  as.character(na.omit(as.character(vehicle.R848.sig.genes$X1[vehicle.R848.sig.genes$lab=="UP"])))
      up_genes <- as.character(na.omit(convert(Vehicle.markers1, fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
      geneList <- T_eset[,1]
      names(geneList) <- convert(names(geneList), fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)
      geneList <- geneList[!is.na(names(geneList))]
      Vehicle.markers.kegg <- FastKEGG(genes = up_genes,
                                              geneList = geneList,
                                              default.universe = F,
                                              organism = 'mmu',
                                              db = mmu.annot,
                                              pvalueCutoff = list(enrichKEGG = 0.05,
                                                                  enrichMKEGG = 0.05),
                                              qvalueCutoff = 0.05,
                                              cnet.showCategory = 10,
                                              verbose = T,
                                              save.path = "KEGG",
                                              names = "Vehicle.markers")
      Vehicle.markers.go <- FastGO(genes = up_genes,
                                          geneList = geneList,
                                          organism="mouse",
                                          default.universe = F,
                                          classlevel = 2:2,
                                          OrgDb  = "org.Mm.eg.db",
                                          keyType = NULL,
                                          pAdjustMethod = "BH",
                                          pvalueCutoff = 0.05,
                                          qvalueCutoff  = 0.05,
                                          cnet.showCategory = 5,
                                          verbose = TRUE,
                                          save.path = "GO",
                                          names = "Vehicle.markers")
      Vehicle.markers.ms <- FastMS(genes=up_genes,
                                          geneList=geneList,
                                          default.universe = F,
                                          pAdjustMethod = "BH",
                                          db = mmu.annot,
                                          organism = "mouse",
                                          pvalueCutoff = list(enrichMS=0.05,
                                                              gseMS = 0.05),
                                          qvalueCutoff  = 0.05,
                                          cnet.showCategory = 5,
                                          verbose = TRUE,
                                          save.path = "Molecular Signatures",
                                          names = "Vehicle.markers")
      object <- vehicle.R848.sig.genes.kegg$enrichKEGG
      FastEnrichPlot(object,
                     pvalueCutoff = NULL,
                     qvalueCutoff = NULL,
                     id.col = "Description",
                     select = NULL,
                     x.title = "Kegg pathway",
                     y.title = "-log10(FDR)",
                     size = 15,
                     short.cutoff = 46,
                     topshow = 10)
      
      #up_genes <-  as.character(na.omit(as.character(vehicle.SR3.sig.genes$X1[vehicle.SR3.sig.genes$lab=="UP"])))
      up_genes <- as.character(na.omit(convert(Sofa.markers1, fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
      geneList <- T_eset[,1]
      names(geneList) <- convert(names(geneList), fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)
      geneList <- geneList[!is.na(names(geneList))]
      Sofa.markers.kegg <- FastKEGG(genes = up_genes,
                                              geneList = geneList,
                                              default.universe = F,
                                              organism = 'mmu',
                                              db = mmu.annot,
                                              pvalueCutoff = list(enrichKEGG = 0.05,
                                                                  enrichMKEGG = 0.05),
                                              qvalueCutoff = 0.05,
                                              cnet.showCategory = 10,
                                              verbose = T,
                                              save.path = "KEGG",
                                              names = "Sofa.markers")
      Sofa.markers.go <- FastGO(genes = up_genes,
                                          geneList = geneList,
                                          organism="mouse",
                                          default.universe = F,
                                          classlevel = 2:2,
                                          OrgDb  = "org.Mm.eg.db",
                                          keyType = NULL,
                                          pAdjustMethod = "BH",
                                          pvalueCutoff = 0.05,
                                          qvalueCutoff  = 0.05,
                                          cnet.showCategory = 5,
                                          verbose = TRUE,
                                          save.path = "GO",
                                          names = "Sofa.markers")
      Sofa.markers.ms <- FastMS(genes=up_genes,
                                          geneList=geneList,
                                          default.universe = F,
                                          pAdjustMethod = "BH",
                                          db = mmu.annot,
                                          organism = "mouse",
                                          pvalueCutoff = list(enrichMS=0.05,
                                                              gseMS = 0.05),
                                          qvalueCutoff  = 0.05,
                                          cnet.showCategory = 5,
                                          verbose = TRUE,
                                          save.path = "Molecular Signatures",
                                          names = "Sofa.markers")
      object <- vehicle.SR3.sig.genes.kegg$enrichKEGG
      FastEnrichPlot(object,
                     pvalueCutoff = NULL,
                     qvalueCutoff = NULL,
                     id.col = "Description",
                     select = NULL,
                     x.title = "Kegg pathway",
                     y.title = "-log10(FDR)",
                     size = 15,
                     short.cutoff = 46,
                     topshow = 4)
      
      #up_genes <-  as.character(na.omit(as.character(vehicle.Sofa.sig.genes$X1[vehicle.Sofa.sig.genes$lab=="UP"])))
      up_genes <- as.character(na.omit(convert(R848.markers1, fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
      geneList <- T_eset[,1]
      names(geneList) <- convert(names(geneList), fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)
      geneList <- geneList[!is.na(names(geneList))]
      R848.markers.kegg <- FastKEGG(genes = up_genes,
                                             geneList = geneList,
                                             default.universe = F,
                                             organism = 'mmu',
                                             db = mmu.annot,
                                             pvalueCutoff = list(enrichKEGG = 0.05,
                                                                 enrichMKEGG = 0.05),
                                             qvalueCutoff = 0.05,
                                             cnet.showCategory = 10,
                                             verbose = T,
                                             save.path = "KEGG",
                                             names = "R848.markers")
      R848.markers.go <- FastGO(genes = up_genes,
                                         geneList = geneList,
                                         organism="mouse",
                                         default.universe = F,
                                         classlevel = 2:2,
                                         OrgDb  = "org.Mm.eg.db",
                                         keyType = NULL,
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff  = 0.05,
                                         cnet.showCategory = 5,
                                         verbose = TRUE,
                                         save.path = "GO",
                                         names = "R848.markers")
      R848.markers.ms <- FastMS(genes=up_genes,
                                         geneList=geneList,
                                         default.universe = F,
                                         pAdjustMethod = "BH",
                                         db = mmu.annot,
                                         organism = "mouse",
                                         pvalueCutoff = list(enrichMS=0.05,
                                                             gseMS = 0.05),
                                         qvalueCutoff  = 0.05,
                                         cnet.showCategory = 5,
                                         verbose = TRUE,
                                         save.path = "Molecular Signatures",
                                         names = "R848.markers")
      object <- vehicle.Sofa.sig.genes.kegg$enrichKEGG
      FastEnrichPlot(object,
                     pvalueCutoff = NULL,
                     qvalueCutoff = NULL,
                     id.col = "Description",
                     select = NULL,
                     x.title = "Kegg pathway",
                     y.title = "-log10(FDR)",
                     size = 15,
                     short.cutoff = 46,
                     topshow = 10)
      
      #up_genes <-  as.character(na.omit(as.character(R848.SR3.sig.genes$X1[R848.SR3.sig.genes$lab=="UP"])))
      up_genes <- as.character(na.omit(convert(SR3.markers1, fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
      geneList <- T_eset[,1]
      names(geneList) <- convert(names(geneList), fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)
      geneList <- geneList[!is.na(names(geneList))]
      SR3.markers.kegg <- FastKEGG(genes = up_genes,
                                           geneList = geneList,
                                           default.universe = F,
                                           organism = 'mmu',
                                           db = mmu.annot,
                                           pvalueCutoff = list(enrichKEGG = 0.05,
                                                               enrichMKEGG = 0.05),
                                           qvalueCutoff = 0.05,
                                           cnet.showCategory = 10,
                                           verbose = T,
                                           save.path = "KEGG",
                                           names = "SR3.markers")
      SR3.markers.go <- FastGO(genes = up_genes,
                                       geneList = geneList,
                                       organism="mouse",
                                       default.universe = F,
                                       classlevel = 2:2,
                                       OrgDb  = "org.Mm.eg.db",
                                       keyType = NULL,
                                       pAdjustMethod = "BH",
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff  = 0.05,
                                       cnet.showCategory = 5,
                                       verbose = TRUE,
                                       save.path = "GO",
                                       names = "SR3.markers")
      SR3.markers.ms <- FastMS(genes=up_genes,
                                       geneList=geneList,
                                       default.universe = F,
                                       pAdjustMethod = "BH",
                                       db = mmu.annot,
                                       organism = "mouse",
                                       pvalueCutoff = list(enrichMS=0.05,
                                                           gseMS = 0.05),
                                       qvalueCutoff  = 0.05,
                                       cnet.showCategory = 5,
                                       verbose = TRUE,
                                       save.path = "Molecular Signatures",
                                       names = "SR3.markers")
      object <- R848.SR3.sig.genes.kegg$enrichKEGG
      FastEnrichPlot(object,
                     pvalueCutoff = NULL,
                     qvalueCutoff = NULL,
                     id.col = "Description",
                     select = NULL,
                     x.title = "Kegg pathway",
                     y.title = "-log10(FDR)",
                     size = 15,
                     short.cutoff = 46,
                     topshow = 10)
      save(Vehicle.markers.kegg,Vehicle.markers.go,Vehicle.markers.ms,
           Sofa.markers.kegg,Sofa.markers.go,Sofa.markers.ms,
           R848.markers.kegg,R848.markers.go,R848.markers.ms,
           SR3.markers.kegg,SR3.markers.go,SR3.markers.ms
           ,file = "./Results/mouse/functionEnrichment.RData")
      load(file="./Results/mouse/functionEnrichment.RData")
      #vehicle.Sofa.sig.genes.kegg$enrichKEGG@result %>% head()
      kegg1 <- Vehicle.markers.kegg$enrichKEGG@result %>% 
         mutate(sample=rep("Vehicle",length(ID))) %>% 
         filter(.,qvalue < 0.05) 
      kegg2 <- Sofa.markers.kegg$enrichKEGG@result %>% 
        mutate(sample=rep("Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05)
      kegg3 <- R848.markers.kegg$enrichKEGG@result %>% 
        mutate(sample=rep("R848",length(ID))) %>% 
        filter(.,qvalue < 0.05)
      kegg4 <- SR3.markers.kegg$enrichKEGG@result %>% 
        mutate(sample=rep("SR3",length(ID))) %>% 
        filter(.,qvalue < 0.05)
      kegg <- kegg1 %>% bind_rows(.,kegg2) %>% 
                        bind_rows(.,kegg3) %>% 
                        bind_rows(.,kegg4) 
       head(kegg)      
       ggplot(data=kegg,mapping=aes(x=factor(sample,levels = c("Vehicle","Sofa","R848","SR3")),y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
          geom_point() +
          theme_cowplot(font_size = 12) +
         scale_x_discrete(limits=c("Vehicle","Sofa","R848","SR3"),labels=c("Vehicle","Sofa","R848","SR3")) +
          labs(x="sample")
       
       
       go_bp1 <- Vehicle.markers.go$GOEnrichment$BP@result %>% 
         mutate(sample=rep("Vehicle",length(ID))) %>% 
         filter(.,qvalue < 0.05) %>% 
         slice_min(qvalue,n=10)
       go_bp2 <- Sofa.markers.go$GOEnrichment$BP@result %>% 
         mutate(sample=rep("Sofa",length(ID))) %>% 
         filter(.,qvalue < 0.05) %>% 
         slice_min(qvalue,n=10)
       go_bp3 <- R848.markers.go$GOEnrichment$BP@result %>% 
         mutate(sample=rep("R848",length(ID))) %>% 
         filter(.,qvalue < 0.05) %>% 
         slice_min(qvalue,n=10)
       go_bp4 <- SR3.markers.go$GOEnrichment$BP@result %>% 
         mutate(sample=rep("SR3",length(ID))) %>% 
         filter(.,qvalue < 0.05) %>% 
         slice_min(qvalue,n=10)
       
       go_bp <- go_bp1 %>% bind_rows(.,go_bp2) %>% 
         bind_rows(.,go_bp3) %>% 
         bind_rows(.,go_bp4) 
       head(go_bp)      
       ggplot(data=go_bp,mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
         geom_point() +
         theme_cowplot(font_size = 12) +
         labs(x="sample")
       
       hallmark1 <- Vehicle.markers.ms$Data$enrichMS$h@result %>% 
         mutate(sample=rep("Vehicle",length(ID))) %>% 
         filter(.,qvalue < 0.05) %>% 
         slice_min(qvalue,n=10)
       hallmark2 <- Sofa.markers.ms$Data$enrichMS$h@result %>% 
         mutate(sample=rep("Sofa",length(ID))) %>% 
         filter(.,qvalue < 0.05) %>% 
         slice_min(qvalue,n=10)
       hallmark3 <- R848.markers.ms$Data$enrichMS$h@result %>% 
         mutate(sample=rep("R848",length(ID))) %>% 
         filter(.,qvalue < 0.05) %>% 
         slice_min(qvalue,n=10)
       hallmark4 <- SR3.markers.ms$Data$enrichMS$h@result %>% 
         mutate(sample=rep("SR3",length(ID))) %>% 
         filter(.,qvalue < 0.05) %>% 
         slice_min(qvalue,n=10)
      
       hallmark <- hallmark1 %>% bind_rows(.,hallmark2) %>% 
         bind_rows(.,hallmark3) %>% 
         bind_rows(.,hallmark4)
       head(hallmark)      
       ggplot(data=hallmark,mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
         geom_point() +
         theme_cowplot(font_size = 12) +
         scale_x_discrete(limits=c("Vehicle","Sofa","R848","SR3"),labels=c("Vehicle","Sofa","R848","SR3")) +
         labs(x="sample")
    }
 
  }
  
  #巨噬细胞中筛选各组上调的基因，作功能注释
  {
    Marcrophage_sce <- readRDS(file="./Results/mouse/Marcrophage.rds")
    Vehicle.markers <- FindMarkers(Marcrophage_sce,ident.1 = "Vehicle",ident.2 = NULL, group.by="orig.ident", only.pos = T)
    Sofa.markers <- FindMarkers(Marcrophage_sce,ident.1 = "Sofa",ident.2 = NULL, group.by="orig.ident", only.pos = T)
    R848.markers <- FindMarkers(Marcrophage_sce,ident.1 = "R848",ident.2 = NULL, group.by="orig.ident", only.pos = T)
    SR3.markers <- FindMarkers(Marcrophage_sce,ident.1 = "SR3",ident.2 = NULL, group.by="orig.ident", only.pos = T)
    
    Vehicle.markers1 <- Vehicle.markers %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=200) %>% rownames()
    Sofa.markers1 <- Sofa.markers %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=200) %>% rownames()
    R848.markers1 <- R848.markers %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=200) %>% rownames()
    SR3.markers1 <- SR3.markers %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=200) %>% rownames()
    
    #plot
    plot.Vehicle.markers1 <- Vehicle.markers  %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=5) %>% rownames()
    plot.Sofa.markers1 <- Sofa.markers  %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=5) %>% rownames()
    plot.R848.markers1 <- R848.markers  %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=5) %>% rownames()
    plot.SR3.markers1 <- SR3.markers  %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=5) %>% rownames()
    
    plot.markers <- unique(c(plot.Vehicle.markers1,plot.Sofa.markers1,plot.R848.markers1,plot.SR3.markers1))
    dh <- DoHeatmap(Marcrophage_sce, features = plot.markers,group.by = "orig.ident") + NoLegend()
    print(dh)
    
    ####功能注释
    {
      mmu.annot <- getAnnotation(gtf.path = "/Users/biofly/project/shijian/gencode.vM28.annotation.gtf.gz")
      mmu.annot <- getAnnotation2(mmu.annot,organism = "mouse")
      #up_genes <-  as.character(na.omit(as.character(vehicle.R848.sig.genes$X1[vehicle.R848.sig.genes$lab=="UP"])))
      up_genes <- as.character(na.omit(convert(Vehicle.markers1, fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
      geneList <- Marcrophage_sce@assays$RNA@data[,1]
      names(geneList) <- convert(names(geneList), fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)
      geneList <- geneList[!is.na(names(geneList))]
      Vehicle.markers.kegg <- FastKEGG(genes = up_genes,
                                              geneList = geneList,
                                              default.universe = F,
                                              organism = 'mmu',
                                              db = mmu.annot,
                                              pvalueCutoff = list(enrichKEGG = 0.05,
                                                                  enrichMKEGG = 0.05),
                                              qvalueCutoff = 0.05,
                                              cnet.showCategory = 10,
                                              verbose = T,
                                              save.path = "KEGG",
                                              names = "marcrophage.vehicle.markers")
      Vehicle.markers.go <- FastGO(genes = up_genes,
                                          geneList = geneList,
                                          organism="mouse",
                                          default.universe = F,
                                          classlevel = 2:2,
                                          OrgDb  = "org.Mm.eg.db",
                                          keyType = NULL,
                                          pAdjustMethod = "BH",
                                          pvalueCutoff = 0.05,
                                          qvalueCutoff  = 0.05,
                                          cnet.showCategory = 5,
                                          verbose = TRUE,
                                          save.path = "GO",
                                          names = "marcrophage.Vehicle.markers")
      Vehicle.markers.ms <- FastMS(genes=up_genes,
                                          geneList=geneList,
                                          default.universe = F,
                                          pAdjustMethod = "BH",
                                          db = mmu.annot,
                                          organism = "mouse",
                                          pvalueCutoff = list(enrichMS=0.05,
                                                              gseMS = 0.05),
                                          qvalueCutoff  = 0.05,
                                          cnet.showCategory = 5,
                                          verbose = TRUE,
                                          save.path = "Molecular Signatures",
                                          names = "marcrophage.Vehicle.markers")
      object <- vehicle.R848.sig.genes.kegg$enrichKEGG
      FastEnrichPlot(object,
                     pvalueCutoff = NULL,
                     qvalueCutoff = NULL,
                     id.col = "Description",
                     select = NULL,
                     x.title = "Kegg pathway",
                     y.title = "-log10(FDR)",
                     size = 15,
                     short.cutoff = 46,
                     topshow = 10)
      
      #up_genes <-  as.character(na.omit(as.character(vehicle.SR3.sig.genes$X1[vehicle.SR3.sig.genes$lab=="UP"])))
      up_genes <- as.character(na.omit(convert(Sofa.markers1, fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
      geneList <- Marcrophage_sce@assays$RNA@data[,1]
      names(geneList) <- convert(names(geneList), fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)
      geneList <- geneList[!is.na(names(geneList))]
      Sofa.markers.kegg <- FastKEGG(genes = up_genes,
                                              geneList = geneList,
                                              default.universe = F,
                                              organism = 'mmu',
                                              db = mmu.annot,
                                              pvalueCutoff = list(enrichKEGG = 0.05,
                                                                  enrichMKEGG = 0.05),
                                              qvalueCutoff = 0.05,
                                              cnet.showCategory = 10,
                                              verbose = T,
                                              save.path = "KEGG",
                                              names = "marcrophage.Sofa.markers")
      Sofa.markers.go <- FastGO(genes = up_genes,
                                          geneList = geneList,
                                          organism="mouse",
                                          default.universe = F,
                                          classlevel = 2:2,
                                          OrgDb  = "org.Mm.eg.db",
                                          keyType = NULL,
                                          pAdjustMethod = "BH",
                                          pvalueCutoff = 0.05,
                                          qvalueCutoff  = 0.05,
                                          cnet.showCategory = 5,
                                          verbose = TRUE,
                                          save.path = "GO",
                                          names = "marcrophage.Sofa.markers")
      Sofa.markers.ms <- FastMS(genes=up_genes,
                                          geneList=geneList,
                                          default.universe = F,
                                          pAdjustMethod = "BH",
                                          db = mmu.annot,
                                          organism = "mouse",
                                          pvalueCutoff = list(enrichMS=0.05,
                                                              gseMS = 0.05),
                                          qvalueCutoff  = 0.05,
                                          cnet.showCategory = 5,
                                          verbose = TRUE,
                                          save.path = "Molecular Signatures",
                                          names = "marcrophage.Sofa.markers")
      object <- vehicle.SR3.sig.genes.kegg$enrichKEGG
      FastEnrichPlot(object,
                     pvalueCutoff = NULL,
                     qvalueCutoff = NULL,
                     id.col = "Description",
                     select = NULL,
                     x.title = "Kegg pathway",
                     y.title = "-log10(FDR)",
                     size = 15,
                     short.cutoff = 46,
                     topshow = 4)
      
      #up_genes <-  as.character(na.omit(as.character(vehicle.Sofa.sig.genes$X1[vehicle.Sofa.sig.genes$lab=="UP"])))
      up_genes <- as.character(na.omit(convert(R848.markers1, fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
      names(geneList) <- convert(names(geneList), fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)
      geneList <- geneList[!is.na(names(geneList))]
      R848.markers.kegg <- FastKEGG(genes = up_genes,
                                             geneList = geneList,
                                             default.universe = F,
                                             organism = 'mmu',
                                             db = mmu.annot,
                                             pvalueCutoff = list(enrichKEGG = 0.05,
                                                                 enrichMKEGG = 0.05),
                                             qvalueCutoff = 0.05,
                                             cnet.showCategory = 10,
                                             verbose = T,
                                             save.path = "KEGG",
                                             names = "marcrophage.R848.markers")
      R848.markers.go <- FastGO(genes = up_genes,
                                         geneList = geneList,
                                         organism="mouse",
                                         default.universe = F,
                                         classlevel = 2:2,
                                         OrgDb  = "org.Mm.eg.db",
                                         keyType = NULL,
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff  = 0.05,
                                         cnet.showCategory = 5,
                                         verbose = TRUE,
                                         save.path = "GO",
                                         names = "marcrophage.R848.markers")
      R848.markers.ms <- FastMS(genes=up_genes,
                                         geneList=geneList,
                                         default.universe = F,
                                         pAdjustMethod = "BH",
                                         db = mmu.annot,
                                         organism = "mouse",
                                         pvalueCutoff = list(enrichMS=0.05,
                                                             gseMS = 0.05),
                                         qvalueCutoff  = 0.05,
                                         cnet.showCategory = 5,
                                         verbose = TRUE,
                                         save.path = "Molecular Signatures",
                                         names = "marcrophage.R848.markers")
      object <- vehicle.Sofa.sig.genes.kegg$enrichKEGG
      FastEnrichPlot(object,
                     pvalueCutoff = NULL,
                     qvalueCutoff = NULL,
                     id.col = "Description",
                     select = NULL,
                     x.title = "Kegg pathway",
                     y.title = "-log10(FDR)",
                     size = 15,
                     short.cutoff = 46,
                     topshow = 10)
      
      #up_genes <-  as.character(na.omit(as.character(R848.SR3.sig.genes$X1[R848.SR3.sig.genes$lab=="UP"])))
      up_genes <- as.character(na.omit(convert(SR3.markers1, fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
      names(geneList) <- convert(names(geneList), fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)
      geneList <- geneList[!is.na(names(geneList))]
      SR3.markers.kegg <- FastKEGG(genes = up_genes,
                                           geneList = geneList,
                                           default.universe = F,
                                           organism = 'mmu',
                                           db = mmu.annot,
                                           pvalueCutoff = list(enrichKEGG = 0.05,
                                                               enrichMKEGG = 0.05),
                                           qvalueCutoff = 0.05,
                                           cnet.showCategory = 10,
                                           verbose = T,
                                           save.path = "KEGG",
                                           names = "marcrophage.SR3.markers")
      SR3.markers.go <- FastGO(genes = up_genes,
                                       geneList = geneList,
                                       organism="mouse",
                                       default.universe = F,
                                       classlevel = 2:2,
                                       OrgDb  = "org.Mm.eg.db",
                                       keyType = NULL,
                                       pAdjustMethod = "BH",
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff  = 0.05,
                                       cnet.showCategory = 5,
                                       verbose = TRUE,
                                       save.path = "GO",
                                       names = "marcrophage.SR3.markers")
      SR3.markers.ms <- FastMS(genes=up_genes,
                                       geneList=geneList,
                                       default.universe = F,
                                       pAdjustMethod = "BH",
                                       db = mmu.annot,
                                       organism = "mouse",
                                       pvalueCutoff = list(enrichMS=0.05,
                                                           gseMS = 0.05),
                                       qvalueCutoff  = 0.05,
                                       cnet.showCategory = 5,
                                       verbose = TRUE,
                                       save.path = "Molecular Signatures",
                                       names = "marcrophage.SR3.markers")
      object <- R848.SR3.sig.genes.kegg$enrichKEGG
      FastEnrichPlot(object,
                     pvalueCutoff = NULL,
                     qvalueCutoff = NULL,
                     id.col = "Description",
                     select = NULL,
                     x.title = "Kegg pathway",
                     y.title = "-log10(FDR)",
                     size = 15,
                     short.cutoff = 46,
                     topshow = 10)

     save(
           Vehicle.markers.kegg,Vehicle.markers.go,Vehicle.markers.ms,
           Sofa.markers.kegg,Sofa.markers.go,Sofa.markers.ms,
           R848.markers.kegg,R848.markers.go,R848.markers.ms,
           SR3.markers.kegg,SR3.markers.go,SR3.markers.ms
           ,file = "./Results/mouse/marcrophage.functionEnrichment.RData")
      #vehicle.Sofa.sig.genes.kegg$enrichKEGG@result %>% head()
      kegg1 <- Vehicle.markers.kegg$enrichKEGG@result %>% 
        mutate(sample=rep("Vehicle",length(ID))) %>% 
        filter(.,qvalue < 0.05) 
      kegg2 <- Sofa.markers.kegg$enrichKEGG@result %>% 
        mutate(sample=rep("Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05)
      kegg3 <- R848.markers.kegg$enrichKEGG@result %>% 
        mutate(sample=rep("R848",length(ID))) %>% 
        filter(.,qvalue < 0.05)
      kegg4 <- SR3.markers.kegg$enrichKEGG@result %>% 
        mutate(sample=rep("SR3",length(ID))) %>% 
        filter(.,qvalue < 0.05)
     
      kegg <- kegg1 %>% bind_rows(.,kegg2) %>% 
        bind_rows(.,kegg3) %>% 
        bind_rows(.,kegg4) 
      head(kegg)      
      ggplot(data=kegg,mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
        geom_point() +
        theme_cowplot(font_size = 12) +
        labs(x="sample")
      
      
      go_bp1 <- Vehicle.markers.go$GOEnrichment$BP@result %>% 
        mutate(sample=rep("Vehicle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_bp2 <- Sofa.markers.go$GOEnrichment$BP@result %>% 
        mutate(sample=rep("Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_bp3 <- R848.markers.go$GOEnrichment$BP@result %>% 
        mutate(sample=rep("R848",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_bp4 <- SR3.markers.go$GOEnrichment$BP@result %>% 
        mutate(sample=rep("SR3",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
     
      go_bp <- go_bp1 %>% bind_rows(.,go_bp2) %>% 
        bind_rows(.,go_bp3) %>% 
        bind_rows(.,go_bp4) 
      head(go_bp)      
      ggplot(data=go_bp,mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
        geom_point() +
        theme_cowplot(font_size = 12) +
        labs(x="sample")
      
      go_mf1 <- vehicle.Sofa.sig.genes.go$GOEnrichment$MF@result %>% 
        mutate(sample=rep("Sofa vs Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_mf2 <- vehicle.R848.sig.genes.go$GOEnrichment$MF@result %>% 
        mutate(sample=rep("R848 vs Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_mf3 <- vehicle.SR3.sig.genes.go$GOEnrichment$MF@result %>% 
        mutate(sample=rep("SR3 vs Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_mf4 <- SR3.markers.go$GOEnrichment$MF@result %>% 
        mutate(sample=rep("R848 vs Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
     
      go_mf <- go_mf1 %>% bind_rows(.,go_mf2) %>% 
        bind_rows(.,go_mf3) %>% 
        bind_rows(.,go_mf4) %>% 
        bind_rows(.,go_mf5) %>% 
        bind_rows(.,go_mf6) 
      head(go_mf)      
      ggplot(data=go_mf,mapping=aes(x=sample,y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
        geom_point() +
        theme_cowplot(font_size = 12)
      
      
      go_cc1 <- vehicle.Sofa.sig.genes.go$GOEnrichment$CC@result %>% 
        mutate(sample=rep("Sofa vs Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_cc2 <- vehicle.R848.sig.genes.go$GOEnrichment$CC@result %>% 
        mutate(sample=rep("R848 vs Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_cc3 <- vehicle.SR3.sig.genes.go$GOEnrichment$CC@result %>% 
        mutate(sample=rep("SR3 vs Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_cc4 <- Sofa.R848.sig.genes.go$GOEnrichment$CC@result %>% 
        mutate(sample=rep("R848 vs Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_cc5 <- Sofa.SR3.sig.genes.go$GOEnrichment$CC@result %>% 
        mutate(sample=rep("SR3 vs Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_cc6 <-  R848.SR3.sig.genes.go$GOEnrichment$CC@result %>% 
        mutate(sample=rep("SR3 vs R848",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_cc <- go_cc1 %>% bind_rows(.,go_cc2) %>% 
        bind_rows(.,go_cc3) %>% 
        bind_rows(.,go_cc4) %>% 
        bind_rows(.,go_cc5) %>% 
        bind_rows(.,go_cc6) 
      head(go_cc)  
      
      ggplot(data=go_cc,mapping=aes(x=sample,y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
        geom_point() +
        theme_cowplot(font_size = 12)
      
      
      hallmark1 <- Vehicle.markers.ms$Data$enrichMS$h@result %>% 
        mutate(sample=rep("Vehicle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      hallmark2 <- Sofa.markers.ms$Data$enrichMS$h@result %>% 
        mutate(sample=rep("Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      hallmark3 <- R848.markers.ms$Data$enrichMS$h@result %>% 
        mutate(sample=rep("R848",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      hallmark4 <- SR3.markers.ms$Data$enrichMS$h@result %>% 
        mutate(sample=rep("SR3",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
     
      hallmark <- hallmark1 %>% bind_rows(.,hallmark2) %>% 
        bind_rows(.,hallmark3) %>% 
        bind_rows(.,hallmark4) 
      head(hallmark)      
      ggplot(data=hallmark,mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
        geom_point() +
        theme_cowplot(font_size = 12) +
        labs(x="sample")
    }
  }
  
  #中性粒细胞中筛选各组上调的基因，作功能注释
  {
    sce_merge_res1 <- readRDS("./Results/data/sce_merge_res1.rds")
    #细胞打标签
    sce_merge_res1$sce@meta.data$celltype <- as.character(plyr::mapvalues(sce_merge_res1$sce@meta.data$seurat_clusters,from = c(0:4,6:16),
                                                                          to=c("Epithelial_cell","CD8+T_cell","Neutrophil_cell","Neutrophil_cell","Marcrophage",
                                                                               "CD8+T_cell","NK","CD4+T_cell","Stromal_cell","Dendritic_cell","Epithelial_cell",
                                                                               "CD8+T_cell","B_cell","Dendritic_cell","Plasma","Mast_cell")))
    sce_merge_res1$sce@meta.data$celltype1 <- as.character(plyr::mapvalues(sce_merge_res1$sce@meta.data$seurat_clusters,from = c(0:4,6:16),
                                                                           to=c("Epithelial_cell","T_cell","Neutrophil_cell","Neutrophil_cell","Marcrophage",
                                                                                "T_cell","NK","T_cell","Stromal_cell","Dendritic_cell","Epithelial_cell",
                                                                                "T_cell","B_cell","Dendritic_cell","Plasma","Mast_cell")))
    Neutrophil.sce<- subset(sce_merge_res1$sce,celltype == "Neutrophil_cell")
    Vehicle.markers <- FindMarkers(Neutrophil.sce,ident.1 = "Vehicle",ident.2 = NULL, group.by="orig.ident", only.pos = T)
    Sofa.markers <- FindMarkers(Neutrophil.sce,ident.1 = "Sofa",ident.2 = NULL, group.by="orig.ident", only.pos = T)
    R848.markers <- FindMarkers(Neutrophil.sce,ident.1 = "R848",ident.2 = NULL, group.by="orig.ident", only.pos = T)
    SR3.markers <- FindMarkers(Neutrophil.sce,ident.1 = "SR3",ident.2 = NULL, group.by="orig.ident", only.pos = T)
    
    
    Vehicle.markers1 <- Vehicle.markers %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=200) %>% rownames()
    Sofa.markers1 <- Sofa.markers %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=200) %>% rownames()
    R848.markers1 <- R848.markers %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=200) %>% rownames()
    SR3.markers1 <- SR3.markers %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=200) %>% rownames()

    
    #plot
    plot.Vehicle.markers1 <- Vehicle.markers  %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=10) %>% rownames()
    plot.Sofa.markers1 <- Sofa.markers  %>% filter(.,p_val_adj < 0.05)%>% dplyr::slice_max(avg_log2FC,n=10) %>% rownames()
    plot.R848.markers1 <- R848.markers  %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=10) %>% rownames()
    plot.SR3.markers1 <- SR3.markers  %>% filter(.,p_val_adj < 0.05) %>% dplyr::slice_max(avg_log2FC,n=10) %>% rownames()
    
    plot.markers <- unique(c(plot.Vehicle.markers1,plot.Sofa.markers1,plot.R848.markers1,plot.SR3.markers1))
    dh <- DoHeatmap(Neutrophil.sce, features = plot.markers,group.by = "orig.ident") + NoLegend()
    print(dh)
    
    ####功能注释
    {
      mmu.annot <- getAnnotation(gtf.path = "/Users/biofly/project/shijian/gencode.vM28.annotation.gtf.gz")
      mmu.annot <- getAnnotation2(mmu.annot,organism = "mouse")
      up_genes <- as.character(na.omit(convert(Vehicle.markers1, fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
      geneList <- Neutrophil.sce@assays$RNA@data[,1]
      names(geneList) <- convert(names(geneList), fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)
      geneList <- geneList[!is.na(names(geneList))]
      Vehicle.markers.kegg <- FastKEGG(genes = up_genes,
                                              geneList = geneList,
                                              default.universe = F,
                                              organism = 'mmu',
                                              db = mmu.annot,
                                              pvalueCutoff = list(enrichKEGG = 0.05,
                                                                  enrichMKEGG = 0.05),
                                              qvalueCutoff = 0.05,
                                              cnet.showCategory = 10,
                                              verbose = T,
                                              save.path = "KEGG",
                                              names = "neutrophil.Vehicle.markers")
      Vehicle.markers.go <- FastGO(genes = up_genes,
                                          geneList = geneList,
                                          organism="mouse",
                                          default.universe = F,
                                          classlevel = 2:2,
                                          OrgDb  = "org.Mm.eg.db",
                                          keyType = NULL,
                                          pAdjustMethod = "BH",
                                          pvalueCutoff = 0.05,
                                          qvalueCutoff  = 0.05,
                                          cnet.showCategory = 5,
                                          verbose = TRUE,
                                          save.path = "GO",
                                          names = "neutrophil.Vehicle.markers")
      Vehicle.markers.ms <- FastMS(genes=up_genes,
                                          geneList=geneList,
                                          default.universe = F,
                                          pAdjustMethod = "BH",
                                          db = mmu.annot,
                                          organism = "mouse",
                                          pvalueCutoff = list(enrichMS=0.05,
                                                              gseMS = 0.05),
                                          qvalueCutoff  = 0.05,
                                          cnet.showCategory = 5,
                                          verbose = TRUE,
                                          save.path = "Molecular Signatures",
                                          names = "neutrophil.Vehicle.markers")
      object <- vehicle.R848.sig.genes.kegg$enrichKEGG
      FastEnrichPlot(object,
                     pvalueCutoff = NULL,
                     qvalueCutoff = NULL,
                     id.col = "Description",
                     select = NULL,
                     x.title = "Kegg pathway",
                     y.title = "-log10(FDR)",
                     size = 15,
                     short.cutoff = 46,
                     topshow = 10)

      up_genes <- as.character(na.omit(convert(Sofa.markers1, fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
      Sofa.markers.kegg <- FastKEGG(genes = up_genes,
                                              geneList = geneList,
                                              default.universe = F,
                                              organism = 'mmu',
                                              db = mmu.annot,
                                              pvalueCutoff = list(enrichKEGG = 0.05,
                                                                  enrichMKEGG = 0.05),
                                              qvalueCutoff = 0.05,
                                              cnet.showCategory = 10,
                                              verbose = T,
                                              save.path = "KEGG",
                                              names = "neutrophil.Sofa.markers")
      Sofa.markers.go <- FastGO(genes = up_genes,
                                          geneList = geneList,
                                          organism="mouse",
                                          default.universe = F,
                                          classlevel = 2:2,
                                          OrgDb  = "org.Mm.eg.db",
                                          keyType = NULL,
                                          pAdjustMethod = "BH",
                                          pvalueCutoff = 0.05,
                                          qvalueCutoff  = 0.05,
                                          cnet.showCategory = 5,
                                          verbose = TRUE,
                                          save.path = "GO",
                                          names = "neutrophil.Sofa.markers")
      Sofa.markers.ms <- FastMS(genes=up_genes,
                                          geneList=geneList,
                                          default.universe = F,
                                          pAdjustMethod = "BH",
                                          db = mmu.annot,
                                          organism = "mouse",
                                          pvalueCutoff = list(enrichMS=0.05,
                                                              gseMS = 0.05),
                                          qvalueCutoff  = 0.05,
                                          cnet.showCategory = 5,
                                          verbose = TRUE,
                                          save.path = "Molecular Signatures",
                                          names = "neutrophil.Sofa.markers")
      object <- vehicle.SR3.sig.genes.kegg$enrichKEGG
      FastEnrichPlot(object,
                     pvalueCutoff = NULL,
                     qvalueCutoff = NULL,
                     id.col = "Description",
                     select = NULL,
                     x.title = "Kegg pathway",
                     y.title = "-log10(FDR)",
                     size = 15,
                     short.cutoff = 46,
                     topshow = 4)
      
      #up_genes <-  as.character(na.omit(as.character(vehicle.Sofa.sig.genes$X1[vehicle.Sofa.sig.genes$lab=="UP"])))
      up_genes <- as.character(na.omit(convert(R848.markers1, fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
      R848.markers.kegg <- FastKEGG(genes = up_genes,
                                             geneList = geneList,
                                             default.universe = F,
                                             organism = 'mmu',
                                             db = mmu.annot,
                                             pvalueCutoff = list(enrichKEGG = 0.05,
                                                                 enrichMKEGG = 0.05),
                                             qvalueCutoff = 0.05,
                                             cnet.showCategory = 10,
                                             verbose = T,
                                             save.path = "KEGG",
                                             names = "neutrophil.R848.markers")
      R848.markers.go <- FastGO(genes = up_genes,
                                         geneList = geneList,
                                         organism="mouse",
                                         default.universe = F,
                                         classlevel = 2:2,
                                         OrgDb  = "org.Mm.eg.db",
                                         keyType = NULL,
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff  = 0.05,
                                         cnet.showCategory = 5,
                                         verbose = TRUE,
                                         save.path = "GO",
                                         names = "neutrophil.R848.markers")
      R848.markers.ms <- FastMS(genes=up_genes,
                                         geneList=geneList,
                                         default.universe = F,
                                         pAdjustMethod = "BH",
                                         db = mmu.annot,
                                         organism = "mouse",
                                         pvalueCutoff = list(enrichMS=0.05,
                                                             gseMS = 0.05),
                                         qvalueCutoff  = 0.05,
                                         cnet.showCategory = 5,
                                         verbose = TRUE,
                                         save.path = "Molecular Signatures",
                                         names = "neutrophil.R848.markers")
      object <- vehicle.Sofa.sig.genes.kegg$enrichKEGG
      FastEnrichPlot(object,
                     pvalueCutoff = NULL,
                     qvalueCutoff = NULL,
                     id.col = "Description",
                     select = NULL,
                     x.title = "Kegg pathway",
                     y.title = "-log10(FDR)",
                     size = 15,
                     short.cutoff = 46,
                     topshow = 10)
      
      #up_genes <-  as.character(na.omit(as.character(R848.SR3.sig.genes$X1[R848.SR3.sig.genes$lab=="UP"])))
      up_genes <- as.character(na.omit(convert(SR3.markers1, fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
      SR3.markers.kegg <- FastKEGG(genes = up_genes,
                                           geneList = geneList,
                                           default.universe = F,
                                           organism = 'mmu',
                                           db = mmu.annot,
                                           pvalueCutoff = list(enrichKEGG = 0.05,
                                                               enrichMKEGG = 0.05),
                                           qvalueCutoff = 0.05,
                                           cnet.showCategory = 10,
                                           verbose = T,
                                           save.path = "KEGG",
                                           names = "neutrophil.SR3.markers")
      SR3.markers.go <- FastGO(genes = up_genes,
                                       geneList = geneList,
                                       organism="mouse",
                                       default.universe = F,
                                       classlevel = 2:2,
                                       OrgDb  = "org.Mm.eg.db",
                                       keyType = NULL,
                                       pAdjustMethod = "BH",
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff  = 0.05,
                                       cnet.showCategory = 5,
                                       verbose = TRUE,
                                       save.path = "GO",
                                       names = "neutrophil.SR3.markers")
      SR3.markers.ms <- FastMS(genes=up_genes,
                                       geneList=geneList,
                                       default.universe = F,
                                       pAdjustMethod = "BH",
                                       db = mmu.annot,
                                       organism = "mouse",
                                       pvalueCutoff = list(enrichMS=0.05,
                                                           gseMS = 0.05),
                                       qvalueCutoff  = 0.05,
                                       cnet.showCategory = 5,
                                       verbose = TRUE,
                                       save.path = "Molecular Signatures",
                                       names = "neutrophil.SR3.markers")
      
      save(
        Vehicle.markers.kegg,Vehicle.markers.go,Vehicle.markers.ms,
        Sofa.markers.kegg,Sofa.markers.go,Sofa.markers.ms,
        R848.markers.kegg,R848.markers.go,R848.markers.ms,
        SR3.markers.kegg,SR3.markers.go,SR3.markers.ms
        ,file = "./Results/mouse/neutrophil.functionEnrichment.RData")
      #vehicle.Sofa.sig.genes.kegg$enrichKEGG@result %>% head()
      kegg1 <- Vehicle.markers.kegg$enrichKEGG@result %>% 
        mutate(sample=rep("Vehicle",length(ID))) %>% 
        filter(.,qvalue < 0.05) 
      kegg2 <- Sofa.markers.kegg$enrichKEGG@result %>% 
        mutate(sample=rep("Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05)
      kegg3 <- R848.markers.kegg$enrichKEGG@result %>% 
        mutate(sample=rep("R848",length(ID))) %>% 
        filter(.,qvalue < 0.05)
      kegg4 <- SR3.markers.kegg$enrichKEGG@result %>% 
        mutate(sample=rep("SR3",length(ID))) %>% 
        filter(.,qvalue < 0.05)
      
      kegg <- kegg1 %>% 
        bind_rows(.,kegg2) %>% 
        bind_rows(.,kegg3) %>% 
        bind_rows(.,kegg4) 
      head(kegg)      
      ggplot(data=kegg,mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
        geom_point() +
        theme_cowplot(font_size = 12) +
        labs(x="sample")
      
      
      go_bp1 <- Vehicle.markers.go$GOEnrichment$BP@result %>% 
        mutate(sample=rep("Vehicle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_bp2 <- Sofa.markers.go$GOEnrichment$BP@result %>% 
        mutate(sample=rep("Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_bp3 <- R848.markers.go$GOEnrichment$BP@result %>% 
        mutate(sample=rep("R848",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_bp4 <- SR3.markers.go$GOEnrichment$BP@result %>% 
        mutate(sample=rep("SR3",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      
      go_bp <- go_bp1 %>% 
        bind_rows(.,go_bp2) %>% 
        bind_rows(.,go_bp3) %>% 
        bind_rows(.,go_bp4) 
      head(go_bp)      
      ggplot(data=go_bp,mapping=aes(x=factor(sample,levels = c("Vehicle","Sofa","R848","SR3")),y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
        geom_point() +
        theme_cowplot(font_size = 12) +
        scale_x_discrete(limits = c("Vehicle","Sofa","R848","SR3"),labels = c("Vehicle","Sofa","R848","SR3")) +
        labs(x="sample")
      
      go_mf1 <- vehicle.Sofa.sig.genes.go$GOEnrichment$MF@result %>% 
        mutate(sample=rep("Sofa vs Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_mf2 <- vehicle.R848.sig.genes.go$GOEnrichment$MF@result %>% 
        mutate(sample=rep("R848 vs Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_mf3 <- vehicle.SR3.sig.genes.go$GOEnrichment$MF@result %>% 
        mutate(sample=rep("SR3 vs Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_mf4 <- Sofa.R848.sig.genes.go$GOEnrichment$MF@result %>% 
        mutate(sample=rep("R848 vs Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_mf5 <- Sofa.SR3.sig.genes.go$GOEnrichment$MF@result %>% 
        mutate(sample=rep("SR3 vs Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_mf6 <-  R848.SR3.sig.genes.go$GOEnrichment$MF@result %>% 
        mutate(sample=rep("SR3 vs R848",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_mf <- go_mf2 %>% 
        bind_rows(.,go_mf3) %>% 
        bind_rows(.,go_mf4) %>% 
        bind_rows(.,go_mf5) %>% 
        bind_rows(.,go_mf6) 
      head(go_mf)      
      ggplot(data=go_mf,mapping=aes(x=sample,y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
        geom_point() +
        theme_cowplot(font_size = 12)
      
      go_cc1 <- vehicle.Sofa.sig.genes.go$GOEnrichment$CC@result %>% 
        mutate(sample=rep("Sofa vs Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_cc2 <- vehicle.R848.sig.genes.go$GOEnrichment$CC@result %>% 
        mutate(sample=rep("R848 vs Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_cc3 <- vehicle.SR3.sig.genes.go$GOEnrichment$CC@result %>% 
        mutate(sample=rep("SR3 vs Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_cc4 <- Sofa.R848.sig.genes.go$GOEnrichment$CC@result %>% 
        mutate(sample=rep("R848 vs Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_cc5 <- Sofa.SR3.sig.genes.go$GOEnrichment$CC@result %>% 
        mutate(sample=rep("SR3 vs Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_cc6 <-  R848.SR3.sig.genes.go$GOEnrichment$CC@result %>% 
        mutate(sample=rep("SR3 vs R848",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      go_cc <- go_cc2  %>% 
        bind_rows(.,go_cc3) %>% 
        bind_rows(.,go_cc4) %>% 
        bind_rows(.,go_cc5) %>% 
        bind_rows(.,go_cc6) 
      head(go_cc)  
      
      ggplot(data=go_cc,mapping=aes(x=sample,y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
        geom_point() +
        theme_cowplot(font_size = 12)
      
      hallmark1 <- Vehicle.markers.ms$Data$enrichMS$h@result %>% 
        mutate(sample=rep("Vehichle",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      hallmark2 <- Sofa.markers.ms$Data$enrichMS$h@result %>% 
        mutate(sample=rep("Sofa",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      hallmark3 <- R848.markers.ms$Data$enrichMS$h@result %>% 
        mutate(sample=rep("R848",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      hallmark4 <- SR3.markers.ms$Data$enrichMS$h@result %>% 
        mutate(sample=rep("SR3",length(ID))) %>% 
        filter(.,qvalue < 0.05) %>% 
        slice_min(qvalue,n=10)
      
      hallmark <- hallmark1 %>% 
        bind_rows(.,hallmark2) %>% 
        bind_rows(.,hallmark3) %>% 
        bind_rows(.,hallmark4) 
      head(hallmark)  
      hallmark$sample <- factor(hallmark$sample,levels=c("Vehicle","Sofa","R848","SR3"))
      ggplot(data=hallmark,mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=Description,color=-log(as.numeric(qvalue)),size=parse_ratio(GeneRatio))) +
        geom_point() +
        theme_cowplot(font_size = 12) +
        scale_x_discrete(limits = c("Vehicle","Sofa","R848","SR3"),labels = c("Vehicle","Sofa","R848","SR3")) +
        labs(x="sample")
    }
  }
  
  #M1和M2分析
  {
    library(CytoTRACE)
    library(monocle3)
    Marcrophage_sce <- readRDS(file="./Results/mouse/Marcrophage.rds")
    Marcrophage_sce$M_type <- ifelse(Marcrophage_sce$seurat_clusters==1,"M1","M2")
    p1 <- dittoDimPlot(Marcrophage_sce, reduction.use  = "umap", var="seurat_clusters", size = 1, do.label = TRUE, 
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "Marcrophage") + theme(legend.text = element_text(face="bold",size = 12))
    p2 <- dittoDimPlot(Marcrophage_sce, reduction.use  = "umap", var="orig.ident", size = 1, do.label = TRUE, 
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "myeloid cells") + theme(legend.text = element_text(face="bold",size = 12))
    p3 <- dittoDimPlot(Marcrophage_sce, reduction.use  = "umap", var="M_type", size = 1, do.label = TRUE, 
                       do.ellipse = TRUE, legend.size = 5, shape.legend.size=5, labels.size=4, do.raster=T, raster.dpi = 500,legend.show = T, main = "myeloid cells") + theme(legend.text = element_text(face="bold",size = 12))
    p1 | p2 | p3
    #Cytotrace
    count_mat <-  as.matrix(GetAssayData(Marcrophage_sce, slot='counts',assay = "RNA"))
    results <- CytoTRACE(count_mat,ncores = 10)
    # 使用自己的umap
    cord <- Marcrophage_sce@reductions$umap@cell.embeddings
    pheno <- as.character(Marcrophage_sce$M_type)
    names(pheno) <- colnames(Marcrophage_sce)
    plotCytoTRACE(results, phenotype = pheno, emb=cord[names(pheno),])
    Marcrophage_sce[['cytotrace_score']] = results$CytoTRACE[colnames(Marcrophage_sce)]
    FeaturePlot(Marcrophage_sce, features = 'cytotrace_score') + 
      scale_colour_gradientn(colours = c("lightgrey","#3288BD" ,"#ABDDA4","#FEE08B", "#F46D43","#9E0142"))
    plotCytoGenes(results, numOfGenes = 20)
    saveRDS(results,file="/Users/biofly/project/shijian/Results/mouse/Marcrophage_CytoTrace.rds")
    results <- readRDS(file="/Users/biofly/project/shijian/Results/mouse//Marcrophage_CytoTrace.rds")
    #monocle3
    cds <- SeuratToMonocle3(seurat_object=Marcrophage_sce,scale_all = FALSE,assay = "RNA",slot="counts",reduction_for_projection = "pca",UMAP_cluster_slot = NULL)
    cds <- learn_graph(cds)
    cds <- order_cells(cds, reduction_method = "UMAP")
    p <- plot_cells(cds = cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE,label_leaves=F,label_branch_points = F)
    ggsave(p,filename = "figure6_monocle3.pdf")
    plot_genes_in_pseudotime(cds,color_cells_by="pseudotime",min_expr=NULL)
    
    # 评估M1细胞抗原呈递能力
    genes <- c("HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DRB1","HLA-DRB5","HLA-DRA","HLA-A","HLA-B","HLA-C") #人的
    convertMouseGeneList1(genes)
    convert_human_to_mouse_symbols(genes)
    # 小鼠
    M1_seurat <- subset(Marcrophage_sce,M_type=="M1")
    genes <- c("H2-D1","H2-K1","H2-L","H-2Q","H2","H2-T","H2-Eb1","H2-Ea","H2-Bl")
    a <- FetchData(M1_seurat,vars = genes,slot = "data")  #"H2-D1,H2-K1,H2-Eb1"
    a <- a %>% tibble::rownames_to_column("cell") 
    a$sample <- stringr::str_split(a$cell,"_",simplify = T)[,1]
    head(a)
    a %>% 
      dplyr::group_by(sample) %>% 
      dplyr::summarise(`H2-D1`= mean(`H2-D1`),`H2-K1`=mean(`H2-K1`),`H2-Eb1`=mean(`H2-Eb1`)) %>% 
      tibble::column_to_rownames("sample") -> heatmap.data
    pheatmap::pheatmap(t(heatmap.data[c("Vehicle","Sofa","R848","SR3"),]),scale="row",cluster_rows = F,cluster_cols = F,filename = "figure6_heatmap.pdf")
    #M2 vs M1差异基因
    diff_genes <- FindMarkers(Marcrophage_sce, ident.1 ="M1", ident.2=NULL,group.by = "M_type")
    M2_pos_genes <- FindMarkers(Marcrophage_sce,ident.1 ="M2",ident.2=NULL,group.by = "M_type",only.pos = T)
    write.table(M2_pos_genes,file="./Results/data/mouse_M2_pos_genes.txt",sep="\t",quote=F,row.names=T,col.names=F)
    mmu.annot <- getAnnotation(gtf.path = "/Users/biofly/project/shijian/gencode.vM28.annotation.gtf.gz")
    mmu.annot <- getAnnotation2(mmu.annot,organism = "mouse")
    up_genes <- as.character(na.omit(convert(rownames(M2_pos_genes),fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)))
    geneList <- Marcrophage_sce@assays$RNA@data[,1]
    names(geneList) <- convert(names(geneList), fromtype = "SYMBOL",totype = "ENSEMBL",db=mmu.annot)
    geneList <- geneList[!is.na(names(geneList))]
    kegg_result <- FastKEGG(genes = up_genes,
                                     geneList = geneList,
                                     default.universe = F,
                                     organism = 'mmu',
                                     db = mmu.annot,
                                     pvalueCutoff = list(enrichKEGG = 0.05,
                                                         enrichMKEGG = 0.05),
                                     qvalueCutoff = 0.05,
                                     cnet.showCategory = 10,
                                     verbose = T,
                                     save.path = "KEGG",
                                     names = "M2_M1_upGenes")
    go_results <- FastGO(genes = up_genes,
                                 geneList = geneList,
                                 organism="mouse",
                                 default.universe = F,
                                 classlevel = 2:2,
                                 OrgDb  = "org.Mm.eg.db",
                                 keyType = NULL,
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff  = 0.05,
                                 cnet.showCategory = 5,
                                 verbose = TRUE,
                                 save.path = "GO",
                                 names = "M2_M1_upGenes")
     ms_results <- FastMS(genes=up_genes,
                           geneList=geneList,
                           default.universe = F,
                           pAdjustMethod = "BH",
                           db = mmu.annot,
                           organism = "mouse",
                           pvalueCutoff = list(enrichMS=0.05,
                                               gseMS = 0.05),
                           qvalueCutoff  = 0.05,
                           cnet.showCategory = 5,
                           verbose = TRUE,
                           save.path = "Molecular Signatures",
                           names = "M2_M1_upGenes")
     
     #小鼠数据M2上调基因的功能注释，在加上REACTOME、PID、BIOCARTA都做一下
     #Reactome enrichment analysis
     library(ReactomePA)
     up_genes_id <- convert(up_genes,fromtype = "ENSEMBL",totype = "ENTREZID",db=mmu.annot)
     up_genes_id <- up_genes_id[!is.na(up_genes_id)]
     reactome_results <- enrichPathway(gene=up_genes_id, organism = "mouse",pvalueCutoff = 0.05, qvalueCutoff = 0.05,readable=TRUE)
     dotplot(reactome_results)
     #biocarta和PID富集分析
     names(geneList) <- convert(names(geneList), fromtype = "ENSEMBL", totype = "ENTREZID", db=mmu.annot)
     geneList <- geneList[!is.na(names(geneList))]
     msigdf.mouse %>% filter(.,category_subcode=="cp.biocarta") %>% select(.,geneset,mouse.symbol) -> df_biocarta
     msigdf.mouse %>% filter(.,category_subcode=="cp.pid") %>% select(.,geneset,mouse.symbol) -> df_pid
     geneList <- Marcrophage_sce@assays$RNA@data[,1]
     biocartaObject <- enricher(rownames(M2_pos_genes),universe = names(geneList), pvalueCutoff = 0.05, qvalueCutoff  = 0.05, pAdjustMethod = "BH", TERM2GENE = df_biocarta)
     pidObject <- enricher(rownames(M2_pos_genes),universe = names(geneList), pvalueCutoff = 0.05, qvalueCutoff  = 0.05, pAdjustMethod = "BH", TERM2GENE = df_pid)
     dotplot(biocartaObject)
     dotplot(pidObject)
  }
  C2 <- fread("/Users/biofly/project/shijian/M2_M1_upGenes/Molecular Signatures/M2_M1_upGenes_c2_Molecular Signatures enrichment.csv")
  reactom_select <- c("REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION","REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
  "REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS","REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA","REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE",
  "REACTOME_PROTEIN_LOCALIZATION","REACTOME_INFECTION_WITH_MYCOBACTERIUM_TUBERCULOSIS","REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES","REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION",
  "REACTOME_METABOLISM_OF_CARBOHYDRATES","REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2",
  "REACTOME_ADAPTIVE_IMMUNE_SYSTEM","REACTOME_SIGNALING_BY_INTERLEUKINS",
  "REACTOME_HEMOSTASIS","REACTOME_IRON_UPTAKE_AND_TRANSPORT","REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
  "REACTOME_TRANSPORT_OF_SMALL_MOLECULES","REACTOME_CELLULAR_RESPONSES_TO_STIMULI","REACTOME_NEUTROPHIL_DEGRANULATION","REACTOME_INNATE_IMMUNE_SYSTEM")
  p <- C2[Fastgrep2(reactom_select,C2$V1),] %>% na.omit() %>% as.data.frame() %>% mutate(qq=-log(p.adjust)) %>% arrange(.,as.numeric(qq)) %>% 
    ggplot(mapping=aes(x=factor(V1,levels=V1),y=qq)) +
      geom_bar(stat = 'identity',fill="#53903b") + 
      coord_flip() +
      labs(y="-log(p.adjusted)",x="REATOME pathway") + 
      theme_bw()
  ggsave(p,filename = "./code/tianjin_Supple/singlecell_mouse/pictures/reactome_results.pdf",width = 14,height=7)
}

