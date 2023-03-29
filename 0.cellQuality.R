# @title 每个样本的质量控制、去双细胞、细胞类型注释、以及结果的展示
# @author shi jian
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
