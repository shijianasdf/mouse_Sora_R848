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
