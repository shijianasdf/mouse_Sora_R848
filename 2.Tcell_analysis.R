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
