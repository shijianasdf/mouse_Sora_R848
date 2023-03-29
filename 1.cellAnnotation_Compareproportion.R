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
