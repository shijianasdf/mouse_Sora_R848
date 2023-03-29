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
    FeaturePlot(Marcrophage_sce, features = 'cytotrace_score') + scale_colour_gradientn(colours = c("lightgrey","#3288BD" ,"#ABDDA4","#FEE08B", "#F46D43","#9E0142"))
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
