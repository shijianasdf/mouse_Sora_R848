# @author shijian  
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
  
