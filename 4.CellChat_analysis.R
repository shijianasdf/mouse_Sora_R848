#' @description cellchat analysis pipeline
#' @author  shi jian
sce_merge_res1 <- readRDS("./Results/data/mouse_sce_merge_res1.rds")
sce_merge_res1$sce@meta.data %>% head()
table(sce_merge_res1$sce@meta.data$celltype)
p1 <- dittoDimPlot(sce_merge_res1$sce, reduction.use  = "umap",var="celltype",size = 1,do.label = TRUE, labels.highlight = F,do.ellipse = F,legend.size = 9,shape.legend.size=9,labels.size=7,do.raster=T,raster.dpi = 500) + 
            theme(legend.text = element_text(face="bold",size = 16))
sce_merge_res1$sce@meta.data$orig.ident %>% unique()
Vehicle.seurat <- subset(sce_merge_res1$sce,orig.ident=="Vehicle")
R848.seurat <- subset(sce_merge_res1$sce,orig.ident=="R848")
Sofa.seurat <- subset(sce_merge_res1$sce,orig.ident=="Sofa")
SR3.seurat <- subset(sce_merge_res1$sce,orig.ident=="SR3")
Vehicle.cellchat <- FastCellChat(data.input=Vehicle.seurat@assays$RNA@data,meta=Vehicle.seurat@meta.data,
             species=c("Human","Mouse")[2],group.by="celltype",search=NULL,project=F)
R848.cellchat <- FastCellChat(data.input=R848.seurat@assays$RNA@data,meta=R848.seurat@meta.data,
             species=c("Human","Mouse")[2],group.by="celltype",search=NULL,project=F)
Sofa.cellchat <- FastCellChat(data.input=Sofa.seurat@assays$RNA@data,meta=Sofa.seurat@meta.data,
             species=c("Human","Mouse")[2],group.by="celltype",search=NULL,project=F)
SR3.cellchat <- FastCellChat(data.input=SR3.seurat@assays$RNA@data,meta=SR3.seurat@meta.data,
             species=c("Human","Mouse")[2],group.by="celltype",search=NULL,project=F)
save(Vehicle.cellchat,R848.cellchat,Sofa.cellchat,SR3.cellchat,file="./code/tianjin_Supple/cellchat_result.rda")

object.list <- list(Vehicle = Vehicle.cellchat, Sofa=Sofa.cellchat, R848 = R848.cellchat, SR3=SR3.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
df.net <- subsetCommunication(cellchat)
df.net$Vehicle %>% filter(.,source=="Neutrophil_cell") %>% as.data.frame() %>% 
                   dplyr::group_by(source,target) %>% 
                   dplyr::summarise(num=n()) %>% 
                   mutate(sample=rep("Vehicle",length(source))) -> a1
df.net$Sofa %>% filter(.,source=="Neutrophil_cell") %>% as.data.frame() %>% 
  dplyr::group_by(source,target) %>% 
  dplyr::summarise(num=n()) %>% 
  mutate(sample=rep("Sofa",length(source))) -> a2
df.net$R848 %>% filter(.,source=="Neutrophil_cell") %>% as.data.frame() %>% 
  dplyr::group_by(source,target) %>% 
  dplyr::summarise(num=n()) %>% 
  mutate(sample=rep("R848",length(source))) -> a3
df.net$SR3 %>% filter(.,source=="Neutrophil_cell") %>% as.data.frame() %>% 
  dplyr::group_by(source,target) %>% 
  dplyr::summarise(num=n()) %>% 
  mutate(sample=rep("SR3",length(source))) -> a4
a <- rbind.data.frame(a1,a2,a3,a4)
a <- as.data.frame(a)
head(a)
library(igraph)
library(ggraph)
library(tidygraph)
library(ktplots)
library(ggrepel)
colnames(a) <- c("from","to","num","sample")
aa.list <- split(a,a$to)
plot.data <- do.call(rbind.data.frame,lapply(aa.list,function(x){
  pos <- which.max(x$num)
  x[pos,]
}))
# plot using ggraph
levels <- c("CD4+T_cell","CD8+T_cell","Dendritic_cell","Marcrophage",
            "B_cell","Epithelial_cell","Mast_cell","NK","Neutrophil_cell",
            "Plasma","Stromal_cell")
graph <- as_tbl_graph(plot.data) %>% 
  mutate(Popularity = centrality_degree(mode = 'out'))
graph <- permute(graph,match(V(graph)$name,levels))
node_group_colors = rep("black",11)#dittoColors()[1:11]
edge_group_colors = c(Vehicle="#51855b",Sofa="#b8dcb1",R848="#f4c289",SR3="#ec8633")
p <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_arc(aes(colour = sample),arrow = arrow(length = unit(4, 'mm')),width=2,alpha=1) +
  scale_edge_color_manual(name="Sample",values = edge_group_colors) + 
  geom_node_point(pch = 19, size=8,alpha=0.8) +
  theme_void() + coord_fixed() + scale_size_continuous(limits = c(0,1)) + 
  scale_shape_manual(values = c(ligand = 19, receptor = 15)) +
  scale_color_manual(values = node_group_colors) +
  geom_text_repel(aes(x = x, y = y, label = name), segment.square = TRUE, segment.inflect = TRUE,
                  segment.size = 0.2, force = 0.5, size = 4, force_pull = 0)
ggsave(p,filename = "./cellchat_network.pdf",width=9,height=9)
