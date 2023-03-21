#' @description  run cellphonedb cellphonedb目前只支持人类，所以要把小鼠基因转成人类基因
#' @author shi jian 
library(dplyr)
library(ggsignif)
library(data.table)
library(Seurat)
library(ggplot2)
library(tidytext)
library(forcats)
library(stringr)
library(nichenetr)
library(magrittr)
#Sometimes the main Ensembl BioMart server can break, which leads to this error. 
#It's not a problem with your code, and hopefully the Ensembl team will fix the issue shortly.
sce_merge_res1 <- readRDS("./Results/data/mouse_sce_merge_res1.rds")
sce_merge_res1$sce@meta.data %>% head()
table(sce_merge_res1$sce@meta.data$celltype)
p1 <- dittoDimPlot(sce_merge_res1$sce, reduction.use  = "umap",var="celltype",size = 1,do.label = TRUE, labels.highlight = F,do.ellipse = F,legend.size = 9,shape.legend.size=9,labels.size=7,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 16))
p1
##########注意cellphonedb是一个样本一个样本的跑##########
#制作每个样本的细胞注释类型文件
cellAnno_outDir <- "code/tianjin_Supple/singlecell_mouse/cellAnno_outDir"
if(!file.exists(cellAnno_outDir)){dir.create(cellAnno_outDir,recursive = T)}
metadata <- sce_merge_res1$sce@meta.data %>% tibble::rownames_to_column("cell") %>% dplyr::select(.,orig.ident,cell,celltype)
samples_cellAnno <- split(metadata,f=metadata$orig.ident)
lapply(samples_cellAnno,function(x){
  name <- x %>% dplyr::select(.,orig.ident) %>% distinct()
  temp <- x %>% dplyr::select(.,c(cell,celltype))
  write.table(temp,file = paste0("./code/tianjin_Supple/singlecell_mouse/cellAnno_outDir/cellphonedb_metadata_",name,".txt"),sep="\t",quote=F,row.names=F,col.names=T)
  return(temp)
})
#制作表达文件
makeCellPhoneExpr <- function(mouse.seurat,sample,outDir){
  expr.matrix <- mouse.seurat@assays$RNA@data
  expr.matrix %>% as.data.frame() %>% tibble::rownames_to_column("Gene") -> expr.matrix
  musGenes <- expr.matrix$Gene
  #小鼠基因转成人类基因
  humGenes <- convertMouseGeneList2(musGenes) 
  pos <- which(!is.na(humGenes$Homo.sapiens))
  expr.matrix1 <- expr.matrix[pos,] 
  expr.matrix1$Gene <- unlist(humGenes[pos,]$HUM)
  write.table(expr.matrix1,file = paste0(outDir,"/cellphonedb_expr_",sample,".txt"),sep="\t",quote=F,row.names=F,col.names=T)
  return(expr.matrix1)
}
cellPhone_exprDir <- "code/tianjin_Supple/singlecell_mouse/cellPhone_exprDir"
if(!file.exists(cellPhone_exprDir)){dir.create(cellPhone_exprDir,recursive = T)}
samples <- unique(sce_merge_res1$sce@meta.data$orig.ident)
for(i in 1:length(samples)){
  subset_seurat <- subset(sce_merge_res1$sce,orig.ident==samples[i])
  makeCellPhoneExpr(subset_seurat,samples[i],cellPhone_exprDir)
}
#运行cellphonedb程序
cd /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse
mkdir cellphonedb_output
cd cellphonedb_output
mkdir R848
cellphonedb method statistical_analysis --counts-data gene_name \
          --output-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/R848 \
          /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellAnno_outDir/cellphonedb_metadata_R848.txt \
          /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellPhone_exprDir/cellphonedb_expr_R848.txt \
          --threads 10
#气泡图
cellphonedb plot dot_plot --pvalues-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/R848/pvalues.txt \
                  --means-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/R848/means.txt \
                  --output-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/R848 \
                  --output-name R848.dotplot.pdf
#热图
cellphonedb plot heatmap_plot --pvalues-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/R848/pvalues.txt \
                              --output-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/R848 \
                              --pvalue 0.05 \
                              --count-name R848.heatmap_count.pdf \
                              --log-name R848.heatmap_log_count.pdf \
                              --count-network-name R848.count_network.txt \
                              --interaction-count-name R848.interaction_count.txt /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellAnno_outDir/cellphonedb_metadata_R848.txt 
cd /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output
mkdir Sofa
cellphonedb method statistical_analysis --counts-data gene_name \
--output-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/Sofa \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellAnno_outDir/cellphonedb_metadata_Sofa.txt \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellPhone_exprDir/cellphonedb_expr_Sofa.txt \
--threads 10

cd /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output
mkdir Vehicle
cd Vehicle
cellphonedb method statistical_analysis --counts-data gene_name \
--output-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/Vehicle \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellAnno_outDir/cellphonedb_metadata_Vehicle.txt \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellPhone_exprDir/cellphonedb_expr_Vehicle.txt \
--threads 10

cd /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output
mkdir SR3
cd SR3
cellphonedb method statistical_analysis --counts-data gene_name \
--output-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/SR3 \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellAnno_outDir/cellphonedb_metadata_SR3.txt \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellPhone_exprDir/cellphonedb_expr_SR3.txt \
--threads 10
# output
# deconvoluted.txt：每对受体-配体提供了附加信息
# mean.txt：每对受体-配体相互作用的表达平均值，所有相互作用伙伴的平均值：平均值是指相应相互作用的细胞类型对中个体伙伴平均表达值的总平均值。
#           如果其中一个平均值为 0，则总平均值设置为 0。
# pvalues.txt: 每对受体-配体的P值，p 值是指在每个相互作用的细胞类型对中相互作用的配体-受体对的富集
# significant_means.txt: 每对受体-配体显著性结果的平均表达量值。所有交互伙伴的显著平均计算。如果 p.value < 0.05，则该值将是平均值。或者，该值设置为 0。
#重要的是要记住，相互作用不是对称的。换言之，当测试簇 X_Y 之间的配体/受体对 A_B 时，伙伴 A 的表达被认为在第一个簇 (X) 内，而伙伴 B 在第二个簇 (Y) 内的表达。
#因此，X_Y 和 Y_X 代表不同的比较，将具有不同的 p 值和均值。
Vehicle_cellphone_pvalues <- fread("./code/tianjin_Supple/singlecell_mouse/cellphonedb_output/Vehicle/pvalues.txt")
Sofa_cellphone_pvalues <- fread("./code/tianjin_Supple/singlecell_mouse/cellphonedb_output/Sofa/pvalues.txt")
R848_cellphone_pvalues <- fread("./code/tianjin_Supple/singlecell_mouse/cellphonedb_output/R848/pvalues.txt")
SR3_cellphone_pvalues <- fread("./code/tianjin_Supple/singlecell_mouse/cellphonedb_output/SR3/pvalues.txt")
#plot.data <- data.frame(sample=character(0),celltype_=character(0),celltype_interact=character(0),direction=character(0),num=numeric(0))
makeplotdata <- function(x,cell_type,sample){
  Vehicle_cellphone_pvalues <- x
  Vehicle_cellphone_pvalues <- as.data.frame(Vehicle_cellphone_pvalues)
  pos <- grep(cell_type,colnames(Vehicle_cellphone_pvalues))
  Vehicle_cellphone <- Vehicle_cellphone_pvalues[ ,colnames(Vehicle_cellphone_pvalues)[pos]]
  # a <- Vehicle_cellphone_pvalues[ ,colnames(Vehicle_cellphone_pvalues)[pos1]]
  # b <- Vehicle_cellphone_pvalues[ ,colnames(Vehicle_cellphone_pvalues)[pos2]]
  rownames(Vehicle_cellphone) <- Vehicle_cellphone_pvalues$interacting_pair
  Vehicle_cellphone <- ifelse(Vehicle_cellphone < 0.05,1,0)
  n1 <- colSums(Vehicle_cellphone)
  p <- data.frame(celltype_pair=names(n1),num=n1)
  p$direction <- ifelse(grepl(paste0(cell_type,"\\|"),p$celltype_pair),"sender","recevier")
  d <- stringr::str_split(p$celltype_pair,"\\|",simplify = T)
  p$celltype <- NA
  p$celltype[!grepl(cell_type,d[,1])] <- d[,1][!grepl(cell_type,d[,1])]
  p$celltype[!grepl(cell_type,d[,2])] <- d[,2][!grepl(cell_type,d[,2])]
  p$celltype[d[,1]==cell_type & d[,2]==cell_type] <- cell_type
  h <- p[which(p$celltype_pair == paste0(cell_type,"|",cell_type)),]
  h$direction <- "recevier"
  p <- rbind(p,h)
  p %>% dplyr::group_by(celltype) %>% dplyr::summarise(num=sum(num)) %>% mutate(direction=rep("both",length(celltype))) -> l
  l$num[which(l$celltype == cell_type)] <- l$num[which(l$celltype == cell_type)]/2
  p %>% dplyr::select(.,celltype,num,direction) %>% bind_rows(l) %>% mutate(sample=rep(sample,length(celltype))) -> plot.data
  return(plot.data)
}
a1 <- makeplotdata(Vehicle_cellphone_pvalues,"Neutrophil_cell","Vehicle")
a2 <- makeplotdata(Sofa_cellphone_pvalues,"Neutrophil_cell","Sofa")
a3 <- makeplotdata(R848_cellphone_pvalues,"Neutrophil_cell","R848")
a4 <- makeplotdata(SR3_cellphone_pvalues,"Neutrophil_cell","SR3")
a <- rbind.data.frame(a1,a2,a3,a4)
a
#只画sender那一行，互作最多的用药组，网络图展示
library(igraph)
library(ggraph)
library(tidygraph)
library(ktplots)
library(ggrepel)
a %>% dplyr::filter(.,direction == "sender") %>%   #取最多互作的样本可视化
      mutate(from=rep("Neutrophil_cell",length(celltype))) %>% 
      dplyr::rename(.,to=celltype) %>% 
      select(.,c(5,1,2,3,4)) -> aa 
aa.list <- split(aa,aa$to)
plot.data <- do.call(rbind.data.frame,lapply(aa.list,function(x){
  pos <- which.max(x$num)
  x[pos,]
}))
# plot using ggraph
graph <- as_tbl_graph(plot.data) %>% 
  mutate(Popularity = centrality_degree(mode = 'out'))
levels <- c("CD4+T_cell","CD8+T_cell","Dendritic_cell","Marcrophage",
            "B_cell","Mast_cell","NK","Neutrophil_cell",
            "Plasma","Stromal_cell","Epithelial_cell")
graph <- permute(graph,match(V(graph)$name,levels))
# .gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# nn = length(unique(E(graph)$sample))
# edge_group_colors = .gg_color_hue(nn)
# nnode = length(unique(V(graph)$name))
# node_group_colors = .gg_color_hue(nnode)
# edge_group_colors = c("Activating" = "#e15759", "Chemotaxis" = "#59a14f", "Inhibitory" = "#4e79a7", "   Intracellular trafficking" = "#9c755f", "DC_development" = "#B07aa1")
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
ggsave(p,filename = "./cellphonedb_network.pdf",width=9,height=9)
b1 <- makeplotdata(Vehicle_cellphone_pvalues,"Marcrophage","Vehicle")
b2 <- makeplotdata(Sofa_cellphone_pvalues,"Marcrophage","Sofa")
b3 <- makeplotdata(R848_cellphone_pvalues,"Marcrophage","R848")
b4 <- makeplotdata(SR3_cellphone_pvalues,"Marcrophage","SR3")
b <- rbind.data.frame(b1,b2,b3,b4)
b
library(ggrepel)
p1 <- ggplot(data=a,mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=num,fill=sample)) + 
  geom_bar(stat="identity",position=position_dodge()) +
  facet_grid(celltype~direction,scales = "free_y") +
  scale_fill_manual(values = dittoColors()) + 
  geom_text(mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=num+2,label=num)) +
  #scale_y_continuous(expand = c(0,1)) +
  labs(x="sample",y="interaction numbers",title = "Neutrophil") +
  theme_bw() +
  theme(axis.ticks.length=unit(0.5,'cm'),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size=20)) + 
  theme(strip.text = element_text(size = 8,face = "bold"),panel.margin = unit(0, "lines"),strip.placement = "outside") +
  coord_flip()
ggsave(p1,filename = "figure7_neutrophil_cellphonedb.pdf",width=16,height=12)
print(p1)
#pie图展示
#生成待可视化的数据
tt <- Marcrophage_sce@meta.data %>% 
  mutate(seurat_clusters=paste0("cluster",as.character(seurat_clusters))) %>% 
  dplyr::group_by(seurat_clusters) %>% 
  dplyr::summarise(proportion=table(orig.ident)/length(orig.ident),names=names(table(orig.ident))) %>% 
  mutate(proportion=as.numeric(proportion)) %>% 
  mutate(labels=scales::percent(proportion)) 
plot.data  <- data.frame(sample=character(0),proportion=numeric(0),alteration=character(0),
                         celltype=character(0),direction=character(0))    
a %>% head()
plot.list <- split(a,list(a$celltype,a$direction))
plot.list <- lapply(plot.list,function(x){
                #x <- plot.list[[1]]
                pos <- which(x$sample == "Vehicle")
                x$num <- x$num - x$num[pos]
                x %<>% mutate(alteration=if_else(num == 0,"no",if_else(num > 0,"gain","loss"))) %>% 
                       mutate(num=abs(num)) %>% 
                       mutate(proportion=num/sum(num)) %>% 
                       mutate(label=scales::percent(proportion)) %>% 
                       mutate(label1=if_else(alteration=="no",label,
                                             if_else(alteration=="gain",paste0("+",label),paste0("-",label))))
                x
              })   
plot.data <- do.call(rbind.data.frame,plot.list)
plot.data %>% filter(.,sample!= "Vehicle",alteration != "no") -> plot.data
cols <- dittoColors()
scales::show_col(cols[1:20])
library(ggpattern) #给ggplot2图片加入阴影
p2 <- ggplot(data = plot.data, mapping=aes(x="",y=proportion,fill=sample)) +
        geom_col_pattern(aes(pattern=alteration),
                         pattern_alpha=0.5,
                         pattern_density=0.1,
                         pattern_colour="#A7727D", #9DC08B BAD7E9 lightblue
                         pattern_size=0.4) +
        geom_text(mapping=aes(label = label),position = position_stack(vjust = 0.5),size=3) +
        scale_fill_manual(values = c(SR3="#ff8000",R848="#ffc080",Sofa="#addead")) +
        scale_pattern_discrete(choices = c("stripe", "crosshatch")) +  #"crosshatch"
        #scale_pattern_density_discrete(range = c(0.01, 0.3)) +
        coord_polar(theta = "y") +
        facet_grid(direction~celltype) +
        theme_void()
ggsave(p2,filename="./multipie_plot.pdf",width = 14,height=10)
p2 <- ggplot(data = plot.data, mapping=aes(x="",y=proportion,fill=sample,color=alteration)) +
        geom_bar(stat="identity",width=0.1) +
        geom_text(mapping=aes(label = label),position = position_stack(vjust = 0.5),size=3) +
        scale_fill_manual(values = c(SR3="#ff8000",R848="#ffc080",Sofa="#addead")) +
        scale_color_manual(values=c(gain=cols[5],loss=cols[11])) + 
        #scale_fill_manual(values = c(gain=cols[3],loss=cols[9])) +
        #scale_color_manual(values=c(SR3="#ff8000",R848="#ffc080",Sofa="#addead")) + 
        coord_polar(theta = "y") +
        facet_grid(direction~celltype) +
        theme_void()
p1 <- ggplot(data=b,mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=num,fill=sample)) + 
        geom_bar(stat="identity",position=position_dodge()) +
        facet_grid(celltype~direction,scales = "free_y") +
        scale_fill_manual(values = dittoColors()) + 
        geom_text(mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=num+2,label=num)) +
        #scale_y_continuous(expand = c(0,1)) +
        labs(x="sample",y="interaction numbers",title = "Macrophage") +
        theme_bw() +
        theme(axis.ticks.length=unit(0.5,'cm'),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size=20)) + 
        theme(strip.text = element_text(size = 10,face = "bold"),panel.margin = unit(0, "lines"),strip.placement = "outside") +
        coord_flip()
print(p1)

#这个图中按照sender那一列，SR3互作数目与第二个互作最多的组做差，按照差异程度从高到低排序做一个柱状图
cellNum_diff <- function(a,celltype1){
  b <- a %>% filter(.,direction=="sender",celltype==celltype1)
  d <- b$num[which(b$sample == "SR3")]
  b1 <- b %>% filter(.,sample!="SR3") 
  s <- b1$num[which.max(b1$num)]     
  e <- d-s
  temp.df.i <- data.frame(celltype=celltype1,diff=e)
  return(temp.df.i)
}
do.call(rbind.data.frame,lapply(unique(a$celltype),function(x){cellNum_diff(a,x)})) %>% 
  arrange(.,desc(diff)) -> order
p <- do.call(rbind.data.frame,lapply(unique(a$celltype),function(x){cellNum_diff(a,x)})) %>% 
  arrange(.,desc(diff)) %>% 
   ggplot(mapping=aes(x=factor(celltype,levels=order$celltype),y=diff,fill=celltype)) +
    geom_bar(stat="identity") +
    theme_bw() +
    labs(title="Nuetrophil",x="celltypes") +
    scale_fill_manual(values=dittoColors()) +
    coord_flip()
ggsave(p,filename="figure7_diff_interactions.pdf")    



#评估树突状细胞抗原呈递能力
genes <- c("HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DRB1","HLA-DRB5","HLA-DRA","HLA-A","HLA-B","HLA-C") #人的
convertMouseGeneList1(genes)
convert_human_to_mouse_symbols(genes)
#小鼠
Dendritic_seurat <- subset(sce_merge_res1$sce,celltype=="Dendritic_cell")
genes <- c("H2-D1","H2-K1","H2-L","H-2Q","H2","H2-T","H2-Eb1","H2-Ea","H2-Bl")
a <- FetchData(Dendritic_seurat,vars = genes,slot = "data")  #"H2-D1,H2-K1,H2-Eb1"
a <- a %>% tibble::rownames_to_column("cell") 
a$sample <- stringr::str_split(a$cell,"_",simplify = T)[,1]
head(a)
a %>% 
   dplyr::group_by(sample) %>% 
   dplyr::summarise(`H2-D1`= mean(`H2-D1`),`H2-K1`=mean(`H2-K1`),`H2-Eb1`=mean(`H2-Eb1`)) %>% 
   tibble::column_to_rownames("sample") -> heatmap.data
pheatmap::pheatmap(t(heatmap.data[c("Vehicle","Sofa","R848","SR3"),]), cluster_rows = F,cluster_cols = F,filename = "figure7_dc_pheatmap.pdf")


#1. 小鼠的数据：四组中CD8+T cytotoxic评分的比较；
library(nichenetr)
sce_merge_res1 <- readRDS("./Results/data/mouse_sce_merge_res1.rds")
CD8_seurat <- subset(sce_merge_res1$sce,celltype=="CD8+T_cell")
eset <- CD8_seurat@assays$RNA@data
T_signatures <- fread("./code/tianjin_Supple/singlecell_mouse/T_marker.csv",header = T)
Cytotoxic_signature <- T_signatures %>% filter(.,term=="Cytotoxic signature") %>% dplyr::select(.,gene) %>% as.data.frame()
Cytotoxic_signature <- Cytotoxic_signature$gene
Exhaustion_signature <- T_signatures %>% filter(.,term=="Exhaustion signature") %>% dplyr::select(.,gene) %>% as.data.frame()
Exhaustion_signature <- Exhaustion_signature$gene
signature <- list(Exhaustion=as.character(na.omit(convert_human_to_mouse_symbols(Exhaustion_signature))),Cytotoxic=as.character(na.omit(convert_human_to_mouse_symbols(Cytotoxic_signature))))
signatures <- c(Cytotoxic_signature,T_signatures,list(TLS_Nature=IOBR::signature_tme$TLS_Nature),Interferon_signature,list(CYT=c("PRF1","GZMA")))
CD8_T_ssgsea <- calculate_sig_score_ssgsea(pdata = NULL,
                                           eset=eset,
                                           signature=signature,
                                           mini_gene_count=2,
                                           column_of_sample="ID",
                                           adjust_eset = FALSE)
CD8_T_ssgsea <- CD8_T_ssgsea %>% mutate(sample=stringr::str_split(ID,"_",simplify = T)[,1]) 
CD8_T_ssgsea %<>% reshape2::melt(.,id.vars=c("ID","Index","sample"))
CD8_T_ssgsea %>% 
  dplyr::group_by(sample,variable) %>% 
  dplyr::summarise(mean=mean(value)) %>% 
  right_join(.,CD8_T_ssgsea,by=c("sample","variable")) -> plot.data
split(plot.data,plot.data$variable) -> plot.data.list
p1 <- ggplot(data=plot.data.list[[1]],mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=value,fill=mean)) +
  geom_violin(scale = "width", adjust = 2, trim = TRUE,alpha=0.7) +
  ggsignif::geom_signif(comparisons = GetComb(c("Vehicle","Sofa","R848","SR3")),test = t.test,step_increase = 0.1) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), name = "Exhaustion signatures score") + 
  labs(x="sample",y="ssgsea score",title="Exhaustion signatures score compared between 4 groups within CD8+ T")+
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 
p2 <- ggplot(data=plot.data.list[[2]],mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=value,fill=mean)) +
  geom_violin(scale = "width", adjust = 2, trim = TRUE,alpha=0.7)+
  ggsignif::geom_signif(comparisons = GetComb(c("Vehicle","Sofa","R848","SR3")),test = t.test,step_increase = 0.1) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), name = "Cytotoxic signatures score") + 
  labs(x="sample",y="ssgsea score",title="Cytotoxic signatures score compared between 4 groups within CD8+ T")+
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 
p <- p1 + p2
ggsave(p,filename="original_Cytotoxic_Tcell_Exhaustion.pdf",width=12,height=8)
p <- CD8_T_ssgsea %>% 
  dplyr::group_by(sample,variable) %>% 
  dplyr::summarise(mean=mean(value)) %>% 
  right_join(.,CD8_T_ssgsea,by=c("sample","variable")) %>% 
  ggplot(mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=value,fill=mean)) +
  geom_violin(scale = "width", adjust = 2, trim = TRUE)+
  facet_wrap(vars(variable),scales = "free_y") + 
  ggsignif::geom_signif(comparisons = GetComb(c("Vehicle","Sofa","R848","SR3")),test = t.test,step_increase = 0.1) +
  scale_fill_gradient(low = "#F2E80D", high = "#E65224") +
  labs(x="sample",y="ssgsea score",title="signatures score compared between 4 groups within CD8+ T")+
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 
ggsave(p,filename = "figure7_Tcell_Exhaustion.pdf")
#热图展示，并加显著性P值
heat.mat <- CD8_T_ssgsea %>% 
  dplyr::group_by(sample,variable) %>% 
  dplyr::summarise(mean=mean(value)) %>% 
  reshape2::dcast(.,sample~variable)
heat.mat %<>% tibble::column_to_rownames("sample") %>% t() 
heat.mat <- heat.mat[,c("Vehicle","Sofa","R848","SR3")]
pheatmap::pheatmap(heat.mat,cluster_rows = F,cluster_cols = F, 
                   filename = "CD8+T_Exhaustion_Cytotoxic1.pdf",width=7,height=4)
cols <- colorRampPalette(colors = c("blue","white","red"))(100)
plotHeat(heat.mat,cluster_rows = F,cluster_cols = F,color =cols, width=7,height=4)
plotHeat(heat.mat,cluster_rows = F,cluster_cols = F,color =cols, filename = "CD8+T_Exhaustion_Cytotoxic2.pdf",width=7,height=4)

#换一套signatire集合重新计算
library(IOBR)
Cytotoxic_signature <- IOBR::signature_tme[grep("Cytotoxic|cytotoxic|exhaustion|Exhaustion|Exhausted",names(IOBR::signature_tme))][c(1:3,5)]
Cytotoxic_signature$Cytotoxic_cells_Danaher_et_al <-  as.character(na.omit(convert_human_to_mouse_symbols(Cytotoxic_signature$Cytotoxic_cells_Danaher_et_al)))
Cytotoxic_signature$Cytotoxic_cells_Bindea_et_al <-  as.character(na.omit(convert_human_to_mouse_symbols(Cytotoxic_signature$Cytotoxic_cells_Bindea_et_al)))
Cytotoxic_signature$T_cell_exhaustion_Peng_et_al <-  as.character(na.omit(convert_human_to_mouse_symbols(Cytotoxic_signature$T_cell_exhaustion_Peng_et_al)))
Cytotoxic_signature$Exhausted_CD8_Danaher_et_al <- as.character(na.omit(convert_human_to_mouse_symbols(Cytotoxic_signature$Exhausted_CD8_Danaher_et_al)))
Cytotoxic_signature$Exhausted_CD8_Zemin_et_al <- as.character(na.omit(convert_human_to_mouse_symbols(c("CTLA4","PDCD1","CXCL13","ENTPD1","LAG3","LAYN","TIGIT","BATF","HAVCR2","TNFRSF9","GZMB","CCL5","BAG3"))))
#Cytotoxic_signature$Exhausted_CD8_Zemin_et_al1 <- as.character(na.omit(convert_human_to_mouse_symbols(c("ENTPD1", "LAYN", "ITGAE", "BATF"))))
CD8_T_ssgsea1 <- calculate_sig_score_ssgsea(pdata = NULL,
                                           eset=eset,
                                           signature=Cytotoxic_signature,
                                           mini_gene_count=1,
                                           column_of_sample="ID",
                                           adjust_eset = FALSE)
CD8_T_ssgsea1 <- CD8_T_ssgsea1 %>% mutate(sample=stringr::str_split(ID,"_",simplify = T)[,1]) 
CD8_T_ssgsea1 %<>% reshape2::melt(.,id.vars=c("ID","Index","sample"))
CD8_T_ssgsea1 %>% 
  dplyr::group_by(sample,variable) %>% 
  dplyr::summarise(mean=mean(value)) %>% 
  right_join(.,CD8_T_ssgsea1,by=c("sample","variable")) -> plot.data
plot.data.list <- split(plot.data,plot.data$variable)

p1 <- ggplot(data=plot.data.list[[1]],mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=value,fill=mean)) +
  geom_violin(scale = "width", adjust = 2, trim = TRUE)+
  #facet_wrap(vars(variable),scales = "free_y") + 
  ggsignif::geom_signif(comparisons = GetComb(c("Vehicle","Sofa","R848","SR3")),test = t.test,step_increase = 0.1) +
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "white", "red"))(100), name = "Cytotoxic_CD8_Danaher_et_al") + 
  #ggrepel::geom_text_repel(data=annotation.data,mapping=aes(x=sample,y=max_value,label=mean)) +
  labs(x="sample",y="ssgsea score",title="Cytotoxic signatures score compared between 4 groups within CD8+ T")+
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 
p2 <- ggplot(data=plot.data.list[[2]],mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=value,fill=mean)) +
  geom_violin(scale = "width", adjust = 2, trim = TRUE)+
  #facet_wrap(vars(variable),scales = "free_y") + 
  ggsignif::geom_signif(comparisons = GetComb(c("Vehicle","Sofa","R848","SR3")),test = t.test,step_increase = 0.1) +
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "white", "red"))(100), name = "Exhausted_CD8_Danaher_et_al") + 
  #ggrepel::geom_text_repel(data=annotation.data,mapping=aes(x=sample,y=max_value,label=mean)) +
  labs(x="sample",y="ssgsea score",title="Exhausted signatures score compared between 4 groups within CD8+ T")+
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 
p3 <- ggplot(data=plot.data.list[[3]],mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=value,fill=mean)) +
  geom_violin(scale = "width", adjust = 2, trim = TRUE)+
  #facet_wrap(vars(variable),scales = "free_y") + 
  ggsignif::geom_signif(comparisons = GetComb(c("Vehicle","Sofa","R848","SR3")),test = t.test,step_increase = 0.1) +
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "white", "red"))(100), name = "Cytotoxic_CD8_Bindea_et_al") + 
  #ggrepel::geom_text_repel(data=annotation.data,mapping=aes(x=sample,y=max_value,label=mean)) +
  labs(x="sample",y="ssgsea score",title="Cytotoxic signatures score compared between 4 groups within CD8+ T")+
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 
p4 <- ggplot(data=plot.data.list[[4]],mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=value,fill=mean)) +
  geom_violin(scale = "width", adjust = 2, trim = TRUE) +
  #facet_wrap(vars(variable),scales = "free_y") + 
  ggsignif::geom_signif(comparisons = GetComb(c("Vehicle","Sofa","R848","SR3")),test = t.test,step_increase = 0.1) +
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "white", "red"))(100), name = "T_cell_exhaustion_Peng_et_al") + 
  #ggrepel::geom_text_repel(data=annotation.data,mapping=aes(x=sample,y=max_value,label=mean)) +
  labs(x="sample",y="ssgsea score",title="Exhausted signatures score compared between 4 groups within CD8+ T")+
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 
p5 <- ggplot(data=plot.data.list[[5]],mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=value,fill=mean)) +
          geom_violin(scale = "width", adjust = 2, trim = TRUE) +
          #facet_wrap(vars(variable),scales = "free_y") + 
          ggsignif::geom_signif(comparisons = GetComb(c("Vehicle","Sofa","R848","SR3")),test = t.test,step_increase = 0.1) +
          scale_fill_gradientn(colours = colorRampPalette(c("blue", "white", "red"))(100), name = "Exhausted_CD8_Zemin_et_al") + 
          #ggrepel::geom_text_repel(data=annotation.data,mapping=aes(x=sample,y=max_value,label=mean)) +
          labs(x="sample",y="ssgsea score",title="Exhausted signatures score compared between 4 groups within CD8+ T")+
          ggpubr::theme_pubr() + 
          theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 
p <- (p1 | p2 | p3) / (p4 | p5)
ggsave(p,filename = "CD8+T_Exhaustion_Cytotoxic_violinplot1.pdf",width = 14,height = 14)


annotation.data <- CD8_T_ssgsea1 %>% 
                      dplyr::group_by(sample,variable) %>% 
                      dplyr::summarise(mean=mean(value)) %>% 
                      right_join(.,CD8_T_ssgsea1,by=c("sample","variable")) %>% 
                      dplyr::group_by(sample,variable,mean) %>% 
                      dplyr::summarise(max_value=max(value)+0.1)
p <- CD8_T_ssgsea1 %>% 
  dplyr::group_by(sample,variable) %>% 
  dplyr::summarise(mean=mean(value)) %>% 
  right_join(.,CD8_T_ssgsea1,by=c("sample","variable")) %>% 
  ggplot(mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=value,fill=mean)) +
  geom_violin(scale = "width", adjust = 2, trim = TRUE)+
  facet_wrap(vars(variable),scales = "free_y") + 
  ggsignif::geom_signif(comparisons = GetComb(c("Vehicle","Sofa","R848","SR3")),test = t.test,step_increase = 0.1) +
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "white", "red"))(100), name = "z值") + 
  ggrepel::geom_text_repel(data=annotation.data,mapping=aes(x=sample,y=max_value,label=mean)) +
  labs(x="sample",y="ssgsea score",title="signatures score compared between 4 groups within CD8+ T")+
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 
ggsave(p,filename = "CD8+T_Exhaustion_Cytotoxic_violinplot.pdf",width = 10,height = 10)



#四组间neutrophil与DC cellphonedb ligand-receptor互作的结果展示；
#cellphonedb画图(中性粒细胞作为sender细胞)
Vehicle_cellphone_pvalues <- fread("/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/Vehicle/pvalues.txt")
Sofa_cellphone_pvalues <- fread("/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/Sofa/pvalues.txt")
R848_cellphone_pvalues <- fread("/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/R848/pvalues.txt")
SR3_cellphone_pvalues <- fread("/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/SR3/pvalues.txt")
Vehicle_cellphone_means <- fread("/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/Vehicle/significant_means.txt")
Sofa_cellphone_means <- fread("/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/Sofa/significant_means.txt")
R848_cellphone_means <- fread("/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/R848/significant_means.txt")
SR3_cellphone_means <- fread("/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output/SR3/significant_means.txt")
makeplotdata1 <- function(x,y,sample){
  Vehicle_cellphone_pvalues <- x
  Vehicle_cellphone_means <- y
  Vehicle_cellphone_pvalues <- as.data.frame(Vehicle_cellphone_pvalues)
  pos <- grep("Neutrophil_cell\\|",colnames(Vehicle_cellphone_pvalues))
  Vehicle_cellphone <- Vehicle_cellphone_pvalues[ ,colnames(Vehicle_cellphone_pvalues)[pos]]
  rownames(Vehicle_cellphone) <- Vehicle_cellphone_pvalues$interacting_pair
  Vehicle_cellphone1 <- Vehicle_cellphone %>% 
    tibble::rownames_to_column("interaction_pairs") %>% 
    reshape2::melt(.,id.vars=c("interaction_pairs")) %>% 
    data.table::setnames(.,old=c("interaction_pairs","variable","value"),new=c("ligand_receptor","cell_cell","p.value")) %>% 
    filter(.,p.value <= 0.05)
  Vehicle_cellphone_means <- as.data.frame(Vehicle_cellphone_means)
  pos <- grep("Neutrophil_cell\\|",colnames(Vehicle_cellphone_means))
  Vehicle_cellphone_means1 <- Vehicle_cellphone_means[ ,colnames(Vehicle_cellphone_means)[pos]]
  rownames(Vehicle_cellphone_means1) <- Vehicle_cellphone_means$interacting_pair
  Vehicle_cellphone_means1 <- Vehicle_cellphone_means1 %>% 
    tibble::rownames_to_column("interaction_pairs") %>% 
    reshape2::melt(.,id.vars=c("interaction_pairs")) %>% 
    data.table::setnames(.,old=c("interaction_pairs","variable","value"),new=c("ligand_receptor","cell_cell","mean.value")) %>% 
    filter(.,mean.value > 0)
  Vehicle_cellphone1 %>% inner_join(.,Vehicle_cellphone_means1,by=c("ligand_receptor","cell_cell")) %>% mutate(sample=rep(sample,length(ligand_receptor))) -> plot.data
  return(plot.data)
}
a1 <- makeplotdata1(Vehicle_cellphone_pvalues,Vehicle_cellphone_means,"Vehicle")
a2 <- makeplotdata1(Sofa_cellphone_pvalues,Sofa_cellphone_means,"Sofa")
a3 <- makeplotdata1(R848_cellphone_pvalues,R848_cellphone_means,"R848")
a4 <- makeplotdata1(SR3_cellphone_pvalues,SR3_cellphone_means,"SR3")
a <- rbind.data.frame(a1,a2,a3,a4)
#基因转换成小鼠
p <- a %>% filter(.,cell_cell=="Neutrophil_cell|Dendritic_cell",ligand_receptor != "C5AR1_RPS19") %>% 
        ggplot(mapping = aes(x=factor(sample,levels = c("Vehicle","Sofa","R848","SR3")),y=ligand_receptor,size=-log10(p.value+0.00000001),color=mean.value)) +
          geom_point() + 
          scale_color_gradient(low = "#F2E80D", high = "#E65224") +
          theme_bw() +
          labs(x="sample",title = "Neutrophil_cell_Dendritic_cell",size="-log10(p.value)")
ggsave(p,filename="figure7_Neutrophil_cell_Dendritic_cell.pdf")


