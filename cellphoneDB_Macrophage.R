#' 观察巨噬细胞（M1,M2）cellphoneDB
#' @author shijian
library(dplyr)
library(ggsignif)
library(data.table)
library(Seurat)
library(ggplot2)
library(tidytext)
library(forcats)
library(stringr)
library(nichenetr)
sce_merge_res1 <- readRDS("./Results/data/mouse_sce_merge_res1.rds")
Marcrophage_sce <- readRDS(file="./Results/mouse/Marcrophage.rds")
Marcrophage_sce$celltype2 <- as.character(plyr::mapvalues(Marcrophage_sce$seurat_clusters,from=c(0:1),to=c("M2","M1")))
M1_Marcrophage <- subset(Marcrophage_sce,celltype2=="M1")
M2_Marcrophage <- subset(Marcrophage_sce,celltype2=="M2")
# meta_data <- sce_merge_res1$sce@meta.data %>% tibble::rownames_to_column("cell")  %>% 
#               left_join(.,tibble::rownames_to_column(as.data.frame(Marcrophage_sce@meta.data[,"celltype2"],row.names=rownames(Marcrophage_sce@meta.data)),var="cell"),by="cell") %>% 
#               dplyr::rename(.,celltype2='Marcrophage_sce@meta.data[, "celltype2"]')
unique(sce_merge_res1$sce$celltype)
sce_merge_res1$sce$celltype2 <- ifelse(sce_merge_res1$sce$celltype=="Marcrophage",
                                       ifelse(rownames(sce_merge_res1$sce@meta.data) %in% colnames(Marcrophage_sce),ifelse(rownames(sce_merge_res1$sce@meta.data) %in% colnames(M1_Marcrophage),"M1","M2"),"remove"),
                                       sce_merge_res1$sce$celltype)
merge_sce <- subset(sce_merge_res1$sce,celltype2 != "remove")
p1 <- dittoDimPlot(merge_sce, reduction.use  = "umap",var="celltype2",size = 1,do.label = TRUE, labels.highlight = F,do.ellipse = F,legend.size = 9,shape.legend.size=9,labels.size=7,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 16))
p1
##########注意cellphonedb是一个样本一个样本的跑##########
#制作每个样本的细胞注释类型文件
cellAnno_outDir <- "code/tianjin_Supple/singlecell_mouse/cellAnno_outDir1"
if(!file.exists(cellAnno_outDir)){dir.create(cellAnno_outDir,recursive = T)}
metadata <- merge_sce@meta.data %>% tibble::rownames_to_column("cell") %>% dplyr::select(.,orig.ident,cell,celltype2)
samples_cellAnno <- split(metadata,f=metadata$orig.ident)
lapply(samples_cellAnno,function(x){
  name <- x %>% dplyr::select(.,orig.ident) %>% distinct()
  temp <- x %>% dplyr::select(.,c(cell,celltype2))
  write.table(temp,file = paste0("./code/tianjin_Supple/singlecell_mouse/cellAnno_outDir1/cellphonedb_metadata_",name,".txt"),sep="\t",quote=F,row.names=F,col.names=T)
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
cellPhone_exprDir <- "code/tianjin_Supple/singlecell_mouse/cellPhone_exprDir1"
if(!file.exists(cellPhone_exprDir)){dir.create(cellPhone_exprDir,recursive = T)}
samples <- unique(merge_sce@meta.data$orig.ident)
for(i in 1:length(samples)){
  subset_seurat <- subset(merge_sce,orig.ident==samples[i])
  makeCellPhoneExpr(subset_seurat,samples[i],cellPhone_exprDir)
}
#运行cellphonedb程序
cd /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse
mkdir cellphonedb_output1
cd cellphonedb_output1
mkdir R848
cellphonedb method statistical_analysis --counts-data gene_name \
--output-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output1/R848 \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellAnno_outDir1/cellphonedb_metadata_R848.txt \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellPhone_exprDir1/cellphonedb_expr_R848.txt \
--threads 10
cd /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output1
mkdir Sofa
cellphonedb method statistical_analysis --counts-data gene_name \
--output-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output1/Sofa \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellAnno_outDir1/cellphonedb_metadata_Sofa.txt \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellPhone_exprDir1/cellphonedb_expr_Sofa.txt \
--threads 10
cd /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output1
mkdir Vehicle
cd Vehicle
cellphonedb method statistical_analysis --counts-data gene_name \
--output-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output1/Vehicle \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellAnno_outDir1/cellphonedb_metadata_Vehicle.txt \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellPhone_exprDir1/cellphonedb_expr_Vehicle.txt \
--threads 10
cd /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output1
mkdir SR3
cd SR3
cellphonedb method statistical_analysis --counts-data gene_name \
--output-path /Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output1/SR3 \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellAnno_outDir1/cellphonedb_metadata_SR3.txt \
/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellPhone_exprDir1/cellphonedb_expr_SR3.txt \
--threads 10
# output
# deconvoluted.txt：每对受体-配体提供了附加信息
# mean.txt：每对受体-配体相互作用的表达平均值，所有相互作用伙伴的平均值：平均值是指相应相互作用的细胞类型对中个体伙伴平均表达值的总平均值。
#           如果其中一个平均值为 0，则总平均值设置为 0。
# pvalues.txt: 每对受体-配体的P值，p 值是指在每个相互作用的细胞类型对中相互作用的配体-受体对的富集
# significant_means.txt: 每对受体-配体显著性结果的平均表达量值。所有交互伙伴的显著平均计算。如果 p.value < 0.05，则该值将是平均值。或者，该值设置为 0。
#重要的是要记住，相互作用不是对称的。换言之，当测试簇 X_Y 之间的配体/受体对 A_B 时，伙伴 A 的表达被认为在第一个簇 (X) 内，而伙伴 B 在第二个簇 (Y) 内的表达。
#因此，X_Y 和 Y_X 代表不同的比较，将具有不同的 p 值和均值。
Vehicle_cellphone_pvalues <- fread("/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output1/Vehicle/pvalues.txt")
Sofa_cellphone_pvalues <- fread("/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output1/Sofa/pvalues.txt")
R848_cellphone_pvalues <- fread("/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output1/R848/pvalues.txt")
SR3_cellphone_pvalues <- fread("/Users/biofly/project/shijian/code/tianjin_Supple/singlecell_mouse/cellphonedb_output1/SR3/pvalues.txt")
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
  d <- str_split(p$celltype_pair,"\\|",simplify = T)
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
a1 <- makeplotdata(Vehicle_cellphone_pvalues,"M1","Vehicle")
a2 <- makeplotdata(Sofa_cellphone_pvalues,"M1","Sofa")
a3 <- makeplotdata(R848_cellphone_pvalues,"M1","R848")
a4 <- makeplotdata(SR3_cellphone_pvalues,"M1","SR3")
a <- rbind.data.frame(a1,a2,a3,a4)
a
b1 <- makeplotdata(Vehicle_cellphone_pvalues,"M2","Vehicle")
b2 <- makeplotdata(Sofa_cellphone_pvalues,"M2","Sofa")
b3 <- makeplotdata(R848_cellphone_pvalues,"M2","R848")
b4 <- makeplotdata(SR3_cellphone_pvalues,"M2","SR3")
b <- rbind.data.frame(b1,b2,b3,b4)
b
p1 <- ggplot(data=a,mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=num,fill=sample)) + 
  geom_bar(stat="identity",position=position_dodge()) +
  facet_grid(celltype~direction,scales = "free_y") +
  scale_fill_manual(values = dittoColors()) + 
  labs(x="sample",y="interaction numbers",title = "M1") +
  theme_bw() +
  theme(axis.ticks.length=unit(0.5,'cm'),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size=20)) + 
  theme(strip.text = element_text(size = 8,face = "bold")) +
  coord_flip()
ggsave(p1,filename = "M1_cellphonedb.pdf",width=16,height=12)
print(p1)
p1 <- ggplot(data=b,mapping=aes(x=factor(sample,levels=c("Vehicle","Sofa","R848","SR3")),y=num,fill=sample)) + 
  geom_bar(stat="identity",position=position_dodge()) +
  facet_grid(celltype~direction,scales = "free_y") +
  scale_fill_manual(values = dittoColors()) + 
  labs(x="sample",y="interaction numbers",title = "M2") +
  theme_bw() +
  theme(axis.ticks.length=unit(0.5,'cm'),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size=20)) + 
  theme(strip.text = element_text(size = 10,face = "bold")) +
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
      ggplot(mapping=aes(x=factor(celltype,levels=order$celltype),y=diff,fill=celltype)) +
      geom_bar(stat="identity") +
      theme_bw() +
      labs(title="M1",x="celltype",y="different interaction pairs in SR3") +
      scale_fill_manual(values=dittoColors()) +
      theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
ggsave(p,filename = "figure6_different_interactions.pdf")
