library(Seurat)
library(ggplot2)
library(data.table)
library(tidyverse)
library(dplyr)
library(Matrix)
library(circlize)
library(readr)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(monocle)
library(gplots)
library(RColorBrewer)
library(pheatmap)

#----------------------数据导入-------------------
#
#90cell yan data
a = as.matrix(fread(".\\data\\yan_data.csv")) #有时要去除最后的无效基因
b = fread(".\\data\\yan_genedata.csv")
genename = a[,1]
a = a[,-1]
rownames(a) = genename
c = fread(".\\yan\\yan_celldata.csv",header = T)
label = c$cell_type1
sce_meta = data.frame(row.names = colnames(a),celltype = label)
yan = CreateSeuratObject(counts = a, meta.data = sce_meta,project = "Yan")
#yan = CreateSeuratObject(counts = a, meta.data = sce_meta,project = "Yan",min.cells = 3, min.features = 100)
#-----------------------质控-----------------------
#寻找线粒体基因

yan[["percent.mt"]] = PercentageFeatureSet(yan,pattern = "^MT-") #求线粒体基因百分比
table(yan[["percent.mt"]])
#小提琴可视化，展示三个图，每个图可根据group.by=''分群
VlnPlot(object = yan, features = c("nFeature_RNA","nCount_RNA","percent.mt"))
#过滤，去除表达gene数小于500和大于6000的细胞，去除线粒体RNA比例高于20%的细胞
#yan = subset(yan,subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt <20 )
#plot，图上方显示相关系数
#plot1 = FeatureScatter(yan, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(yan, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2

plot22 = VlnPlot(object = yan,
        features = c("C9orf152","RPS11","ELMO2","CREB3L1","PNMA1","MMP2","TMEM216","LOC653712","C10orf90","ZHX3"),
        ncol = 5,
        #group.by = 'singleR',
        pt.size = 0.4,
        ) 

#---------------------------------scGCT插补--------------------------------
a = as.matrix(fread(".\\yan_imputed.csv"))
a = a[,-1]
rownames(a) = genename
#----------------------标准化、高变基因-----------------------
yan <- NormalizeData(yan,scale.factor =10000)
yan <- FindVariableFeatures(yan,nfeatures = 1500)
yan <- ScaleData(yan)
#输出特征方差图
top10 = head(VariableFeatures(yan),10)
plot4 = VariableFeaturePlot(yan)
plot5 = LabelPoints(plot = plot4,points = top10,repel = TRUE)
plot5

#-------------------------降维-------------------------------
# 运行PCA分析并计算投影
yan <- RunPCA(yan)
#查看前五个PC-主成分
print(yan[["pca"]],dim = 1:5,nfeatures = 5)
#可视化
VizDimLoadings(object = yan,dim = 1:4,reduction = "pca",nfeatures = 20) #两个主成分，20个基因
DimPlot(yan,reduction = "pca",group.by = "celltype")
DimHeatmap(yan,dims = 1:10,cells = 90, balanced = TRUE,nfeatures = 30,ncol = 5) #两个主成分热图，一行两个图
yan = JackStraw(object = yan, num.replicate = 100)
yan = ScoreJackStraw(object = yan, dim = 1:20)
JackStrawPlot(object = yan,dims = 1:20,reduction = "pca")
ElbowPlot(yan,reduction = "pca") #根据ElbowPlot确定PC数量
#-------------------------聚类-------------------------------
yan <- FindNeighbors(yan, dims = 1:20)
yan <- FindClusters(yan,)
yan <-RunUMAP(yan,dims = 1:20)
yan <-RunTSNE(yan,perplexity=10)
#-----------------------可视化-------------------------------
DimPlot(yan, reduction = "umap", group.by = "seurat_clusters") # Visualize custom clusters
DimPlot(yan, reduction = "umap", group.by = "celltype",label = T,pt.size = 1)
DimPlot(yan, reduction = "tsne", group.by = "celltype",label = T) # Visualize custom clusters
DimPlot(yan, reduction = "pca", group.by = "celltype") # Visualize custom clusters

#--------------------差异基因表达---------------------------
Idents(yan) <- label
yan.markers <- FindAllMarkers(object = yan,logfc.threshold=0.05)
#取每种细胞类型按log2FC值从大到小排序的前五个差异基因，用来绘制热图
top5 <- yan.markers%>%group_by(cluster)%>%top_n(n = 5, wt = abs(avg_log2FC))
mat<-GetAssayData(yan, slot = "counts")
gene_features <- top5
cluster_info <- sort(yan@active.ident)
mat <- as.matrix(mat[top5$gene, names(cluster_info)])
##DoHeatmap(object = data, features = top5$gene) + NoLegend()

#pheatmap<- pheatmap(exp,cluster_cols=F,annotation_col=phe,show_colnames=F,scale='row',show_rownames=F)
#计算得到的差异基因在调整后的P值<1e-5且log2FC>2的基因个数
FC <- yan.markers$avg_log2FC[yan.markers$p_val_adj < 1e-5]
names(FC) <- rownames(yan.markers[yan.markers$p_val_adj < 1e-5,])
DEs_uncorr <- list()
#DEuncorr[] <- rownames(pbmc.markers[pbmc.markers$p_val_adj < 1e-5 & pbmc.markers$avg_log2FC > 2,])
DEs_uncorr <- yan.markers$gene[yan.markers$p_val_adj < 1e-5 & yan.markers$avg_log2FC > 0]
#657
library(ComplexHeatmap)
Heatmap(mat,                          #数据矩阵
        width = unit(8, "cm"),        #图片宽度
        height = unit(7.5, "cm"),     #图片高度
        name = "Yan",
        border = TRUE,
        #rect_gp = gpar(col = "blue", lwd = 0.01),  
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_gp = gpar(family = 'Broman',fontsize = 7.5, fontface = "bold"),
        column_split = cluster_info,
        top_annotation = HeatmapAnnotation(
          Cluster = cluster_info,
          annotation_legend_param =(cluster=list(anno_block(gp = gpar(fill = col),
                                                            labels = levels(cluster_info))
          ))),
        column_title = 'Cell Type',
        column_title_gp = gpar(family = 'Broman', fontsize = 15),
        row_title = "Gene Makers",
        row_title_gp = gpar(family = 'Broman', fontsize = 15),)

#火山图
EnhancedVolcano(top5,
                lab=top5$gene,
                x = 'avg_log2FC',
                y ='p_val')
#展示感兴趣基因
FeaturePlot(yan, features = c('BOD1','RGS2','TUBB8','MBD3L2','LDHA','KRT18'),
            ncol = 5)

#------------------------细胞轨迹推断monocle---------------------------
expr_matrix<-as(as.matrix(yan@assays[["RNA"]]@counts),'sparseMatrix')
pd<-yan@meta.data
#pd$celltype<-bh@active.ident
fd<-data.frame(gene_short_name=row.names(yan),row.names = row.names(yan))

pd<-new('AnnotatedDataFrame',data=pd)
fd<-new('AnnotatedDataFrame',data=fd)
yan_monocle <- newCellDataSet(expr_matrix,featureData = fd,phenoData = pd,
                             lowerDetectionLimit = 0.5,
                             expressionFamily = negbinomial.size())
# 计算主要图形
yan_monocle <- estimateSizeFactors(yan_monocle)
yan_monocle <- estimateDispersions(yan_monocle)
#选择聚类的差异表达基因

library(gplots)
library(RColorBrewer)
library(pheatmap)

#选择dpf方法的差异表达基因
expressed_genes <- row.names(subset(fData(yan_monocle))) #过滤gene，num_cells_expressed>=10
diff_test_res <- differentialGeneTest(yan_monocle[expressed_genes,],
                                      fullModelFormulaStr = "~celltype")
##差异表达基因作为轨迹构建的基因，差异基因的选择标准是qval<0.01，decreasing=F表示数值按增加排序
deg<-subset(diff_test_res,qval<0.0001)
deg<-deg[order(deg$qval,decreasing=F),]
ordering_genes <- row.names (deg)

##将得到的基因列表嵌入对象
yan_monocle <- setOrderingFilter(yan_monocle, ordering_genes)#基因被储存在bh_monocle@featureData@data[["use_for_ordering"]]
plot_ordering_genes(yan_monocle)
#黑色的点用来构建轨迹的差异基因，灰色的点是背景基因，红色的线是根据基因表达大小和离散度分布的趋势（找到的基因属于离散度较高的基因）

#降维

bh_monocle <- reduceDimension(bh_monocle,max_components = 2,
                              method = 'DDRTree')

#pbmc_small_monocle <- preprocessCDS(pbmc_small_monocle, num_dim = 10)
#pbmc_small_monocle <- reduceDimension(pbmc_small_monocle,max_components = 2,
#                            method = 'DDRTree')
#pbmc_small_monocle <- clusterCells(pbmc_small_monocle)
#pbmc_small_monocle <- learnGraph(pbmc_small_monocle)

# 对单元进行排序
bh_monocle <- orderCells(bh_monocle)
plot_cell_trajectory(bh_monocle, color_by = "seurat_clusters")
#plot_cell_trajectory(bh_monocle, color_by = "singleR")
plot_cell_trajectory(bh_monocle, color_by = "Pseudotime")
plot_cell_trajectory(bh_monocle, color_by = "celltype")
plot_cell_trajectory(bh_monocle, color_by = "celltype") +
  facet_wrap(~State, nrow = 1)




