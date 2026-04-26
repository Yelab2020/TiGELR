MyColorCode <- c( '#48a854', '#ffe119', '#526cc9', '#f58231','#d83d53',
                  '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4',
                  '#469990', '#dcbeff', '#e8a767', '#f1e7ae', '#cc5500',
                  '#aaffc3', '#808000', '#ffd8b1', '#6c51a3', '#d3d3d3')
colors <- c('#dc8e97','#e3d1db','#74a893','#ac9141','#5ac6e9','#e5c06e','#7587b1','#c7deef','#e97371','#e1a4c6',
            '#916ba6','#cb8f82','#7db3af','#d2e0ac','#f3f3b0','#a5d177','#e0bc58','#64abc0','#fab37f','#e98741',
            '#8fc0dc','#967568','#f2d3ca','#eebd85','#82c785','#edeaa4','#cdaa9f','#794976','#bcacd3','#889b5d',
            '#4e9592','#dbad5f','#64ae79','#64ae79','#ac5092','#f37d74','#4a6eb6','#bbe3ee','#a5d177','#5b9951',
            '#fbd9ae','#cf622b','#90cdc1','#f4eb9b','#5f9d58','#f5949f','#c1cf9a','#70afab','#a8c375','#8e95c9',
            '#cf622b','#efcfee','#5384c4','#cacae3','#e45878','#996da8','#3b99d4','#edf5e1')

colors_tumor_niche <- c( '#48a854', '#ffe119', '#f58231','#d83d53',
                         '#911eb4', '#f032e6', '#bfef45', '#fabed4',
                         '#469990', '#dcbeff', '#f1e7ae', '#cc5500',
                         '#aaffc3')

colors_MP <- c('#ffe119', 'red', '#82c785','#911eb4',
               '#d83d53', '#42d4f4', '#f032e6', '#bfef45', '#fabed4')

CellType_colors <- c( '#ffe119','#48a854','#7db3af','#bfef45','#dcbeff','#e97371','#911eb4','#f2d3ca',
                              '#42d4f4', '#f032e6','#e8a767', '#f1e7ae', '#cc5500')


library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)  
library(magrittr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(harmony)
library(tidydr)
library(readr)


cosmx <- readr::read_rds('merge_data_final.rds.gz') 

Epi <- subset(cosmx, subtype == 'EpithelialCells')
Epi_sub <- subset(Epi, tissue =='T')

set.seed(123)
VGENES=rownames(Epi_sub)
VGENES=setdiff(VGENES,VGENES[grep("^SystemControl",VGENES)])
VGENES = VGENES[!grepl("Neg|NegProbe", VGENES)]
Epi_sub <- Epi_sub[VGENES,]
Epi_sub <- NormalizeData(Epi_sub)
Epi_sub <- FindVariableFeatures(object = Epi_sub , mean.function = ExpMean, dispersion.function = LogVMR)
Epi_sub <- ScaleData(object = Epi_sub)
Epi_sub <- RunPCA(Epi_sub, npcs = 50)
ElbowPlot(Epi_sub)
WholeTissueHarmony <- Epi_sub %>% 
  RunHarmony("Run_Tissue_name", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(WholeTissueHarmony, 'harmony')
WholeTissueHarmony <- WholeTissueHarmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>%    # dims = 1:30
  FindNeighbors(reduction = "harmony", dims = 1:20)  # dims = 1:30
WholeTissueHarmony <- FindClusters(WholeTissueHarmony, resolution = 0.8, algorithm = 1) 


features = c("LYZ","TIMP1","AREG","CEACAM6","S100A4","S100P","CD44","MYC","CCND1")

WholeTissueHarmony <- AddModuleScore(WholeTissueHarmony,
                                     features = list(features),
                                     name = "Tumor_Score")

FeaturePlot(WholeTissueHarmony, raster=FALSE, features = 'Tumor_Score1',cols = c("#7BAFDE", '#ffff00', "#DC050C"))+
  ggtitle("Tumor Marker Score")+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(WholeTissueHarmony)%in%rownames(niche))
WholeTissueHarmony$niche <- niche[colnames(WholeTissueHarmony), 1]

library(scCustomize)
Plot_Density_Joint_Only(WholeTissueHarmony, features = features)
Stacked_VlnPlot(WholeTissueHarmony, group.by = 'niche',features = features, x_lab_rotate = TRUE)


DimPlot(WholeTissueHarmony, cols = colors, group.by = 'CellType', raster=FALSE)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))


DimPlot(WholeTissueHarmony, cols = colors, group.by = 'seurat_clusters', raster=FALSE)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))


DotPlot(WholeTissueHarmony, group.by = 'CellType',features = features) +
  RotatedAxis()+
  scale_color_gradientn(colors = c("#7BAFDE", '#ffff00', "#DC050C")) +
  scale_size_continuous(range = c(1, 10))

#### Spatial ####
MSI_T <- subset(cosmx, subset = sample == 'MSS_T6')

MSI_T$CellType <- ifelse(MSI_T$CellType %in% c('EpithelialCells', 'TumorCells'), as.character(MSI_T$CellType),'Others')
MSI_T$CellType <- factor(MSI_T$CellType, levels = c('EpithelialCells', 'TumorCells', 'Others'))
simplified_colors <- c('Others' = "grey",'TumorCells' = "red",'EpithelialCells' = "blue")
DefaultBoundary(MSI_T[["A_CRC.1"]]) <- "segmentation"
ImageDimPlot(MSI_T, group.by = 'CellType', size = 1, cols = simplified_colors, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cosmx$CellType <- factor(cosmx$CellType, levels =  c('EndothelialCells','EpithelialCells','Glia','Fibroblasts','MastCells','DC','Macrophages',
                                                                     'Neutrophils','BCells','PlasmaCells','TCells','TumorCells') )

niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(niche))

cosmx$niche <- niche[colnames(cosmx), 1]

df <- cosmx@meta.data[,c('CellType','niche')]
df$niche <- paste0("Niche", df$niche)
prop_table <- prop.table(table(df$niche, df$CellType), margin = 1)
prop_matrix <- as.data.frame.matrix(t(prop_table))

cell_types <- c('EndothelialCells','EpithelialCells','Glia','Fibroblasts','MastCells','DC','Macrophages',
                'Neutrophils','BCells','PlasmaCells','TCells','TumorCells')
cell_colors <- setNames(CellType_colors[1:length(cell_types)], cell_types)
row_ha <- rowAnnotation(
  CellType = factor(cell_types),
  col = list(CellType = cell_colors),
  show_annotation_name = FALSE,  
  show_legend = F          
)

library(circlize)
library(viridis)


Heatmap(
  prop_matrix,
  name = "Composition", 
  # col = colorRampPalette(c("#7BAFDE", "#FFFF00", "#DC050C"))(100),
  col = colorRamp2(seq(0, 1, length.out = 100), viridis(100)),
  row_names_side = "left",
  left_annotation = row_ha,
  show_row_names = T,       
  row_title = NULL,               
  # column_names_side = "top",
  column_names_rot = 90,
  cluster_rows = FALSE,
  cluster_columns = TRUE, 
  show_column_dend = TRUE,  
  rect_gp = gpar(col = "gray", lwd = 0.5),  

  heatmap_legend_param = list(
    at = c(0, 0.5, 1), 
    labels = c("0", "0.5", "1")
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (prop_matrix[i, j] > 0.3) {
      grid.text("*", x, y, gp = gpar(fontsize = 16)) 
    }
  }
)


#### niche-celltype barplot ####

ClusterFreq <- cosmx@meta.data[,c("CellType","niche")] %>% table %>%
  data.frame() %>% set_colnames(c("CellType","niche","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(CellType,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=CellType,value=Per,-niche)
ClusterPer$niche=factor(ClusterPer$niche,levels = names(table(cosmx$niche)))

ClusterPer$CellType <- factor(ClusterPer$CellType, 
                                      levels = c('EndothelialCells','EpithelialCells','Glia','Fibroblasts','MastCells','DC','Macrophages',
                                                 'Neutrophils','BCells','PlasmaCells','TCells','TumorCells'))

ggplot(data = ClusterPer, aes(x = niche, fill=CellType,y=Per)) + #log2(as.numeric(CopyNumber)))
  geom_bar(stat = "identity",width = 0.6)+
  scale_y_continuous(expand=c(0,0))+labs(y="Percentage(%)")+
  # coord_flip()+
  scale_fill_manual(values=CellType_colors)+
  theme(axis.text=element_text(color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black"),
        axis.title.x = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10),
        legend.position = "none")




MP <- read_csv('MP.csv')
gene_list_by_cluster <- as.list(MP)

go_results <- list()
for (cl in names(gene_list_by_cluster)) {
  genes <- gene_list_by_cluster[[cl]]
  gene_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_ids) >= 10) {  
    ego <- enrichGO(gene          = gene_ids$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "ENTREZID",
                    ont           = "ALL",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.01,
                    readable      = TRUE)
    if (!is.null(ego) && nrow(ego@result) > 0) {
      ego@result$Cluster <- cl
      go_results[[cl]] <- ego
    } else {
      message(paste0("No GO enrichment found for cluster ", cl))
    }
  } else {
    message(paste0("Too few genes mapped for cluster ", cl))
  }
}

top_go_df <- do.call(rbind, lapply(go_results, function(res) {
  if (!is.null(res) && nrow(res@result) > 0) {
    res_df <- res@result %>%
      arrange(p.adjust) %>%
      slice_head(n = 10) %>%
      select(Description, Cluster = Cluster, p.adjust)
    return(res_df)
  }
}))
top_go_df$Cluster <- paste0('MP ',top_go_df$Cluster)

top_go_df$log10p <- -log10(top_go_df$p.adjust)

heatmap_df <- top_go_df %>%
  pivot_wider(names_from = Cluster, values_from = log10p, values_fill = 0) %>%
  select(-p.adjust) %>%
  group_by(Description) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)), .groups = "drop")

mat <- as.data.frame(heatmap_df)
rownames(mat) <- mat$Description
mat <- as.matrix(mat[,-1])
rownames(mat) <- stringr::str_wrap(rownames(mat), width = 50) 

desired_order <- paste0('MP ',c(11,14,9,10,4,6,13,2,1,7,12,16,8,17,3,5))
Heatmap(mat,
        name = "-log10(p.adjust)",  # legend title
        col = colorRampPalette(c("#7BAFDE", "#FFFF00", "#DC050C"))(100),
        column_names_rot = 0,
        column_order = desired_order, 
        cluster_columns = F,
        cluster_rows = TRUE,
        show_column_dend = F, 
        show_row_dend = F,
        rect_gp = gpar(col = "gray", lwd = 0.5)
)


test <- apply(as.data.frame(gene_list_by_cluster), 2, function(x){
  obj <- AddModuleScore(sc,features = list(x))
  obj@meta.data[,"Cluster1"]
})
test <- as.data.frame(test)

for (i in 1:17) {
  sc[[paste0("MP", i)]] <- test[, i]
}


Modules <- sc@meta.data[, c(11:27)]
max_idx <- max.col(Modules, ties.method = "first")
max_vals <- Modules[cbind(seq_len(nrow(Modules)), max_idx)]
second_vals <- apply(Modules, 1, function(x) sort(x, decreasing = TRUE)[2])
sc$MP <- ifelse(max_vals * 0.9 > second_vals,
                colnames(Modules)[max_idx],
                "Unassigned")

sc$MP[sc$MP == 'MP1'] <- 'MYC'
sc$MP[sc$MP == 'MP2'] <- 'EMT-III'
sc$MP[sc$MP == 'MP3'] <- 'PDAC-classical'
sc$MP[sc$MP == 'MP4'] <- 'Stress'
sc$MP[sc$MP == 'MP5'] <- 'PDAC-classical'
sc$MP[sc$MP == 'MP6'] <- 'EMT-II'
sc$MP[sc$MP == 'MP7'] <- 'Hypoxia'
sc$MP[sc$MP == 'MP8'] <- 'Unassigned'
sc$MP[sc$MP == 'MP9'] <- 'Cell Cycle'
sc$MP[sc$MP == 'MP10'] <- 'MYC'
sc$MP[sc$MP == 'MP11'] <- 'Cell Cycle'
sc$MP[sc$MP == 'MP12'] <- 'Unassigned'
sc$MP[sc$MP == 'MP13'] <- 'PDAC-classical'
sc$MP[sc$MP == 'MP14'] <- 'Interferon'
sc$MP[sc$MP == 'MP15'] <- 'Translation init.' 
sc$MP[sc$MP == 'MP16'] <- 'Unassigned'
sc$MP[sc$MP == 'MP17'] <- 'Unassigned'

sc$MP[sc$MP == 'Translation init.'] <- 'Unassigned' 

set.seed(42)
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(object = sc , mean.function = ExpMean, dispersion.function = LogVMR)
sc <- ScaleData(object =  sc)
sc <- RunPCA(object =  sc,  npcs = 30, verbose = FALSE)   # npcs = 50
WholeTissueHarmony <- sc %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(WholeTissueHarmony, 'harmony')
WholeTissueHarmony <- WholeTissueHarmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>%    # dims = 1:30
  FindNeighbors(reduction = "harmony", dims = 1:20)
DimPlot(WholeTissueHarmony, group.by = 'MP', cols = colors_MP)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))
DimPlot(WholeTissueHarmony, group.by = 'Module', cols = MyColorCode)

WholeTissueHarmony <- readr::read_rds('ref_tumor_mp.rds.gz')
WholeTissueHarmony <- subset(WholeTissueHarmony, MP != 'Unassigned')


sc <- readr::read_rds('ref_tumor_mp.rds.gz')
Idents(sc) <- 'MP'
DimPlot(sc, group.by = 'MP', cols = colors_MP)
markers <- FindAllMarkers(sc,
                          only.pos = TRUE, 
                          min.pct = 0.2,
                          logfc.threshold = 0.25)  # logfc.threshold = 0.25

genes <- markers$gene[markers$cluster == 'EMT-II' & markers$p_val_adj < 0.0001 & markers$avg_log2FC > 1.5]
genes[genes %in% df$`MP13 EMT-II`]

genes <- c("LDHB","HSPE1","FABP5",
           "CEACAM7","PIGR","CA2","PLAC8","TFF3",
           "EGR1","FOSB","ATF3","BTG2", 
           "LAMC2","PLAUR","LAMB3","PMEPA1","FLNA",
           "S100P","CRIP1","MMP7",
           "NDRG1","ERO1A","VEGFA","PLOD2","P4HA1",
           "TOP2A","CENPF","MKI67","ASPM","TPX2",
           "HLA-DRA","HLA-DRB1","CD74", "SAA1","HLA-DPA1"
)


desired_order <- c("MYC", "PDAC-classical", "Stress","EMT-II", "EMT-III","Hypoxia","Cell Cycle","Interferon")
Idents(sc) <- factor(sc$MP, levels = desired_order)
sc$MP <- factor(sc$MP, levels = desired_order)

library(Scillus)
colors_MP_order <- c('#42d4f4','#f032e6', '#bfef45', 'red', '#82c785','#911eb4','#ffe119','#d83d53')
plot_heatmap(dataset = sc,
             markers = genes,
             sort_var = c("MP"),
             anno_var = c("MP"),
             hm_limit = c(0, 1, 7),
             hm_colors = c("#4575b4", "white", "#d73027"),
             anno_colors = list(colors_MP_order
             ))




ref.obj="ref_tumor_mp.rds.gz"
sc.celltype="MP"

cosmx <- readr::read_rds('merge_data_final.rds.gz')
mal <- subset(cosmx, subset = CellType == "TumorCells")


library(corrplot)  
df <- mal@meta.data[,c('Mean.PanCK','Mean.G','Mean.Membrane','Mean.CD45','Mean.DAPI')]
cor_matrix <- cor(df, method = "pearson") 
corrplot(cor_matrix, method = "color", type = "full")
mal=CosMxInSituType(mal,ref.obj,celltype=sc.celltype,immunofluorescence=c("Membrane","PanCK","DAPI"),species="Human") 

mal$MP <- mal$CellType

set.seed(123)
VGENES=rownames(mal)
VGENES=setdiff(VGENES,VGENES[grep("^SystemControl",VGENES)])
VGENES = VGENES[!grepl("Neg|NegProbe", VGENES)]
mal <- mal[VGENES,]
mal <- SCTransform(mal, assay = "Nanostring", clip.range = c(-10, 10))
mal <- RunPCA(mal, npcs = 30)
WholeTissueHarmony <- mal %>%
  RunHarmony("Run_Tissue_name", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(WholeTissueHarmony, 'harmony')
WholeTissueHarmony <- WholeTissueHarmony %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%    # dims = 1:30
  FindNeighbors(reduction = "harmony", dims = 1:20)  # dims = 1:30
readr::write_rds(WholeTissueHarmony,'malignant.rds.gz',compress = 'gz')

WholeTissueHarmony <- readr::read_rds('malignant.rds.gz')
colnames(WholeTissueHarmony@meta.data)[colnames(WholeTissueHarmony@meta.data) == "CellType"] <- "Module"

DimPlot(WholeTissueHarmony, group.by = 'Module', cols = colors_MP, raster=FALSE, pt.size = 0.01)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(WholeTissueHarmony)%in%rownames(niche))
WholeTissueHarmony$niche <- niche[colnames(WholeTissueHarmony), 1]

DimPlot(WholeTissueHarmony, group.by = 'niche', cols = colors_tumor_niche, raster=FALSE, pt.size = 0.01)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

table(WholeTissueHarmony$niche, WholeTissueHarmony$Module)



 
sc_tumor <- readr::read_rds("ref_tumor_mp.rds.gz")

library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)

data <- sc_tumor@meta.data[,c('MP','Sample','MMR_Status')]

R_oe <- calTissueDist(data,
                      byPatient = F,
                      colname.cluster = "MP",
                      colname.patient = "Sample",
                      colname.tissue = "MMR_Status",
                      method = "chisq", 
                      min.rowSum = 0) 
R_oe

# col_fun <- colorRamp2(c(min(R_oe, na.rm = TRUE), 1, max(R_oe, na.rm = TRUE)), 
#                       c("#7BAFDE", "#FFFF00", "#DC050C"))
col_fun = colorRamp2(c(min(R_oe, na.rm = TRUE), 1, max(R_oe, na.rm = TRUE)), viridis(3))
ht <- Heatmap(as.matrix(R_oe),
              show_heatmap_legend = TRUE,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              show_column_dend = FALSE,
              row_names_side = 'right',
              show_column_names = TRUE,
              show_row_names = TRUE,
              col = col_fun,
              column_names_rot = 0,  
              column_names_centered = TRUE,
              row_names_gp = gpar(fontsize = 10),
              rect_gp = gpar(col = "gray", lwd = 0.5),
              column_names_gp = gpar(fontsize = 10),
              heatmap_legend_param = list(
                title = "Ro/e",
                direction = "vertical",  
                title_position = "topcenter", 
                at = seq(0, 4, by = 1),
                labels = seq(0, 4, by = 1),
                legend_height = unit(4, "cm"), 
                legend_gp = gpar(fill = col_fun(seq(0, 4, length.out = 100)))
              ),
              cell_fun = function(j, i, x, y, width, height, fill) {
                val <- R_oe[i, j]
                label <- if (val > 1) {
                  "+++"
                } else if (val > 0.8) {
                  "++"
                } else if (val >= 0.2) {
                  "+"
                } else if (val > 0) {
                  "+/-"
                } else {
                  "-"
                }
                grid.text(label, x, y, gp = gpar(fontsize = 10))
              }
)
ht


library(clusterProfiler)

sc_tumor <- readr::read_rds("ref_tumor_mp.rds.gz") %>% subset(MP %in% c('EMT-II','EMT-III'))

Idents(sc_tumor) <- "MP" 

two_groups <- c('EMT-II', 'EMT-III')

markers = FindMarkers(object = sc_tumor, ident.1 = two_groups[1], ident.2 = two_groups[2], test.use='MAST')
markers <- markers %>% tibble::rownames_to_column('gene')

go_organism <- "org.Hs.eg.db" # org.Hs.eg.db/org.Mm.eg.db

ids = bitr(markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=go_organism)
markers = merge(markers, ids, by.x='gene', by.y='SYMBOL')
markers$group <- factor(ifelse(markers$avg_log2FC < 0, -1, 1), levels = c(-1, 1))


up_genes <- subset(markers, group==1)$ENTREZID
down_genes <- subset(markers, group==-1)$ENTREZID

library(msigdbr)
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)


up <- enricher(up_genes, TERM2GENE=hallmark)
up <- as.data.frame(up) %>% dplyr::mutate(group=rep(1, n()))

down <- enricher(down_genes, TERM2GENE=hallmark)
down <- as.data.frame(down) %>% dplyr::mutate(group=rep(-1, n()))

dat <- rbind(down, up)


# Rename group information
sample_names <- c(paste0(two_groups[1]," up-regulated"), paste0(two_groups[2]," up-regulated"))
dat <- dat %>% 
  dplyr::mutate(group_type = factor(ifelse(group == 1, sample_names[1], sample_names[2]), levels = sample_names)) %>% 
  dplyr::mutate(Gene_Number = Count * group)


dat <- dat %>% dplyr::group_by(group_type) %>% dplyr::do(head(., n = 10)) 



dat$Description <- sub("^[^_]+_", "", dat$Description)  
dat$Description <- gsub("_", " ", dat$Description)     

dat$log_padjust <- -log10(dat$p.adjust)
dat$abs_gene_num <- abs(dat$Gene_Number)


dat1 <- subset(dat, group_type == sample_names[1])  # EMT-II up-regulated
dat2 <- subset(dat, group_type == sample_names[2])  # EMT-III up-regulated

# EMT-II up-regulated
p1 <- ggplot(dat1, aes(
  x = log_padjust,
  y = reorder(Description, log_padjust),
  size = abs_gene_num,
  fill = log_padjust
)) +
  geom_point(shape = 21, alpha = 0.8, stroke = 0.5) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "-log10(adj.P)",
    direction = 1
  ) +
  scale_size_continuous(range = c(3, 8), name = "Gene Number") +
  scale_x_continuous(
    limits = c(0, max(dat2$log_padjust, na.rm = TRUE) * 1.05),
    breaks = c(
      0,
      pretty(dat2$log_padjust)
    )
  ) +
  labs(
    x = "-log10(Adjusted P-value)",
    y = "Pathway",
    title = sample_names[1]
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    panel.grid.major.y = element_line(linetype = "dotted"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),  
    axis.ticks = element_line(colour = "black", linewidth = 0.3),              
    axis.ticks.length = unit(0.2, "cm")                                      
  )

# EMT-III up-regulated
p2 <- ggplot(dat2, aes(
  x = log_padjust,
  y = reorder(Description, log_padjust),
  size = abs_gene_num,
  fill = log_padjust
)) +
  geom_point(shape = 21, alpha = 0.8, stroke = 0.5) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "-log10(adj.P)",
    direction = 1
  ) +
  scale_size_continuous(range = c(3, 8), name = "Gene Number") +
  scale_x_continuous(
    limits = c(0, max(dat2$log_padjust, na.rm = TRUE) * 1.05),
    breaks = c(
      0,
      pretty(dat2$log_padjust)
    )
  ) +
  labs(
    x = "-log10(Adjusted P-value)",
    y = "Pathway",
    title = sample_names[2]
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    panel.grid.major.y = element_line(linetype = "dotted"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.2, "cm")
  )

library(patchwork)
combined_plot <- p1 / p2 + 
  plot_layout(ncol = 1, guides = "collect") & 
  theme(legend.position = "right")


# GSEA  
deg_df <- FindMarkers(sc_tumor, 
                      ident.1 = two_groups[1], 
                      ident.2 = two_groups[2],
                      logfc.threshold = -Inf, 
                      min.pct = 0)  

data <- deg_df
data$`gene` <- rownames(data)

gene <- data$gene
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

data_all <- data %>% 
  inner_join(gene, by = c("gene" = "SYMBOL"))

data_all_sort <- data_all %>% 
  arrange(desc(avg_log2FC))

geneList = data_all_sort$avg_log2FC
names(geneList) <- data_all_sort$ENTREZID 


# GSEA
library(msigdbr)
library(enrichplot)
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

gsea <- GSEA(geneList,
             TERM2GENE = hallmark)


geneSetID <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
               "HALLMARK_FATTY_ACID_METABOLISM",
               "HALLMARK_GLYCOLYSIS"
)


gseaplot2(gsea,
          # title = "EPITHELIAL MESENCHYMAL TRANSITION",  
          geneSetID = geneSetID, 
          # color = c("#d83d53", "#e1a4c6","#e3d1db"),  # EMT-II
          color = c("#1965B0", "#7BAFDE", "#c7deef"),  # EMT-III
          base_size = 15, 
          pvalue_table = T) 




library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)

niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(niche))
cosmx$niche <- niche[colnames(cosmx), 1]
mal <- subset(cosmx, subset = CellType == "TumorCells")

sub <- subset(mal, niche%in%c(0,1,3,7,8,9,10,13,14,15))
data <- sub@meta.data[,c('sample','group','niche')]
data$niche <- paste0("Niche", data$niche)
R_oe <- calTissueDist(data,
                      byPatient = F,
                      colname.cluster = "niche",
                      colname.patient = "sample",
                      colname.tissue = "group",
                      method = "chisq", 
                      min.rowSum = 0) 
R_oe

## +++, Ro/e > 1;
## ++, 0.8 < Ro/e ≤ 1;
## +, 0.2 ≤ Ro/e ≤ 0.8;
## +/−, 0 < Ro/e < 0.2;
## −, Ro/e = 0
# col_fun <- colorRamp2(c(min(R_oe, na.rm = TRUE), 1, max(R_oe, na.rm = TRUE)), 
#                       c("#7BAFDE", "#FFFF00", "#DC050C"))
col_fun = colorRamp2(c(min(R_oe, na.rm = TRUE), 1, max(R_oe, na.rm = TRUE)), viridis(3))

ht <- Heatmap(as.matrix(R_oe),
              show_heatmap_legend = TRUE, 
              cluster_rows = T, 
              cluster_columns = F,
              show_row_dend = F,
              row_names_side = 'right', 
              show_column_names = TRUE,
              show_row_names = TRUE,
              col = col_fun,
              column_names_rot = 0,  
              column_names_centered = TRUE,  # 
              row_names_gp = gpar(fontsize = 10),
              rect_gp = gpar(col = "gray", lwd = 0.5),  
              column_names_gp = gpar(fontsize = 10),
              heatmap_legend_param = list(
                title = "Ro/e",
                direction = "horizontal",
                title_position = "leftcenter",
                at = seq(0, 4, by = 1),
                labels = seq(0, 4, by = 1),
                legend_gp = gpar(fill = col_fun(R_oe))
              ),
              cell_fun = function(j, i, x, y, width, height, fill) {
                val <- R_oe[i, j]
                label <- if (val > 1) {
                  "+++"
                } else if (val > 0.8) {
                  "++"
                } else if (val >= 0.2) {
                  "+"
                } else if (val > 0) {
                  "+/-"
                } else {
                  "-"
                }
                grid.text(label, x, y, gp = gpar(fontsize = 10))
              }
)
draw(
  ht,
  heatmap_legend_side = "bottom" 
)



cosmx <- readr::read_rds('merge_data_final.rds.gz')

niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(niche))
cosmx$niche <- niche[colnames(cosmx), 1]
mal <- subset(cosmx, subset = CellType == "TumorCells")

mal <- subset(mal, niche%in%c(0,1,3,7,8,9,10,13,14,15))


WholeTissueHarmony <- readr::read_rds('malignant.rds.gz')

niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(WholeTissueHarmony)%in%rownames(niche))
WholeTissueHarmony$niche <- niche[colnames(WholeTissueHarmony), 1]

sub <- subset(WholeTissueHarmony, niche%in%c(0,1,3,7,8,9,10,13,14,15)) # 
colnames(sub@meta.data)[colnames(sub@meta.data) == "CellType"] <- "Module"

ClusterFreq <- sub@meta.data[,c("Module","niche")] %>% table %>%
  data.frame() %>% set_colnames(c("Module","niche","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(Module,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=Module,value=Per,-niche)
ClusterPer$niche=factor(ClusterPer$niche,levels = names(table(sub$niche)))

p1 <- ggplot(data = ClusterPer, aes(x = niche, fill=Module,y=Per)) + #log2(as.numeric(CopyNumber)))
  geom_bar(stat = "identity",width = 0.6)+
  scale_y_continuous(expand=c(0,0))+labs(y="percentage(%)")+
  # coord_flip()+
  scale_fill_manual(values=colors_MP)+
  theme(axis.text=element_text(color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black",size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10))


ClusterFreq <- sub@meta.data[,c("group","niche")] %>% table %>%
  data.frame() %>% set_colnames(c("group","niche","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(group,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=group,value=Per,-niche)
ClusterPer$niche=factor(ClusterPer$niche,levels = names(table(sub$niche)))
ClusterPer$group=factor(ClusterPer$group,levels = names(table(sub$group)))

p2 <- ggplot(data = ClusterPer, aes(x = niche, fill=group,y=Per)) + #log2(as.numeric(CopyNumber)))
  geom_bar(stat = "identity",width = 0.6)+
  scale_y_continuous(expand=c(0,0))+labs(y="percentage(%)")+
  # coord_flip()+
  scale_fill_manual(values=colors)+
  theme(axis.text=element_text(color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5,size = 15),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), 
        axis.text.y=element_text(color = "black",size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10))


p1+p2


#### Sample-Niche ####
cosmx <- readr::read_rds('merge_data_final.rds.gz') 

niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(niche))
cosmx$niche <- niche[colnames(cosmx), 1]

ClusterFreq <- cosmx@meta.data[,c("niche","sample")] %>% table %>%
  data.frame() %>% set_colnames(c("niche","sample","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(niche,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=niche,value=Per,-sample)
# ClusterPer$sample=factor(ClusterPer$sample,levels = names(table(cosmx$sample)))
ClusterPer$sample=factor(ClusterPer$sample,levels = c(paste0("MSI_N", 1:6),paste0("MSI_T", c(1:13)),paste0("MSS_N", 1:14),paste0("MSS_T", c(1:24))))
ClusterPer$niche=factor(ClusterPer$niche,levels = names(table(cosmx$niche)))

library(dplyr)
library(stringr)
sample_order <- cosmx@meta.data %>%
  select(sample, group) %>%
  distinct() %>%
  mutate(
    sample_num = as.numeric(str_extract(sample, "\\d+"))  
  ) %>%
  arrange(group, sample_num) %>% 
  pull(sample) 

ClusterPer$sample <- factor(ClusterPer$sample, 
                            levels = sample_order)

ClusterPer$group <- cosmx@meta.data[match(ClusterPer$sample, cosmx@meta.data$sample), "group"]

group_colors <- c("MSI_N" = "#E41A1C", "MSI_T" = "#377EB8", "MSS_N" = "#4DAF4A", "MSS_T" = "#984EA3")

p <- ggplot(data = ClusterPer, aes(x = sample, fill = niche, y = Per)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 100, by = 20),
                     limits = c(0, 100)) +  
  labs(y = "Percentage(%)") +
  scale_fill_manual(values = MyColorCode) +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.line = element_line(color = "black"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.margin = unit(c(1,1,2,1), "lines") 
  )

for(i in seq_along(levels(ClusterPer$sample))){
  sample <- levels(ClusterPer$sample)[i]
  group <- unique(ClusterPer$group[ClusterPer$sample == sample])
  p <- p + annotation_custom(
    grob = rectGrob(gp = gpar(fill = group_colors[group], col = NA)),
    xmin = i-0.5, xmax = i+0.5,
    ymin = 101, ymax = 102  
  )
}

p + coord_cartesian(ylim = c(0, 100), clip = "off")






WholeTissueHarmony <- readr::read_rds('malignant.rds.gz')

ClusterFreq <- WholeTissueHarmony@meta.data[,c("MP","sample")] %>% table %>%
  data.frame() %>% set_colnames(c("MP","sample","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(MP,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=MP,value=Per,-sample)
ClusterPer$sample=factor(ClusterPer$sample,levels = names(table(WholeTissueHarmony$sample)))
ClusterPer <- ClusterPer %>% mutate(group = ifelse(grepl("MSI", sample), "MSI_T", "MSS_T"))

sum(ClusterPer$Per[ClusterPer$sample == 'MSI_T5'])

# color <- c("MSI_T" = "#E41A1C", "MSS_T" = "#377EB8")
color <- c("MSI_T" = "#9DD1D0", "MSS_T" = "#E36B69")


p_df <- ClusterPer %>%
  group_by(MP) %>%
  summarise(
    p = wilcox.test(Per ~ group)$p.value
  )
p_df$padj <- p.adjust(p_df$p, method = "BH")




ggplot(ClusterPer, aes(x = MP, y = Per, fill = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
  #             size = 0.5, alpha = 0.6) +
  labs(y = "Percentage (%)", x = "Cell Type") +
  scale_fill_manual(values = color) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test", 
    label = "p.signif", 
    size = 6,
    hide.ns = TRUE  
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"), 
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.title = element_blank()
  )