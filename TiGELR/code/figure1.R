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

#### Figure 1 ####
cosmx <- readr::read_rds('merge_data_final.rds.gz')
cosmx$CellType[cosmx$CellType == 'TumorCells'] <- 'EpithelialCells'

cosmx$CellType <- factor(cosmx$CellType, levels =  c('EndothelialCells','EpithelialCells','Glia','Fibroblasts','MastCells','DC','Macrophages',
                                                                     'Neutrophils','BCells','PlasmaCells','TCells') )

MSI_T <- subset(cosmx, subset = sample == 'MSI_T9')

MSI_T$CellType <- factor(MSI_T$CellType, levels =  c('EndothelialCells','EpithelialCells','Glia','Fibroblasts','MastCells','DC','Macrophages',
                                                                     'Neutrophils','BCells','PlasmaCells','TCells') )

DefaultBoundary(MSI_T[["A_CRC.1"]]) <- "segmentation"
ImageDimPlot(MSI_T, group.by = 'CellType', size = 0.1, cols = CellType_colors, border.size = NA)+
  labs(title = "MSI_T") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

#### Zoomed
range(MSI_T$CenterX_global_px)
range(MSI_T$CenterY_global_px)

CellType_colors <- c(
  EndothelialCells = '#ffe119',
  EpithelialCells  = '#48a854',
  Glia             = '#7db3af',
  Fibroblasts      = '#bfef45',
  MastCells        = '#dcbeff',
  DC               = '#e97371',
  Macrophages      = '#911eb4',
  Neutrophils      = '#f2d3ca',
  BCells           = '#42d4f4',
  PlasmaCells      = '#f032e6',
  TCells           = '#e8a767'
)



sub <- subset(MSI_T, CenterX_global_px > 68000 &  CenterX_global_px < 71000 & CenterY_global_px > 7000 & CenterY_global_px < 10000)

p1 <- ImageDimPlot(sub, fov = "A_CRC.1",  cols = CellType_colors, border.size = NA, group.by = 'CellType')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1

p2 <- ImageDimPlot(MSI_T, fov = "A_CRC.1", cols = CellType_colors, border.size = NA, size = 1, group.by = 'CellType')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3 <- p2 + annotate(
  "rect",
  # xmin = 54000, xmax = 57000, # MSI_N2
  # ymin = 54000, ymax = 57000,
  # xmin = 109000, xmax = 112000, # MSS_N5
  # ymin = 34000, ymax = 37000,
  # xmin = 70000, xmax = 73000, # MSS_T6
  # ymin = 67500, ymax = 70500,
  xmin = 7000, xmax =10000, # MSI_T9
  ymin = 68000, ymax = 71000,
  color = "white",
  fill = NA,
  size = 0.8,
  # linetype = "dashed"  
)
print(p3)

p3+p1




ClusterFreq <- cosmx@meta.data[,c("CellType","sample")] %>% table %>%
  data.frame() %>% set_colnames(c("CellType","sample","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(CellType,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=CellType,value=Per,-sample)

ClusterPer$group <- cosmx@meta.data[match(ClusterPer$sample, cosmx@meta.data$sample), "group"]
ClusterPer$sample <- factor(ClusterPer$sample, levels = names(table(cosmx$sample)))
group_colors <- c("MSI_N" = "#E41A1C", "MSI_T" = "#377EB8", "MSS_N" = "#4DAF4A", "MSS_T" = "#984EA3")

ClusterPer$CellType <- factor(ClusterPer$CellType, 
                                      levels = c('EndothelialCells','EpithelialCells','Glia','Fibroblasts','MastCells','DC','Macrophages',
                                                 'Neutrophils','BCells','PlasmaCells','TCells'))
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


p <- ggplot(data = ClusterPer, aes(x = sample, fill = CellType, y = Per)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +  
  labs(y = "Percentage(%)") +
  scale_fill_manual(values = CellType_colors) +
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


ClusterFreq <- cosmx@meta.data[,c("CellType","group")] %>% table %>%
  data.frame() %>% set_colnames(c("CellType","group","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(CellType,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=CellType,value=Per,-group)
ClusterPer$group=factor(ClusterPer$group,levels = names(table(cosmx$group)))


ClusterPer$CellType <- factor(ClusterPer$CellType, 
                                      levels = c('EndothelialCells','EpithelialCells','Glia','Fibroblasts','MastCells','DC','Macrophages',
                                                 'Neutrophils','BCells','PlasmaCells','TCells'))

ggplot(data = ClusterPer, aes(x = group, fill=CellType,y=Per)) + #log2(as.numeric(CopyNumber)))
  geom_bar(stat = "identity",width = 0.6)+
  scale_y_continuous(expand=c(0,0))+labs(y="Percentage(%)")+
  # coord_flip()+
  scale_fill_manual(values=CellType_colors)+
  theme(axis.text=element_text(color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black"),
        axis.title.x = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10))



ClusterFreq <- cosmx@meta.data[,c("CellType","sample")] %>% filter(grepl("T", sample) & CellType != "EpithelialCells")  %>% table %>%
  data.frame() %>% set_colnames(c("CellType","sample","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(CellType,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=CellType,value=Per,-sample)
ClusterPer$sample=factor(ClusterPer$sample,levels = names(table(cosmx$sample)))
ClusterPer <- ClusterPer %>% mutate(group = ifelse(grepl("MSI", sample), "MSI_T", "MSS_T"))

sum(ClusterPer$Per[ClusterPer$sample == 'MSI_T5'])


color <- c("MSI_T" = "#E41A1C", "MSS_T" = "#377EB8")

ggplot(ClusterPer, aes(x = CellType, y = Per, fill = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
              size = 1, alpha = 0.6) +
  labs(y = "Percentage (%)", x = "Cell Type") +
  scale_fill_manual(values = color) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test", 
    label = "p.signif", 
    hide.ns = TRUE   
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_blank()
  )



ClusterFreq <- cosmx@meta.data[, c("sub_subtype", "sample")] %>% filter(grepl("T", sample) & sub_subtype != "EpithelialCells")  %>% table %>% 
  data.frame() %>% set_colnames(c("sub_subtype","sample","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(sub_subtype,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=sub_subtype,value=Per,-sample)
ClusterPer$sample=factor(ClusterPer$sample,levels = names(table(cosmx$sample)))
ClusterPer <- ClusterPer %>% mutate(group = ifelse(grepl("MSI", sample), "MSI_T", "MSS_T"))

sum(ClusterPer$Per[ClusterPer$sample == 'MSI_T5'])


color <- c("MSI_T" = "#E41A1C", "MSS_T" = "#377EB8")

ggplot(ClusterPer, aes(x = sub_subtype, y = Per, fill = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
              size = 1, alpha = 0.6) +
  labs(y = "Percentage (%)", x = "Cell Type") +
  scale_fill_manual(values = color) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test", 
    label = "p.signif", 
    hide.ns = TRUE  
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_blank()
  )




#### UMAP ####
merge_data <- readr::read_rds('merge_data_final.rds.gz')
merge_data$CellType[merge_data$CellType == 'TumorCells'] <- 'EpithelialCells'

merge_data$CellType <- factor(merge_data$CellType, levels =  c('EndothelialCells','EpithelialCells','Glia','Fibroblasts','MastCells','DC','Macrophages',
                                                                               'Neutrophils','BCells','PlasmaCells','TCells') )

set.seed(123)
VGENES=rownames(merge_data)
VGENES=setdiff(VGENES,VGENES[grep("^SystemControl",VGENES)])
VGENES = VGENES[!grepl("Neg|NegProbe", VGENES)]
merge_data <- merge_data[VGENES,]
merge_data <- SCTransform(merge_data, assay = "Nanostring", clip.range = c(-10, 10))
merge_data <- RunPCA(merge_data, npcs = 50)
WholeTissueHarmony <- merge_data %>% 
  RunHarmony("Run_Tissue_name", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(WholeTissueHarmony, 'harmony')
WholeTissueHarmony <- WholeTissueHarmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>%    # dims = 1:30
  FindNeighbors(reduction = "harmony", dims = 1:20)  # dims = 1:30

DimPlot(WholeTissueHarmony, raster = F, cols = CellType_colors,  group.by = 'CellType')

library(tidydr)
library(ggplot2)

p1 <- DimPlot(WholeTissueHarmony, raster = F, cols = CellType_colors,  group.by = 'CellType') + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        plot.title = element_blank()) 


p2 <- DimPlot(WholeTissueHarmony, raster = F, cols = colors,  group.by = 'sample') + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        plot.title = element_blank()) 

p1+p2




group_colors <- c("MSI_N" = "#E41A1C", "MSI_T" = "#377EB8", "MSS_N" = "#4DAF4A", "MSS_T" = "#984EA3")


p3 <- DimPlot(WholeTissueHarmony, raster = F, cols = group_colors,  group.by = 'group') + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        plot.title = element_blank()) 

p3





markers <- list(
  TCells = c("CD3D", "CD3G", "CD8A","CD3E",'CD2'),
  PlasmaCells = c("MZB1", 'JCHAIN'),
  BCells = c("CD79A","MS4A1", "CD19"),
  Neutrophils = c("CSF3R",'CXCL8','G0S2'),
  Macrophages = c("CD68", "CD163",'APOE','CD14'), 
  DC = c("CD1C",'CCL22','CLEC10A'),
  MastCells = c("CPA3", "KIT"),
  Fibroblasts = c("COL1A1", "COL3A1", "DCN"),
  Glia = c("S100B","NRXN1","SOX2"),
  EpithelialCells = c("EPCAM", "CDH1","KRT19"),
  EndothelialCells = c("PECAM1", "VWF", "CD34")
)

DotPlot(merge_data, group.by = 'CellType',features = markers, scale = TRUE) + 
  RotatedAxis() +
  scale_color_viridis_c() + 
  scale_size_continuous(range = c(0.1, 7)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

