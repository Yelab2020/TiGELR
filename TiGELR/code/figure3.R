MyColorCode <- c( '#48a854', '#ffe119', '#526cc9', '#f58231','#d83d53',
                  '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4',
                  '#469990', '#dcbeff', '#e8a767', '#f1e7ae', '#cc5500',
                  '#aaffc3', '#808000', '#ffd8b1', '#6c51a3', '#d3d3d3')

colors <- c('#dc8e97','#e3d1db','#74a893','#ac9141','#5ac6e9','#e5c06e','#7587b1','#c7deef','#e97371','#e1a4c6',
            '#916ba6','#cb8f82','#7db3af','#d2e0ac','#f3f3b0','#a5d177','#e0bc58','#64abc0','#fab37f','#e98741',
            '#8fc0dc','#967568','#f2d3ca','#eebd85','#82c785','#edeaa4','#cdaa9f','#794976','#bcacd3','#889b5d',
            '#4e9592','#dbad5f','#64ae79','#ac5092','#f37d74','#4a6eb6','#bbe3ee','#5b9951',
            '#fbd9ae','#cf622b','#90cdc1','#f4eb9b','#5f9d58','#f5949f','#c1cf9a','#70afab','#a8c375','#8e95c9',
            '#efcfee','#5384c4','#cacae3','#e45878','#996da8','#3b99d4','#edf5e1')

color_niche <- c('#526cc9','#d83d53','#911eb4','#42d4f4','#dcbeff', '#e8a767')


library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)  
library(magrittr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(harmony)
library(rlang)

cosmx <- readr::read_rds('merge_data_final.rds.gz')

cosmx <- subset(cosmx, tissue == 'T')

niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(niche))
cosmx$niche <- niche[colnames(cosmx), 1]

sub <- subset(cosmx, niche%in%c(2,4,5,6,11,12))

df <- sub@meta.data[,c('sub_subtype','niche')]
df$niche <- paste0("Niche", df$niche)
prop_table <- prop.table(table(df$niche, df$sub_subtype), margin = 1)
prop_matrix <- as.data.frame.matrix(t(prop_table))

row_means <- rowMeans(prop_matrix)
row_sds <- apply(prop_matrix, 1, sd)

prop_matrix_zscaled <- prop_matrix

library(ComplexHeatmap)
library(viridis)
Heatmap(
  prop_matrix_zscaled,
  name = "Z-score",  # Composition
  # col = colorRampPalette(c("#7BAFDE", "#FFFF00", "#DC050C"))(100),
  col = colorRamp2(seq(0, 1, length.out = 100), viridis(100)),
  row_names_side = "right",
  # column_names_side = "top",
  column_names_rot = 90,
  cluster_rows = TRUE,
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

#### enrichment ####
library(clusterProfiler)
library(org.Hs.eg.db)  
library(dplyr)
library(tibble)
library(enrichplot)
library(ComplexHeatmap)

niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(niche))
cosmx$niche <- niche[colnames(cosmx), 1]

sub <- subset(cosmx, niche%in%c(2,4,5,6,11,12))

Idents(sub) <- "niche" 


# COSG
library(COSG)
library(patchwork)
library(viridis)
library(tidyr)

set.seed(123)
VGENES=rownames(sub)
VGENES=setdiff(VGENES,VGENES[grep("^SystemControl",VGENES)])
VGENES = VGENES[!grepl("Neg|NegProbe", VGENES)]
sub <- sub[VGENES,]
sub <- SCTransform(sub, assay = "Nanostring", clip.range = c(-10, 10))

markers <- cosg(
  sub,
  groups='all', # table(Idents(sub))
  assay='SCT',
  slot='data',
  mu=1,       
  remove_lowly_expressed=TRUE,  
  expressed_pct=0.1,            
  n_genes_user=50   
)


marker_names <- markers$names

marker_df <- marker_names %>%
  pivot_longer(
    cols = everything(),
    names_to = "cluster",
    values_to = "gene"
  ) %>%
  filter(!is.na(gene))  

new.cluster.ids <- c("5","4","12","2","11","6")
names(new.cluster.ids) <- levels(sub)
sub <- RenameIdents(sub, new.cluster.ids)



#### Boundary ####
cosmx <- readr::read_rds('merge_data_final.rds.gz')
niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(niche))
cosmx$niche <- niche[colnames(cosmx), 1]
boundary <- as.data.frame(read.csv('DBSCAN_Boundary.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(boundary))
cosmx$space_area <- boundary[colnames(cosmx), 12]

MSS_T <- subset(cosmx, sample == 'MSS_T16')
MSS_T$space_area <- factor(MSS_T$space_area, levels =  c('TumorCore','Boundary','Distant'))
DefaultBoundary(MSS_T[["A_CRC.1"]]) <- "segmentation"
ImageDimPlot(MSS_T, group.by = 'space_area', size = 1, cols = c('#DC050C','#1965B0','#48a854'), border.size = NA)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


sub <- subset(cosmx, tissue =='T')
niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(sub)%in%rownames(niche))
sub$niche <- niche[colnames(sub), 1]
sub <- subset(sub, niche%in%c(2,4,5,6,11,12))

library(reshape2)
library(tidyverse)
library(scales)
library(readxl)
library(dplyr)
library(tidyr)
library(ggalluvial)

plotC <- table(sub@meta.data$space_area, sub@meta.data$niche) %>% melt()
colnames(plotC) <- c("space_area", "niche", "Number")

plotC <- plotC %>%
  group_by(space_area) %>%
  mutate(percent = Number / sum(Number))
plotC$niche <- as.factor(plotC$niche)
plotC$space_area <- factor(plotC$space_area,levels = c("TumorCore", "Boundary", "Distant"))

ggplot(plotC, aes(
  x = space_area,
  y = percent,
  fill = niche,
  stratum = niche,
  alluvium = niche
)) +
  geom_flow(width = 0.5, alpha = 0.4, knot.pos = 0.5) +
  geom_stratum(width = 0.5, color = "black") +
  scale_y_continuous(labels = percent_format()) +
  theme_classic() +
  theme(axis.text=element_text(color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black"),
        axis.title.x = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10)) +
  labs(y = "Niche Proportion (%)", x = NULL) +
  scale_fill_manual(values = color_niche)


#### Boxplot niche-group ####
sub <- subset(cosmx, tissue =='T')
sub <- subset(sub, niche%in%c(2,4,5,6,11,12))
ClusterFreq <- sub@meta.data[,c("niche","sample")] %>% filter(grepl("T", sample))  %>% table %>%
  data.frame() %>% set_colnames(c("niche","sample","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(niche,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=niche,value=Per,-sample)
ClusterPer$sample=factor(ClusterPer$sample,levels = names(table(sub$sample)))
ClusterPer <- ClusterPer %>% mutate(group = ifelse(grepl("MSI", sample), "MSI_T", "MSS_T"))
sum(ClusterPer$Per[ClusterPer$sample == 'MSI_T5'])
ClusterPer$niche <- factor(ClusterPer$niche, levels = as.character(c(2,4,5,6,11,12)))

color <- c("MSI_T" = "#E41A1C", "MSS_T" = "#377EB8")
ggplot(ClusterPer, aes(x = niche, y = Per, fill = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
              size = 0.5, alpha = 1) +
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
    axis.text.x = element_text(angle = 0, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_blank())

#### Boxplot boundary niche-group ####
cosmx <- readr::read_rds('merge_data_final.rds.gz')
niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(niche))
cosmx$niche <- niche[colnames(cosmx), 1]
boundary <- as.data.frame(read.csv('DBSCAN_Boundary.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(boundary))
cosmx$space_area <- boundary[colnames(cosmx), 12]

sub <- subset(cosmx,  space_area =='Boundary'& tissue =='T')
sub <- subset(sub,  subset = sample %in% names(which(table(sub$sample) >= 100)))

sub_boundary <- subset(sub, niche%in%c(2,4,5,6,11,12))


ClusterFreq <- sub_boundary@meta.data[,c("niche","sample")] %>% table %>%
  data.frame() %>% set_colnames(c("niche","sample","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(niche,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=niche,value=Per,-sample)
ClusterPer$sample=factor(ClusterPer$sample,levels = names(table(sub_boundary$sample)))
ClusterPer <- ClusterPer %>% mutate(group = ifelse(grepl("MSI", sample), "MSI_T", "MSS_T"))
sum(ClusterPer$Per[ClusterPer$sample == 'MSI_T5'])
ClusterPer$niche <- factor(ClusterPer$niche, levels = as.character(c(2,4,5,6,11,12)))

# color <- c("MSI_T" = "#E41A1C", "MSS_T" = "#377EB8")
color <- c("MSI_T" = "#9DD1D0", "MSS_T" = "#E36B69")

ggplot(ClusterPer, aes(x = niche, y = Per, fill = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
  #             size = 0.5, alpha = 1) +
  labs(y = "Percentage (%)", x = "Niche", title = "Boundary") +
  scale_fill_manual(values = color) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test", 
    label = "p.format", 
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
    legend.title = element_blank(),
    plot.title = element_text(size = 12,hjust = 0.5)
  )

# FDR
p_df <- ClusterPer %>%
  group_by(niche) %>%
  summarise(
    p = wilcox.test(Per ~ group)$p.value
  )
p_df$padj <- p.adjust(p_df$p, method = "BH")


#
cosmx <- readr::read_rds('merge_data_final.rds.gz')
boundary <- as.data.frame(read.csv('DBSCAN_Boundary.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(boundary))
cosmx$space_area <- boundary[colnames(cosmx), 12]

niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(niche))
cosmx$niche <- niche[colnames(cosmx), 1]

sub <- subset(cosmx, space_area == 'Boundary' & tissue == 'T')

niche <- as.data.frame(read.csv('niche.csv', row.names = 1))
all(colnames(sub)%in%rownames(niche))
sub$niche <- niche[colnames(sub), 1]

sub_1 <- subset(sub, niche == '5')


ClusterFreq <- sub_1@meta.data[,c("sub_subtype","sample")] %>% filter(grepl("T", sample))  %>% table %>%
  data.frame() %>% set_colnames(c("sub_subtype","sample","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(sub_subtype,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=sub_subtype,value=Per,-sample)
ClusterPer$sample=factor(ClusterPer$sample,levels = names(table(sub_1$sample)))
ClusterPer <- ClusterPer %>% mutate(group = ifelse(grepl("MSI", sample), "MSI_T", "MSS_T"))

sum(ClusterPer$Per[ClusterPer$sample == 'MSI_T5'])

color <- c("MSI_T" = "#E41A1C", "MSS_T" = "#377EB8")

ggplot(ClusterPer, aes(x = sub_subtype, y = Per, fill = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
              size = 0.5, alpha = 0.6) +
  labs(y = "Percentage (%)", x = "Cell Type", title = 'Boundary') +
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
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.title = element_blank()
  )