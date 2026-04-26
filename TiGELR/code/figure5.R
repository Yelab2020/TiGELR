
#### cellchat ####
library(CellChat)
library(patchwork)
library(magrittr)
library(Seurat)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 800 * 1024^2) 

#### Part I: Data input & processing and initialization of CellChat object
data <- readr::read_rds('merge_data_final.rds.gz')
data$sub_subtype[data$CellType == 'TumorCells'] <- 'TumorCells'
data$sub_subtype[data$CellType == 'TCells'] <- data$subtype[data$CellType == 'TCells']
table(data$sub_subtype)

tumor <- readr::read_rds('malignant.rds.gz')

data@meta.data[colnames(tumor), 'sub_subtype'] <- tumor@meta.data[colnames(tumor), 'CellType']

MSS <- subset_opt(data, group == 'MSS_T' & space_area %in% c('TumorCore','Boundary')) %>% subset_opt(sub_subtype %in% c('EMT-II','Glia','Tip Cells'))
MSS$sub_subtype <- paste0(MSS$sub_subtype, '_MSS')
MSI <- subset_opt(data, group == 'MSI_T' & space_area %in% c('TumorCore','Boundary')) %>% subset_opt(sub_subtype %in% c('EMT-II','Glia','Tip Cells'))
MSI$sub_subtype <- paste0(MSI$sub_subtype, '_MSI')


# Start
seu <- merge(MSS, MSI)

sample_counts <- table(seu$sample)
valid_samples <- names(sample_counts[sample_counts >= 5])
seu <- subset_opt(seu, sample %in% valid_samples)

samples <- unique(seu$sample)
data_list <- list()
meta_data <- seu@meta.data
for (sample in samples) {
  cells_in_sample <- rownames(meta_data[meta_data$sample == sample, ])
  data_list[[sample]] <- subset_opt(seu, cells = cells_in_sample)
}

Idents(seu) <- 'sub_subtype'

data.input <- NULL
meta<- NULL
spatial.locs <- NULL
spatial.factors <- NULL

for (sample_name in names(data_list)) {
  print(sample_name)
  
  seu_sub <- data_list[[sample_name]]
  Idents(seu_sub) <- 'sub_subtype'
  
  set.seed(123)
  VGENES=rownames(seu_sub)
  VGENES=setdiff(VGENES,VGENES[grep("^SystemControl",VGENES)])
  VGENES = VGENES[!grepl("Neg|NegProbe", VGENES)]
  seu_sub <- seu_sub[VGENES,]
  
  seu_sub <- NormalizeData(seu_sub)
  seu_sub <- FindVariableFeatures(object =seu_sub, mean.function = ExpMean, dispersion.function = LogVMR)
  seu_sub <- ScaleData(object = seu_sub)
  
  data.input1 <- Seurat::GetAssayData(seu_sub, slot = "data", assay = "Nanostring")
  meta1 <- data.frame(labels = Idents(seu_sub), samples = names(data_list[sample_name]))
  
  spatial.locs1 <- Seurat::GetTissueCoordinates(seu_sub, scale = NULL, cols = c("imagerow", "imagecol"))
  spatial.locs1 <- spatial.locs1[, 1:2]
  conversion.factor = 0.12028  # CosMx
  d = computeCellDistance(spatial.locs1)
  spot.size = min(d)*conversion.factor # converting the distance in Pixels to Micrometers
  spatial.factors1 = data.frame(ratio = conversion.factor, tol = spot.size/2)
  
  
  if (is.null(data.input)) {
    data.input <- data.input1
  } else {
    data.input <- cbind(data.input, data.input1)
  }
  
  if (is.null(meta)) {
    meta <- meta1
  } else {
    meta <- rbind(meta, meta1)
  }
  
  if (is.null(spatial.locs)) {
    spatial.locs <- spatial.locs1
  } else {
    spatial.locs <- rbind(spatial.locs, spatial.locs1)
  }
  
  if (is.null(spatial.factors)) {
    spatial.factors <- spatial.factors1
  } else {
    spatial.factors <- rbind(spatial.factors, spatial.factors1)
  }
}

rownames(spatial.locs) <- colnames(data.input)
rownames(spatial.factors) <- samples
meta$labels <- factor(meta$labels, levels = levels(Idents(seu)))
meta$samples <- factor(meta$samples, levels = samples)

unique(meta$labels) # check the cell labels
unique(meta$samples) # check the sample labels


cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels", 
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
cellchat

CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 5) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#### Part II: Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = F, interaction.range = 250, scale.distance = 0.2,  # distance.use = FALSE,过滤掉距离远的信号，但是不添加距离限制
                              contact.dependent = TRUE, contact.range = 10)
# cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
#                               distance.use = TRUE, interaction.range = 250, 
#                               scale.distance = 4.5,
#                               contact.dependent = TRUE, contact.range = 20)
cellchat <- filterCommunication(cellchat, min.cells = 10) # 

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

levels(cellchat@idents)

netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds",color.use = c('red', '#dbad5f', '#7db3af', 'red','#dbad5f',
                                                                                    '#7db3af'))

color.use <- c(
  "EMT-II_MSS" = "red",
  "Tip Cells_MSS" = "#dbad5f",
  "Glia_MSS" = "#7db3af",
  "EMT-II_MSI" = "red",
  "Tip Cells_MSI" = "#dbad5f",
  "Glia_MSI" = "#7db3af"
)

netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds",color.use = c('red', '#dbad5f', '#7db3af', 'red','#dbad5f',
                                                                                    '#7db3af'))


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions",
                 color.use = color.use)


# scatter
num.link <- rowSums(cellchat@net$count) + colSums(cellchat@net$count)-diag(cellchat@net$count)

weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets

netAnalysis_signalingRole_scatter(cellchat, title = names(cellchat), weight.MinMax = weight.MinMax)



# outgoing / incoming
# Compute and visualize the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 12, height = 15,cluster.cols = T,
                                         signaling = c())
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 12, height = 15, cluster.cols = T,
                                         signaling = c())
ht1 + ht2


# bubble
df.net <- subsetCommunication(cellchat, thresh = 0.05)  # 0.05
head(df.net)

unique(df.net$pathway_name)

pairLR.use <- extractEnrichedLR(cellchat, signaling = unique(df.net$pathway_name))



# refine plot
# cosmx
cellchat_st <- readRDS("/Cellchat_MSS_merge.rds")
levels(cellchat_st@idents)
df.net_st <- subsetCommunication(cellchat_st, thresh = 0.05) %>%
  .[(.$target == "EMT-II_MSS") & (.$source %in% c("Glia_MSS")), ]  # "Glia_MSS", "Tip Cells_MSS"
head(df.net_st)
unique(df.net_st$pathway_name)
pairLR.use_st <- extractEnrichedLR(cellchat_st, signaling = unique(df.net_st$pathway_name))

# sc
cellchat <- readRDS("Cellchat_sc_MSS.rds")
levels(cellchat@idents)
df.net <- subsetCommunication(cellchat, thresh = 0.05) %>%
  .[(.$target == "EMT-II") & (.$source %in% c("Glia")), ]   # "Glia", "Tip Cells"
head(df.net)
unique(df.net$pathway_name)
pairLR.use <- extractEnrichedLR(cellchat, signaling = unique(df.net$pathway_name))

pairLR_1 <- intersect(pairLR.use$interaction_name, pairLR.use_st$interaction_name)
pairLR_1 <- as.data.frame(pairLR_1)
colnames(pairLR_1) <- 'interaction_name'

pairLR <- readRDS('common_pairLR.rds')
levels(cellchat@idents)   
p <- netVisual_bubble(cellchat, remove.isolate = FALSE, sources.use = c("Glia_MSS","Tip Cells_MSS"),targets.use = c("EMT-II_MSS"),
                      thresh = 0.001, grid.on=T,color.grid = "black",
                      pairLR.use = pairLR)


