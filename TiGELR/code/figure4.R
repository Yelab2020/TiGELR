library(ggplot2)
library(ggpubr)

cosmx <- readr::read_rds('merge_data_final.rds.gz')

sub <- subset(cosmx, space_area == 'Boundary' & tissue == 'T')

ClusterFreq <- sub@meta.data[,c("sub_subtype","sample")] %>% filter(grepl("T", sample))  %>% table %>%
  data.frame() %>% set_colnames(c("sub_subtype","sample","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(sub_subtype,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=sub_subtype,value=Per,-sample)
ClusterPer$sample=factor(ClusterPer$sample,levels = names(table(sub$sample)))
ClusterPer <- ClusterPer %>% mutate(group = ifelse(grepl("MSI", sample), "MSI_T", "MSS_T"))

sum(ClusterPer$Per[ClusterPer$sample == 'MSI_T5'])

color <- c("MSI_T" = "#9DD1D0", "MSS_T" = "#E36B69")


ggplot(ClusterPer, aes(x = sub_subtype, y = Per, fill = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
  #             size = 0.5, alpha = 0.6) +
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
    axis.ticks = element_line(color = "black"),  # 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.title = element_blank()
  )


sc <- readr::read_rds('MergeFilterMainTypes20241227.rds.gz')

merged_seurat <- readRDS('sc_data.rds')
table(merged_seurat$DefineTypes)
table(merged_seurat$MainTypes)


merged_seurat$new_celltype <- merged_seurat$DefineTypes
merged_seurat$new_celltype[rownames(tumor@meta.data)] <- tumor$MP[rownames(tumor@meta.data)]
table(merged_seurat$new_celltype[merged_seurat$MainTypes == 'EpithelialCells'])
table(merged_seurat$DefineTypes[merged_seurat$MainTypes == 'EpithelialCells'])

merged_seurat <- subset(merged_seurat, MainTypes != 'EpithelialCells')
tumor$DefineTypes <- tumor$MP

tumor <- readr::read_rds("ref_tumor_mp.rds.gz") 
table(tumor$MP)
tumor$Tissue <- 'Tumor'
tumor$DefineTypes <- tumor$MP


merged_seurat_final <- merge(merged_seurat, tumor)
table(merged_seurat_final$DefineTypes)

meta <- readr::read_rds('MergeMetaDataNewMainTypes20241016.rds.gz')
meta <- meta[!is.na(meta$MMR_Status), ]
meta <- meta[meta$Tissue == 'Tumor', ]

sc <- subset(sc, Patient %in% unique(meta$Patient))
sc$MMR <- meta$MMR_Status[match(sc$Patient, meta$Patient)]

data <- sc@meta.data[,c('Patient', 'MMR', 'DefineTypes')]

ClusterFreq <- data[,c("DefineTypes","Patient")] %>% table %>%
  data.frame() %>% set_colnames(c("DefineTypes","Patient","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(DefineTypes,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=DefineTypes,value=Per,-Patient)
ClusterPer$Patient=factor(ClusterPer$Patient,levels = names(table(data$Patient)))

ClusterPer$MMR <- data$MMR[match(ClusterPer$Patient, data$Patient)]

sum(ClusterPer$Per[ClusterPer$Patient == 'C103'])

# color <- c("MSI_T" = "#E41A1C", "MSS_T" = "#377EB8")
color <- c("MSI" = "#9DD1D0", "MSS" = "#E36B69")

ClusterPer <- ClusterPer[ClusterPer$DefineTypes %in% c('Glia','Tip Cells'),]


ggplot(ClusterPer, aes(x = DefineTypes, y = Per, fill = MMR)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
  #             size = 0.5, alpha = 0.6) +
  labs(y = "Percentage (%)", x = "Cell Type") +
  scale_fill_manual(values = color) +
  stat_compare_means(
    aes(group = MMR),
    method = "wilcox.test", 
    label = "p.signif", 
    size = 6,
    hide.ns = TRUE   
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),  # 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.title = element_blank()
  )



#### N-T subtype proportion ####
cosmx <- readr::read_rds('merge_data_final.rds.gz')


boundary <- as.data.frame(read.csv('DBSCAN_Boundary.csv', row.names = 1))
all(colnames(cosmx)%in%rownames(boundary))
cosmx$space_area <- boundary[colnames(cosmx), 12]
table(cosmx$space_area)

sub <- cosmx
sub <- subset(sub, subtype != 'EpithelialCells')

ClusterFreq <- sub@meta.data[,c("sub_subtype","sample")] %>% table %>%
  data.frame() %>% set_colnames(c("sub_subtype","sample","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(sub_subtype,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=sub_subtype,value=Per,-sample)
ClusterPer$sample=factor(ClusterPer$sample,levels = names(table(sub$sample)))
ClusterPer <- ClusterPer %>% mutate(group = ifelse(grepl("N", sample), "N", "T"))


sum(ClusterPer$Per[ClusterPer$sample == 'MSI_T5'])

color <- c("N" = "#0076AA", "T" = "#EE2A25")

ggplot(ClusterPer, aes(x = sub_subtype, y = Per, fill = group)) +
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
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.title = element_blank()
  )



#### celltype correlation ####
library(Seurat)
library(tidyverse)

meta_all <- readr::read_rds('MergeMetaDataNewMainTypes20241016.rds.gz')
meta <- meta[!is.na(meta$MMR_Status), ]
meta <- meta[meta$Tissue == 'Tumor', ]


tumor <- readr::read_rds("ref_tumor_mp.rds.gz") %>% subset(MP %in% c('EMT-II','EMT-III')) 

map_df <- unique(meta_all[, c("Sample", "Tissue")])
sample2tissue <- setNames(as.character(map_df$Tissue), as.character(map_df$Sample))
s2_sample <- as.character(tumor@meta.data$Sample)
tumor@meta.data$Tissue <- unname(sample2tissue[s2_sample])

tumor <- subset(tumor, Tissue == 'Tumor')

tumor$DefineTypes <- tumor$MP
table(tumor$DefineTypes)
df <- tumor@meta.data[,c('Patient', 'MMR_Status', 'DefineTypes')]

meta <- meta[meta$Sample%in%tumor$Sample,]
meta <- meta[,c('Patient', 'MMR_Status', 'DefineTypes')]

data <- rbind(meta, df)

meta_data <- data %>%
  dplyr::select(MMR_Status, Patient, DefineTypes) %>%
  mutate(
    Group = as.character(MMR_Status),
    Sample = as.character(Patient),
    CellType = as.character(DefineTypes)
  )
group_a_data <- meta_data %>%
  dplyr::filter(Group == "MSS", CellType %in% c("Glia", "EMT-II")) %>%
  mutate(across(everything(), as.vector)) %>%
  as.data.frame()
cell_count_df <- group_a_data %>%
  dplyr::count(Sample, CellType) %>%
  tidyr::pivot_wider(
    names_from = CellType, 
    values_from = n, 
    values_fill = list(n = 0) 
  )

cor_test_a <- cor.test(cell_count_df$Glia, cell_count_df$`EMT-II`, method = "spearman")
cat("Group A Correlation (Glia vs EMT):\n",
    "R =", round(cor_test_a$estimate, 3), 
    "P =", format.pval(cor_test_a$p.value, digits = 3), "\n")

lm_model <- lm(`EMT-II` ~ Glia, data = cell_count_df)
summary_model <- summary(lm_model)

r_value <- round(cor_test_a$estimate, 3)
p_value <- ifelse(cor_test_a$p.value < 0.001, "< 0.001", 
                  paste0("= ", format(round(cor_test_a$p.value, 3), nsmall = 3)))
equation <- paste0("y = ", round(coef(lm_model)[2], 3), "x + ", round(coef(lm_model)[1], 3))
r_squared <- paste0("R² = ", round(summary_model$r.squared, 3))

stat_label <- paste0(
  "R = ", r_value, "\n",
  "P ", p_value, "\n"
)

p1 <- ggplot(cell_count_df, aes(x = Glia, y = `EMT-II`)) +
  geom_point(size = 1, color = "#0000FF", alpha = 1, stroke = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "#FF0000", fill = "#FFB6C1", linewidth = 0.5) +  
  labs(x = "Glia Count per Sample",
       y = "EMT-II Count per Sample") +
  annotate("text", 
           x = min(cell_count_df$Glia), 
           y = max(cell_count_df$`EMT-II`),
           label = stat_label,
           hjust = 0, vjust = 1, size = 4, color = "black") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

meta_data <- data %>%
  dplyr::select(MMR_Status, Patient, DefineTypes) %>%
  mutate(
    Group = as.character(MMR_Status),
    Sample = as.character(Patient),
    CellType = as.character(DefineTypes)
  )
group_a_data <- meta_data %>%
  dplyr::filter(Group == "MSS", CellType %in% c("Glia", "EMT-III")) %>%
  mutate(across(everything(), as.vector)) %>%
  as.data.frame()
cell_count_df <- group_a_data %>%
  dplyr::count(Sample, CellType) %>%
  tidyr::pivot_wider(
    names_from = CellType, 
    values_from = n, 
    values_fill = list(n = 0) 
  )

cor_test_a <- cor.test(cell_count_df$Glia, cell_count_df$`EMT-III`, method = "spearman")
cat("Group A Correlation (Glia vs EMT):\n",
    "R =", round(cor_test_a$estimate, 3), 
    "P =", format.pval(cor_test_a$p.value, digits = 3), "\n")

lm_model <- lm(`EMT-III` ~ Glia, data = cell_count_df)
summary_model <- summary(lm_model)

r_value <- round(cor_test_a$estimate, 3)
p_value <- ifelse(cor_test_a$p.value < 0.001, "< 0.001", 
                  paste0("= ", format(round(cor_test_a$p.value, 3), nsmall = 3)))
equation <- paste0("y = ", round(coef(lm_model)[2], 3), "x + ", round(coef(lm_model)[1], 3))
r_squared <- paste0("R² = ", round(summary_model$r.squared, 3))

stat_label <- paste0(
  "R = ", r_value, "\n",
  "P ", p_value, "\n"
)

p2 <- ggplot(cell_count_df, aes(x = Glia, y = `EMT-III`)) +
  geom_point(size = 1, color = "#0000FF", alpha = 1, stroke = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "#FF0000", fill = "#FFB6C1", linewidth = 0.5) +  
  labs(x = "Glia Count per Sample",
       y = "EMT-III Count per Sample") +
  annotate("text", 
           x = min(cell_count_df$Glia), 
           y = max(cell_count_df$`EMT-III`),
           label = stat_label,
           hjust = 0, vjust = 1, size = 4, color = "black") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))


meta_data <- data %>%
  dplyr::select(MMR_Status, Patient, DefineTypes) %>%
  mutate(
    Group = as.character(MMR_Status),
    Sample = as.character(Patient),
    CellType = as.character(DefineTypes)
  )
group_a_data <- meta_data %>%
  dplyr::filter(Group == "MSS", CellType %in% c("Tip Cells", "EMT-II")) %>%
  mutate(across(everything(), as.vector)) %>%
  as.data.frame()
cell_count_df <- group_a_data %>%
  dplyr::count(Sample, CellType) %>%
  tidyr::pivot_wider(
    names_from = CellType, 
    values_from = n, 
    values_fill = list(n = 0) 
  )

cor_test_a <- cor.test(cell_count_df$`Tip Cells`, cell_count_df$`EMT-II`, method = "spearman")
cat("Group A Correlation (Glia vs EMT):\n",
    "R =", round(cor_test_a$estimate, 3), 
    "P =", format.pval(cor_test_a$p.value, digits = 3), "\n")

lm_model <- lm(`EMT-II` ~ `Tip Cells`, data = cell_count_df)
summary_model <- summary(lm_model)

r_value <- round(cor_test_a$estimate, 3)
p_value <- ifelse(cor_test_a$p.value < 0.001, "< 0.001", 
                  paste0("= ", format(round(cor_test_a$p.value, 3), nsmall = 3)))
equation <- paste0("y = ", round(coef(lm_model)[2], 3), "x + ", round(coef(lm_model)[1], 3))
r_squared <- paste0("R² = ", round(summary_model$r.squared, 3))

stat_label <- paste0(
  "R = ", r_value, "\n",
  "P ", p_value, "\n"
)

p3 <- ggplot(cell_count_df, aes(x = `Tip Cells`, y = `EMT-II`)) +
  geom_point(size = 1, color = "#0000FF", alpha = 1, stroke = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "#FF0000", fill = "#FFB6C1", linewidth = 0.5) +  
  labs(x = "Tip Cells Count per Sample",
       y = "EMT-II Count per Sample") +
  annotate("text", 
           x = min(cell_count_df$`Tip Cells`), 
           y = max(cell_count_df$`EMT-II`),
           label = stat_label,
           hjust = 0, vjust = 1, size = 4, color = "black") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))


# Tip Cells与 EMT-III
meta_data <- data %>%
  dplyr::select(MMR_Status, Patient, DefineTypes) %>%
  mutate(
    Group = as.character(MMR_Status),
    Sample = as.character(Patient),
    CellType = as.character(DefineTypes)
  )
group_a_data <- meta_data %>%
  dplyr::filter(Group == "MSS", CellType %in% c("Tip Cells", "EMT-III")) %>%
  mutate(across(everything(), as.vector)) %>%
  as.data.frame()
cell_count_df <- group_a_data %>%
  dplyr::count(Sample, CellType) %>%
  tidyr::pivot_wider(
    names_from = CellType, 
    values_from = n, 
    values_fill = list(n = 0) 
  )

cor_test_a <- cor.test(cell_count_df$`Tip Cells`, cell_count_df$`EMT-III`, method = "spearman")
cat("Group A Correlation (Glia vs EMT):\n",
    "R =", round(cor_test_a$estimate, 3), 
    "P =", format.pval(cor_test_a$p.value, digits = 3), "\n")

lm_model <- lm(`EMT-III` ~ `Tip Cells`, data = cell_count_df)
summary_model <- summary(lm_model)

r_value <- round(cor_test_a$estimate, 3)
p_value <- ifelse(cor_test_a$p.value < 0.001, "< 0.001", 
                  paste0("= ", format(round(cor_test_a$p.value, 3), nsmall = 3)))
equation <- paste0("y = ", round(coef(lm_model)[2], 3), "x + ", round(coef(lm_model)[1], 3))
r_squared <- paste0("R² = ", round(summary_model$r.squared, 3))

stat_label <- paste0(
  "R = ", r_value, "\n",
  "P ", p_value, "\n"
)

p4 <- ggplot(cell_count_df, aes(x = `Tip Cells`, y = `EMT-III`)) +
  geom_point(size = 1, color = "#0000FF", alpha = 1, stroke = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "#FF0000", fill = "#FFB6C1", linewidth = 0.5) +  
  labs(x = "Tip Cells Count per Sample",
       y = "EMT-III Count per Sample") +
  annotate("text", 
           x = min(cell_count_df$`Tip Cells`), 
           y = max(cell_count_df$`EMT-III`),
           label = stat_label,
           hjust = 0, vjust = 1, size = 4, color = "black") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))


(p1+p2)/(p3+p4)







