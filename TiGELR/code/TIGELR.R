library(dplyr)  
library(magrittr)
library(rlang)
library(pROC)
library(dplyr)  
library(randomForest)
library(caret)
library(tibble)
library(tidyr)

set.seed(42)

bulk_data_train <- readRDS('pseudobulk_sc_new.rds')

MMR_train <- as.data.frame(bulk_data_train$MMR) %>% 
  set_names('subtype') %>% 
  rownames_to_column()

genes <- readr::read_rds('LR_genes.rds')

data_train <- as.matrix(bulk_data_train@assays$RNA@counts)
bulk_train <- t(data_train)
bulk_train <- as.data.frame(bulk_train[, intersect(genes, colnames(bulk_train))]) %>% 
  rownames_to_column()

merged_data_train <- merge(bulk_train, MMR_train, by = "rowname")

merged_data_valid <- readRDS('MMR_external_bulk.rds')


X_train_raw <- merged_data_train[, intersect(genes, colnames(bulk_train))]
y_train <- as.factor(merged_data_train$subtype)

X_valid_raw <- merged_data_valid[, intersect(genes, colnames(merged_data_valid))]
y_valid <- as.factor(merged_data_valid$subtype)

preProc <- preProcess(X_train_raw, method = c("center", "scale"))
X_train_proc <- predict(preProc, X_train_raw)
X_valid_proc <- predict(preProc, X_valid_raw)

ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  classProbs = TRUE,
  verboseIter = TRUE,
  summaryFunction = twoClassSummary
)

rf_grid <- expand.grid(mtry = seq(1, ncol(X_train_proc), by = 1))
xgb_grid <- expand.grid(
  nrounds = c(50, 100),
  max_depth = c(3, 5),
  eta = c(0.01, 0.05, 0.1),
  gamma = c(0, 0.1),
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.5
)
glmnet_grid <- expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-3, 2, length = 8))
lasso_grid <- expand.grid(alpha = 1, lambda = 10^seq(-3, 2, length = 10))
svm_grid <- expand.grid(C = 10^seq(-2, 2, length = 5))
knn_grid <- expand.grid(k = seq(1, 21, by = 2))
gbm_grid <- expand.grid(
  n.trees = c(100, 200, 500),
  interaction.depth = c(2, 3, 5),
  shrinkage = c(0.01, 0.05, 0.1),
  n.minobsinnode = c(1, 3, 5)
)

models <- list(
  rf = caret::train(X_train_proc, y_train,
                    method = "rf",
                    trControl = ctrl,
                    tuneGrid = rf_grid,
                    metric = "ROC"),
  xgb = caret::train(X_train_proc, y_train,
                     method = "xgbTree",
                     trControl = ctrl,
                     tuneGrid = xgb_grid,
                     metric = "ROC",
                     nthread = 50),
  glm = caret::train(X_train_proc, y_train,
                     method = "glm",
                     trControl = ctrl,
                     metric = "ROC",
                     family = binomial()),
  lasso = caret::train(X_train_proc, y_train,
                       method = "glmnet",
                       trControl = ctrl,
                       tuneGrid = lasso_grid,
                       metric = "ROC"),
  svm = caret::train(X_train_proc, y_train,
                     method = "svmLinear",
                     trControl = ctrl,
                     tuneGrid = svm_grid,
                     metric = "ROC",
                     probability = TRUE), # probability = TRUE
  knn = caret::train(X_train_proc, y_train,
                     method = "knn",
                     trControl = ctrl,
                     tuneGrid = knn_grid,
                     metric = "ROC"),
  gbm = caret::train(X_train_proc, y_train,
                     method = "gbm",
                     trControl = ctrl,
                     tuneGrid = gbm_grid,
                     metric = "ROC",
                     verbose = FALSE)
)

results_list <- list()
sample_pred_list <- list()

for (model_name in names(models)) {
  print(model_name)
  model <- models[[model_name]]

  pred_prob <- predict(model, X_valid_proc, type = "prob")
  pred_class <- predict(model, X_valid_proc)

  cm <- confusionMatrix(pred_class, y_valid, positive = "MSS")
  
  roc_obj <- roc(
    response = y_valid,
    predictor = pred_prob[, "MSS"],
    levels = c("MSI", "MSS"),
    direction = "<"
  )

  metrics <- data.frame(
    Model = model_name,
    ACC = as.numeric(cm$overall["Accuracy"]),
    Sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    Specificity = as.numeric(cm$byClass["Specificity"]),
    F1 = as.numeric(cm$byClass["F1"]),
    AUC = as.numeric(auc(roc_obj)),
    stringsAsFactors = FALSE
  )
  
  results_list[[model_name]] <- metrics

  sample_pred_df <- data.frame(
    Model = model_name,
    Sample = rownames(X_valid_proc),
    TrueLabel = as.character(y_valid),
    PredProb = pred_prob[, "MSS"],
    stringsAsFactors = FALSE
  )
  
  sample_pred_list[[model_name]] <- sample_pred_df
}

final_results <- bind_rows(results_list)
sample_pred_all <- bind_rows(sample_pred_list)

roc_combined <- list()
for (model in names(models)) {
  model_data <- sample_pred_all %>% filter(Model == model)
  roc_combined[[model]] <- roc(
    response = model_data$TrueLabel,
    predictor = model_data$PredProb,
    levels = c("MSI", "MSS"),
    direction = "<"
  )
}

mean_auc <- sapply(roc_combined, auc)

plot_data <- data.frame()
for (model in names(roc_combined)) {
  roc_obj <- roc_combined[[model]]
  model_data <- data.frame(
    Model = model,
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    stringsAsFactors = FALSE
  )
  plot_data <- rbind(plot_data, model_data)
}


my_palette <- c(
  "glm" = "#E58F8E",   
  "xgb" = "#ecb46c",   
  "rf" = "#559d4f",    
  "svm" = "#6A5ACD",   
  "knn" = "#20B2AA",   
  "lasso" = "#9370DB",  
  "gbm" = "#FF6347"    
)


plot.roc(roc_combined$glm, 
         col = my_palette["glm"], 
         ylim = c(0, 1), 
         legacy.axes = TRUE,
         print.auc = TRUE, 
         print.auc.y = 0.7,
         auc.polygon.col = paste0(my_palette["glm"], "20"),
         lwd = 2)
plot.roc(roc_combined$xgb, 
         add = TRUE, 
         col = my_palette["xgb"],
         print.auc = TRUE, 
         print.auc.y = 0.6,
         auc.polygon.col = paste0(my_palette["xgb"], "20"),
         lwd = 2)
plot.roc(roc_combined$rf, 
         add = TRUE, 
         col = my_palette["rf"],
         print.auc = TRUE, 
         print.auc.y = 0.5,
         auc.polygon.col = paste0(my_palette["rf"], "20"),
         lwd = 2)
plot.roc(roc_combined$svm, 
         add = TRUE, 
         col = my_palette["svm"],
         print.auc = TRUE, 
         print.auc.y = 0.4,
         auc.polygon.col = paste0(my_palette["svm"], "20"),
         lwd = 2)
plot.roc(roc_combined$knn, 
         add = TRUE, 
         col = my_palette["knn"],
         print.auc = TRUE, 
         print.auc.y = 0.3,
         auc.polygon.col = paste0(my_palette["knn"], "20"),
         lwd = 2)
plot.roc(roc_combined$lasso, 
         add = TRUE, 
         col = my_palette["lasso"],
         print.auc = TRUE, 
         print.auc.y = 0.2,
         auc.polygon.col = paste0(my_palette["lasso"], "20"),
         lwd = 2)
plot.roc(roc_combined$gbm, 
         add = TRUE, 
         col = my_palette["gbm"],
         print.auc = TRUE, 
         print.auc.y = 0.1,
         auc.polygon.col = paste0(my_palette["gbm"], "20"),
         lwd = 2)
legend("bottomright", 
       legend = c(paste0("glm (AUC = ", round(mean_auc["glm"], 3), ")"),
                  paste0("xgb (AUC = ", round(mean_auc["xgb"], 3), ")"),
                  paste0("rf (AUC = ", round(mean_auc["rf"], 3), ")"),
                  paste0("svm (AUC = ", round(mean_auc["svm"], 3), ")"),
                  paste0("knn (AUC = ", round(mean_auc["knn"], 3), ")"),
                  paste0("lasso (AUC = ", round(mean_auc["lasso"], 3), ")"),
                  paste0("gbm (AUC = ", round(mean_auc["gbm"], 3), ")")),
       col = my_palette, 
       lwd = 2, 
       cex = 0.9,
       bty = "n")



model_colors <- c(
  "glm" = "#E58F8E",  
  "xgb" = "#ecb46c",  
  "rf" = "#559d4f",   
  "svm" = "#6A5ACD",   
  "knn" = "#20B2AA",  
  "lasso" = "#9370DB",   
  "gbm" = "#FF6347"   
)

auc_plot <- ggplot(final_results, aes(x = Model, y = AUC, fill = Model)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2) +
  scale_fill_manual(values = model_colors) +
  labs(title = "Model Performance Across 5 Folds",
       subtitle = "AUC Distribution",
       y = "AUC") +
  theme_minimal()

acc_plot <- ggplot(final_results, aes(x = Model, y = ACC, fill = Model)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2) +
  scale_fill_manual(values = model_colors) +
  labs(y = "Accuracy") +
  theme_minimal()

f1_plot <- ggplot(final_results, aes(x = Model, y = F1, fill = Model)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2) +
  scale_fill_manual(values = model_colors) +
  labs(y = "F1 Score") +
  theme_minimal()
