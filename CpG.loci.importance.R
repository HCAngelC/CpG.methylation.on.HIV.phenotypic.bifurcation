# ============================================================
# Loci Importance Ranking + Heatmap for methylation CpG sites and DMGs
# Methods: Random Forest, XGBoost, Lasso, Pearson Correlation
# ============================================================

# --- 1. Install & Load Libraries ---
packages <- c("randomForest", "xgboost", "ggplot2",
              "dplyr", "tidyr", "caret")

installed <- rownames(installed.packages())
to_install <- packages[!packages %in% installed]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
lapply(packages, library, character.only = TRUE)

# ============================================================
# 2. Data loading
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
df   <- read.csv(args[1])
y    <- as.factor(df[[args[2]]]) #predict column
X    <- df %>% select(-all_of(args[2])) #feature column

# ============================================================
# 3. TRAIN / TEST SPLIT
# ============================================================
train_idx <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_idx, ];  y_train <- y[train_idx]
X_test  <- X[-train_idx, ]; y_test  <- y[-train_idx]

# ============================================================
# 4. METHOD A — Random Forest (MeanDecreaseGini)
# ============================================================
cat("Training Random Forest...\n")
rf_model <- randomForest(x = X_train, y = y_train,
                         ntree = 300, importance = TRUE)
rf_imp <- importance(rf_model, type = 2)[, 1]
rf_df  <- data.frame(Feature       = names(rf_imp),
                     RF_Importance = as.numeric(rf_imp))

# ============================================================
# 5. METHOD B — XGBoost (Gain)
# ============================================================
cat("Training XGBoost...\n")
dtrain    <- xgb.DMatrix(data  = as.matrix(X_train),
                         label = as.numeric(as.character(y_train))
                         - 1)
xgb_model <- xgboost(data = dtrain, nrounds = 150,
                     objective = "binary:logistic",
                     eval_metric = "logloss",
                     eta = 0.1, max_depth = 5, verbose = 0)
xgb_imp <- xgb.importance(model = xgb_model)
xgb_df  <- data.frame(Feature  = xgb_imp$Feature,
                      XGB_Gain = xgb_imp$Gain)

# ============================================================
# 6. METHOD C — Pearson |Correlation| with Target
# ============================================================
cat("Computing correlations...\n")
y_num    <- as.numeric(y_train) - 1
corr_vec <- sapply(X_train, function(col) abs(cor(col, y_num)))
corr_df  <- data.frame(Feature     = names(corr_vec),
                       Correlation = as.numeric(corr_vec))

# ============================================================
# 7. MERGE + NORMALISE + COMPOSITE SCORE
# ============================================================
normalize <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

importance_df <- rf_df %>%
  full_join(xgb_df,   by = "Feature") %>%
  full_join(corr_df,  by = "Feature") %>%
  replace(is.na(.), 0) %>%
  mutate(across(c(RF_Importance, XGB_Gain, Correlation),
                normalize, .names = "Norm_{.col}")) %>%
  rowwise() %>%
  mutate(Competence_Score = mean(c(Norm_RF_Importance, Norm_XGB_Gain,
                                Norm_Correlation))) %>%
  ungroup() %>%
  arrange(desc(Competence_Score))

# ============================================================
# 8. COMPETENCE SCORE BAR CHART
# ============================================================
top <- importance_df %>%
  slice_head(n = 66) %>%
  mutate(Feature = factor(Feature, levels = rev(Feature)),
         Rank    = row_number())

p_bar <- ggplot(top, aes(x = Competence_Score, y = Feature,
                           fill = Competence_Score)) +
  geom_col(show.legend = FALSE, width = 0.72) +
  geom_text(aes(label = paste0("#", Rank, "  ", round(Competence_Score, 3))),
            hjust = -0.05, size = 3, color = "grey30") +
  scale_fill_gradient(low = "#48b3e8", high = "#f27059") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Competence Score (0–1)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.background    = element_rect(fill = "white", color = NA))

# ============================================================
# 11. SAVE
# ============================================================
ggsave("competence_score_barplot.png",
       plot = p_bar, width = 10, height = 9, dpi = 180, bg = "white")
write.csv(importance_df, "importance_scores.csv", row.names = FALSE)

cat("\nDone! Outputs:\n",
    "  competence_score_barplot.png\n",
    "  importance_scores.csv\n")

#execution
#bashRscript CpG.loci.importance.R my_data.csv target_column


                   
