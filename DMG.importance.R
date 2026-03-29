# ============================================================
# Feature Importance Ranking + Heatmap for 100+ Features
# Methods: Random Forest, XGBoost, Lasso, Pearson Correlation
# ============================================================

# --- 1. Install & Load Libraries ---
packages <- c("randomForest", "xgboost", "glmnet", "ggplot2",
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
y    <- as.factor(df[[args[2]]])
X    <- df %>% select(-all_of(args[2]))

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
# 6. METHOD C — Lasso |Coefficient|
# ============================================================
cat("Training Lasso...\n")
lasso_cv   <- cv.glmnet(x = as.matrix(X_train), y = y_train,
                        alpha = 1, family = "binomial", nfolds = 5)
lasso_coef <- coef(lasso_cv, s = "lambda.min")
lasso_df   <- data.frame(
  Feature        = rownames(lasso_coef)[-1],
  Lasso_Coef_Abs = abs(as.numeric(lasso_coef[-1, 1]))
)

# ============================================================
# 7. METHOD D — Pearson |Correlation| with Target
# ============================================================
cat("Computing correlations...\n")
y_num    <- as.numeric(y_train) - 1
corr_vec <- sapply(X_train, function(col) abs(cor(col, y_num)))
corr_df  <- data.frame(Feature     = names(corr_vec),
                       Correlation = as.numeric(corr_vec))

# ============================================================
# 8. MERGE + NORMALISE + COMPOSITE SCORE
# ============================================================
normalize <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

importance_df <- rf_df %>%
  full_join(xgb_df,   by = "Feature") %>%
  full_join(lasso_df, by = "Feature") %>%
  full_join(corr_df,  by = "Feature") %>%
  replace(is.na(.), 0) %>%
  mutate(across(c(RF_Importance, XGB_Gain, Lasso_Coef_Abs, Correlation),
                normalize, .names = "Norm_{.col}")) %>%
  rowwise() %>%
  mutate(Composite_Score = mean(c(Norm_RF_Importance, Norm_XGB_Gain,
                                  Norm_Lasso_Coef_Abs, Norm_Correlation))) %>%
  ungroup() %>%
  arrange(desc(Composite_Score))

cat("\nTop 15 features:\n")
print(importance_df %>% select(Feature, Composite_Score) %>% head(15))

# ============================================================
# 9. HEATMAP — Top 50 Features × 4 Methods
# ============================================================
top_n        <- 50
top_features <- importance_df$Feature[1:top_n]

heatmap_long <- importance_df %>%
  filter(Feature %in% top_features) %>%
  select(Feature,
         `Random Forest` = Norm_RF_Importance,
         `XGBoost`       = Norm_XGB_Gain,
         `Lasso`         = Norm_Lasso_Coef_Abs,
         `Correlation`   = Norm_Correlation) %>%
  mutate(Feature = factor(Feature, levels = rev(top_features))) %>%
  pivot_longer(-Feature, names_to = "Method", values_to = "Score")

p_heatmap <- ggplot(heatmap_long,
                    aes(x = Method, y = Feature, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(Score > 0.05, round(Score, 2), "")),
            size = 2.4, color = "white", fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#0d0221", "#1b3a6b", "#1a6eb5",
               "#48b3e8", "#f5dd90", "#f27059"),
    values = c(0, 0.15, 0.35, 0.55, 0.75, 1),
    limits = c(0, 1),
    name   = "Normalised\nImportance"
  ) +
  labs(
    title    = paste0("Feature Importance Heatmap  |  Top ",
                      top_n, " of ", n_features, " Features"),
    subtitle = "Scores normalised [0–1] per method · Warmer colour = More Important",
    x        = NULL, y = "Feature",
    caption  = paste0("RF (MeanDecreaseGini) · XGBoost (Gain) · ",
                      "Lasso (|Coef|) · |Pearson r|")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(face = "bold", size = 14),
    plot.subtitle     = element_text(size = 9, color = "grey45"),
    plot.caption      = element_text(size = 7, color = "grey55"),
    axis.text.x       = element_text(face = "bold", size = 11),
    axis.text.y       = element_text(size = 7.5),
    axis.ticks        = element_blank(),
    panel.grid        = element_blank(),
    legend.key.height = unit(1.8, "cm"),
    plot.background   = element_rect(fill = "white", color = NA)
  )

# ============================================================
# 10. COMPOSITE BAR CHART — Top 30 Features
# ============================================================
top30 <- importance_df %>%
  slice_head(n = 30) %>%
  mutate(Feature = factor(Feature, levels = rev(Feature)),
         Rank    = row_number())

p_bar <- ggplot(top30, aes(x = Composite_Score, y = Feature,
                           fill = Composite_Score)) +
  geom_col(show.legend = FALSE, width = 0.72) +
  geom_text(aes(label = paste0("#", Rank, "  ", round(Composite_Score, 3))),
            hjust = -0.05, size = 3, color = "grey30") +
  scale_fill_gradient(low = "#48b3e8", high = "#f27059") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title    = "Top 30 Features — Composite Importance Score",
       subtitle = "Equal-weight average of normalised RF · XGB · Lasso · Correlation",
       x = "Composite Score (0–1)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.background    = element_rect(fill = "white", color = NA))

# ============================================================
# 11. SAVE
# ============================================================
ggsave("feature_importance_heatmap.png",
       plot = p_heatmap, width = 11, height = 16, dpi = 180, bg = "white")
ggsave("feature_importance_barplot.png",
       plot = p_bar, width = 10, height = 9, dpi = 180, bg = "white")
write.csv(importance_df, "feature_importance_scores.csv", row.names = FALSE)

cat("\nDone! Outputs:\n",
    "  feature_importance_heatmap.png\n",
    "  feature_importance_barplot.png\n",
    "  feature_importance_scores.csv\n")

#execution
#bashRscript feature_importance_heatmap.R my_data.csv target_column


                   
