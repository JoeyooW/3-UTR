#deepsea boxplot

# setting
library(ggplot2)
library(dplyr)
library(readr)

# read in
df <- read_tsv("/Users/JOEY/Desktop/UQ/Research/data/DeepSea/DeapSea_VARIANT_dis.tsv")

table(df$Group)

# 画 Boxplot
ggplot(df, aes(x = Group, y = max_score, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  theme_minimal() +
  labs(title = "Disease Impact Score by Variant Group",
       x = "Variant Group", y = "DeepSEA Max Score") +
  theme(legend.position = "none")

# 做统计检验（t-test）
t.test(max_score ~ Group, data = df)


ggplot(df, aes(x = max_score, fill = Group)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Score Distribution by Group", x = "Max Score", y = "Density")


#AUC
library(pROC)
df$label <- ifelse(df$Group == "Dis", 1, 0)
roc_obj <- roc(df$label, df$max_score)
plot(roc_obj, main = "ROC Curve: DeepSEA Max Score")
auc(roc_obj)
