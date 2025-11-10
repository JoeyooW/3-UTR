# =========== Final join: per-variant mean scores -> map back to Mcont IDs ===========
if (!requireNamespace("duckdb", quietly = TRUE)) install.packages("duckdb")
if (!requireNamespace("DBI", quietly = TRUE)) install.packages("DBI")
library(DBI); library(duckdb)

# 1) 文件路径（改成你的）
# big_output: 你要汇总/对齐的结果文件（可以是你前面导出的 filtered_output.tsv，或更早的原始CSV）
#   - 若是你刚导出的聚合表，用 '\t' 分隔；若是原始 AlphaGenome 导出的文件，多半是逗号CSV
big_output   <- "/Users/JOEY/Desktop/UQ/Research/data/Alphagenome/filtered_output.tsv"  # filtered result
input_file <- "/Users/JOEY/Desktop/UQ/Research/data/Alphagenome/TEST_VCF_Both_alphagenome.txt" # input for alphagenome
final_file <- "/Users/JOEY/Desktop/UQ/Research/data/Alphagenome/final_variant_scores.tsv"

big_is_tsv   <- TRUE                               # TRUE=制表符；FALSE=逗号

db_file      <- "alphagenome.duckdb"

delim_big <- if (big_is_tsv) "\\t" else ","  # 给 DuckDB 的分隔符

con <- dbConnect(duckdb(), dbdir = db_file)

# 2) 直接在 DuckDB 里完成：解析 variant_id -> CHROM,POS,REF,ALT；聚合；与输入左连接；导出
sql <- sprintf("
CREATE OR REPLACE TEMP VIEW out_raw AS
SELECT
  regexp_replace(split_part(variant_id, ':', 1), '^chr', '')      AS CHROM_key,
  TRY_CAST(split_part(variant_id, ':', 2) AS BIGINT)              AS POS,
  UPPER(TRIM(split_part(split_part(variant_id, ':', 3), '>', 1))) AS REF,
  UPPER(TRIM(split_part(split_part(variant_id, ':', 3), '>', 2))) AS ALT,
  CAST(mean_raw_score AS DOUBLE)        AS raw_val,
  CAST(mean_quantile_score AS DOUBLE)   AS quant_val
FROM read_csv('%s', delim='\\t', header=TRUE);

CREATE OR REPLACE TEMP VIEW out_agg AS
SELECT
  CHROM_key, POS, REF, ALT,
  -- raw
  AVG(raw_val)                    AS mean_raw,
  STDDEV(raw_val)                 AS sd_raw,
  MAX(raw_val) - MIN(raw_val)     AS range_raw,
  MAX(raw_val)                    AS max_raw,
  -- quantile
  AVG(quant_val)                  AS mean_quant,
  STDDEV(quant_val)               AS sd_quant,
  MAX(quant_val) - MIN(quant_val) AS range_quant,
  MAX(quant_val)                  AS max_quant
FROM out_raw
GROUP BY CHROM_key, POS, REF, ALT;

COPY (
  SELECT
    i.variant_id         AS Mcont_id,
    CAST(i.CHROM AS VARCHAR) AS CHROM_key,
    TRY_CAST(i.POS AS BIGINT)  AS POS,
    UPPER(TRIM(i.REF))   AS REF,
    UPPER(TRIM(i.ALT))   AS ALT,
    o.mean_raw,
    o.sd_raw,
    o.range_raw,
    o.max_raw,
    o.mean_quant,
    o.sd_quant,
    o.range_quant,
    o.max_quant
  FROM read_csv('%s', delim='\\t', header=TRUE) AS i
  LEFT JOIN out_agg AS o
    ON CAST(i.CHROM AS VARCHAR) = o.CHROM_key
   AND TRY_CAST(i.POS AS BIGINT) = o.POS
   AND UPPER(TRIM(i.REF)) = o.REF
   AND UPPER(TRIM(i.ALT)) = o.ALT
) TO '%s' WITH (HEADER, DELIMITER '\\t');
", big_output, input_file, final_file)


dbExecute(con, sql)

# 3) 可选：给点统计反馈（匹配了多少）
stats <- dbGetQuery(con, sprintf("
  WITH T AS (
    SELECT *
    FROM read_csv('%s', delim='\\t', header=TRUE)
  )
  SELECT
    COUNT(*)                                   AS n_input,
    SUM(CASE WHEN mean_raw   IS NOT NULL OR mean_quant IS NOT NULL THEN 1 ELSE 0 END) AS n_matched,
    SUM(CASE WHEN mean_raw   IS NULL AND mean_quant IS NULL THEN 1 ELSE 0 END)        AS n_unmatched
  FROM T
", final_file))

dbDisconnect(con, shutdown = TRUE)

cat('✅ 完成！结果已写入:\n', final_file, '\n')
print(stats)




# ---------manually put the info in------

final_file <- read.csv("/Users/JOEY/Desktop/UQ/Research/data/Alphagenome/final_variant_scores.tsv", sep = "\t")
head(final_file,3)

library(ggplot2)
# mean_raw
ggplot(final_file, aes(x = INFO, y = log2(mean_raw), fill = INFO)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots of Alphagenome score mean", x = "Group", y = "log2(Raw_Score)") +
  scale_fill_manual(values = c("Cont" = "lightblue", "Dis" = "salmon"))


# sd_raw
ggplot(final_file, aes(x = INFO, y = log2(sd_raw), fill = INFO)) +
geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots of Alphagenome score sd", x = "Group", y = "log2(sd_raw)") +
  scale_fill_manual(values = c("Cont" = "lightblue", "Dis" = "salmon"))


# range_raw
ggplot(final_file, aes(x = INFO, y = log2(range_raw), fill = INFO)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots of Alphagenome score range", x = "Group", y = " log2(range_raw)") +
  scale_fill_manual(values = c("Cont" = "lightblue", "Dis" = "salmon"))


# max_raw
ggplot(final_file, aes(x = INFO, y = log2(max_raw), fill = INFO)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots of Alphagenome", x = "Group", y = "log2(max_raw)") +
  scale_fill_manual(values = c("Cont" = "lightblue", "Dis" = "salmon"))


# mean_quant
ggplot(final_file, aes(x = INFO, y = mean_quant, fill = INFO)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots of Alphagenome", x = "Group", y = "Quantile_Score") +
  scale_fill_manual(values = c("Cont" = "lightblue", "Dis" = "salmon"))



# sd_quant
ggplot(final_file, aes(x = INFO, y = sd_quant, fill = INFO)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots of Alphagenome", x = "Group", y = "sd_quant") +
  scale_fill_manual(values = c("Cont" = "lightblue", "Dis" = "salmon"))

# range_quant
ggplot(final_file, aes(x = INFO, y = range_quant, fill = INFO)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots of Alphagenome", x = "Group", y = "range_quant") +
  scale_fill_manual(values = c("Cont" = "lightblue", "Dis" = "salmon"))


# max_quant
ggplot(final_file, aes(x = INFO, y = max_quant, fill = INFO)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots of Alphagenome", x = "Group", y = "max_quant") +
  scale_fill_manual(values = c("Cont" = "lightblue", "Dis" = "salmon"))

#---------- density plot-------



library(ggplot2)
library(dplyr)

# 定义配色方案
colors <- c("Cont" = "cyan", "Dis" = "red")

# 通用绘图函数：绘图 + 输出均值
plot_density <- function(data, variable, log_transform = FALSE, xlab = NULL, title = NULL) {
  
  df <- data %>% select(INFO, all_of(variable)) %>%
    rename(value = all_of(variable))
  
  # log2 转换（可选）
  if (log_transform) {
    df <- df %>% mutate(value = log2(value))
    if (is.null(xlab)) xlab <- paste0("log2(", variable, ")")
  } else {
    if (is.null(xlab)) xlab <- variable
  }
  
  # 计算均值并输出
  means <- df %>%
    group_by(INFO) %>%
    summarise(mean_value = mean(value, na.rm = TRUE),
              sd_value = sd(value, na.rm = TRUE),
              n = n())
  
  cat("\n------------------------------\n")
  cat("Variable:", variable, "\n")
  print(means)
  cat("------------------------------\n")
  
  # 绘制 density 图
  ggplot(df, aes(x = value, fill = INFO)) +
    geom_density(color = "black", alpha = 0.4, adjust = 1.2) +
    geom_vline(data = means, aes(xintercept = mean_value),
               color = "black", linetype = "dashed", size = 0.4) +
    theme_minimal() +
    labs(title = title, x = xlab, y = "Density") +
    theme(text = element_text(size = 12),
          legend.position = "none")
}

# ====== 绘图并打印均值 ======
plot_density(final_file, "mean_raw",  TRUE, "log2(mean_raw)",  "AlphaGenome mean raw score distribution")
plot_density(final_file, "sd_raw",    TRUE, "log2(sd_raw)",    "AlphaGenome sd raw score distribution")
plot_density(final_file, "range_raw", TRUE, "log2(range_raw)", "AlphaGenome range raw score distribution")
plot_density(final_file, "max_raw",   TRUE, "log2(max_raw)",   "AlphaGenome max raw score distribution")

plot_density(final_file, "mean_quant",  FALSE, "mean_quant",  "AlphaGenome mean quantile score distribution")
plot_density(final_file, "sd_quant",    FALSE, "sd_quant",    "AlphaGenome sd quantile score distribution")
plot_density(final_file, "range_quant", FALSE, "range_quant", "AlphaGenome range quantile score distribution")
plot_density(final_file, "max_quant",   FALSE, "max_quant",   "AlphaGenome max quantile score distribution")


