# ============================================================
# Single-file RegVar analysis (INFO-based) -- FINAL (robust Figure1)
# Input: processed_TEST_VCF_Both.csv
# Outputs: enrichment, hit counts, density, scatter matrix,
#          pred_eqtl/gwas by group, motif_RBPs & Rep_info items
#          + auto-save PDF & PNG + multi-panel Figure 1 (a–f)
#          (NO LR calculation; NO magick/pdftools needed)
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(janitor); library(ggrepel); library(GGally)
})

# ---------- paths (EDIT ME) ----------
path_csv <- "/Users/JOEY/Desktop/UQ/Research/data/Romo_Regvar/Romo_Regvar_output/processed_TEST_VCF_Both.csv"
outdir   <- "/Users/JOEY/Desktop/UQ/Research/data/Romo_Regvar/regvar_analysis_final"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------- helpers ----------
HaldaneOR <- function(a,b,c,d) ((a+0.5)/(c+0.5))/((b+0.5)/(d+0.5))
lab_info  <- function(x) factor(x, levels=c(0,1), labels=c("Benign","Pathogenic"))
safe_num  <- function(x) suppressWarnings(as.numeric(x))

# 同步保存 PDF 与 PNG（无需 magick/pdftools）
ggsave_both <- function(file_pdf, plot, width, height, dpi = 300) {
  ggsave(file_pdf, plot = plot, width = width, height = height)
  file_png <- sub("\\.pdf$", ".png", file_pdf, ignore.case = TRUE)
  ggsave(file_png, plot = plot, width = width, height = height, dpi = dpi)
}

# ---------- read ----------
df <- read.csv(path_csv, check.names = FALSE) %>% clean_names()
stopifnot("info" %in% names(df) || "INFO" %in% names(df))
if (!"INFO" %in% names(df)) df <- df %>% mutate(INFO = safe_num(info)) else df$INFO <- safe_num(df$INFO)
stopifnot("var_id" %in% names(df))
tot_path <- sum(df$INFO == 1, na.rm = TRUE)
tot_ben  <- sum(df$INFO == 0, na.rm = TRUE)

# ============================================================
# A) Enrichment: motif / state / motif×state / eCLIP
# ============================================================

# --- A1 motif (motif_RBPs by "_") ---
if ("motif_RBPs" %in% names(df)) {
  motif_long <- df %>%
    filter(!is.na(motif_RBPs)) %>%
    separate_rows(motif_RBPs, sep = "_") %>%
    transmute(var_id, INFO, motif = str_trim(motif_RBPs)) %>%
    filter(motif != "") %>%
    distinct(var_id, motif, .keep_all = TRUE)
  
  motif_counts <- motif_long %>%
    count(motif, INFO, name = "n") %>%
    pivot_wider(names_from = INFO, values_from = n, values_fill = 0) %>%
    rename(Benign = `0`, Pathogenic = `1`)
  
  motif_stats <- motif_counts %>%
    rowwise() %>%
    mutate(
      a = Pathogenic, b = Benign, c = tot_path - a, d = tot_ben - b,
      fisher_p = fisher.test(matrix(c(a,b,c,d), 2, byrow=TRUE))$p.value,
      OR = HaldaneOR(a,b,c,d), log2OR = log2(OR), negLog10P = -log10(fisher_p)
    ) %>% ungroup() %>%
    mutate(bonferroni = p.adjust(fisher_p, "bonferroni"),
           fdr        = p.adjust(fisher_p, "BH")) %>%
    arrange(fdr, fisher_p)
  
  write_csv(motif_stats, file.path(outdir, "A1_motif_stats.csv"))
  
  p_motif <- ggplot(motif_stats, aes(log2OR, negLog10P, label = motif)) +
    geom_point() + geom_vline(xintercept = 0, linetype = 2) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 25) +
    labs(x="log2(OR)", y="-log10(p-value)", title="Motif enrichment (INFO)") +
    theme_minimal()
  ggsave_both(file.path(outdir, "A1_motif_volcano.pdf"), p_motif, 8, 5)
} else {
  motif_long <- tibble(var_id = character(), INFO = numeric(), motif = character())
  message("No motif_RBPs column, skip motif enrichment.")
}

# --- A2 state & motif×state (motif_cat by "__") ---
if (all(c("motif_RBPs","motif_cat") %in% names(df))) {
  rv_pairs <- df %>%
    mutate(
      rbps_list  = ifelse(is.na(motif_RBPs), NA, str_split(motif_RBPs, "_")),
      state_list = ifelse(is.na(motif_cat),   NA, str_split(motif_cat,   "__"))
    ) %>%
    transmute(var_id, INFO, pairs = purrr::map2(rbps_list, state_list, ~{
      x <- .x; y <- .y
      lx <- if (is.null(x)) 0 else length(x)
      ly <- if (is.null(y)) 0 else length(y)
      L  <- max(lx, ly)
      if (L==0) return(tibble(motif=character(), state=character()))
      if (lx<L) x <- c(x, rep(NA_character_, L-lx))
      if (ly<L) y <- c(y, rep(NA_character_, L-ly))
      tibble(motif=x, state=y)
    })) %>%
    unnest(pairs) %>%
    mutate(motif = na_if(str_trim(motif), ""), state = na_if(str_trim(state), "")) %>%
    filter(!is.na(motif), !is.na(state)) %>%
    distinct(var_id, motif, state, .keep_all = TRUE)
  
  # state composition
  rv_state <- rv_pairs %>% distinct(var_id, state, .keep_all = TRUE)
  p_state <- rv_state %>%
    mutate(INFO_f = lab_info(INFO)) %>%
    ggplot(aes(state, fill = INFO_f)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    labs(x="State", y="Proportion", fill="Class", title="Motif state proportion (INFO)") +
    theme_minimal()
  ggsave_both(file.path(outdir, "A2_state_proportion.pdf"), p_state, 8, 5)
  
  # motif×state enrichment
  mxs_counts <- rv_pairs %>%
    unite("motif_state", motif, state, sep = "::", remove = FALSE) %>%
    count(motif_state, INFO, name = "n") %>%
    pivot_wider(names_from = INFO, values_from = n, values_fill = 0) %>%
    rename(Benign = `0`, Pathogenic = `1`)
  
  mxs_stats <- mxs_counts %>%
    rowwise() %>%
    mutate(
      a = Pathogenic, b = Benign, c = tot_path - a, d = tot_ben - b,
      fisher_p = fisher.test(matrix(c(a,b,c,d), 2, byrow=TRUE))$p.value,
      OR = HaldaneOR(a,b,c,d), log2OR = log2(OR), negLog10P = -log10(fisher_p)
    ) %>% ungroup() %>%
    mutate(bonferroni = p.adjust(fisher_p, "bonferroni"),
           fdr        = p.adjust(fisher_p, "BH")) %>%
    arrange(fdr, fisher_p)
  
  write_csv(mxs_stats, file.path(outdir, "A3_motifXstate_stats.csv"))
  
  p_mxs <- ggplot(mxs_stats, aes(log2OR, negLog10P, label = motif_state)) +
    geom_point() + geom_vline(xintercept = 0, linetype = 2) +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 35) +
    labs(x="log2(OR)", y="-log10(p-value)", title="Motif × State enrichment (INFO)") +
    theme_minimal()
  ggsave_both(file.path(outdir, "A3_motifXstate_volcano.pdf"), p_mxs, 9, 6)
} else {
  rv_pairs <- tibble(var_id=character(), INFO=numeric(), motif=character(), state=character())
  message("Missing motif_RBPs or motif_cat, skip state & motif×state.")
}

# --- A3 eCLIP (eclip_tot by "__") ---
if ("eclip_tot" %in% names(df)) {
  eclip_long <- df %>%
    filter(!is.na(eclip_tot)) %>%
    separate_rows(eclip_tot, sep = "__") %>%
    transmute(var_id, INFO, rbp = str_trim(eclip_tot)) %>%
    filter(rbp != "") %>%
    distinct(var_id, rbp, .keep_all = TRUE)
  
  eclip_counts <- eclip_long %>%
    count(rbp, INFO, name = "n") %>%
    pivot_wider(names_from = INFO, values_from = n, values_fill = 0) %>%
    rename(Benign = `0`, Pathogenic = `1`)
  
  eclip_stats <- eclip_counts %>%
    rowwise() %>%
    mutate(
      a = Pathogenic, b = Benign, c = tot_path - a, d = tot_ben - b,
      fisher_p = fisher.test(matrix(c(a,b,c,d), 2, byrow=TRUE))$p.value,
      OR = HaldaneOR(a,b,c,d), log2OR = log2(OR), negLog10P = -log10(fisher_p)
    ) %>% ungroup() %>%
    mutate(bonferroni = p.adjust(fisher_p, "bonferroni"),
           fdr        = p.adjust(fisher_p, "BH")) %>%
    arrange(fdr, fisher_p)
  
  write_csv(eclip_stats, file.path(outdir, "A4_eCLIP_stats.csv"))
  
  p_eclip <- ggplot(eclip_stats, aes(log2OR, negLog10P, label = rbp)) +
    geom_point() + geom_vline(xintercept = 0, linetype = 2) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 25) +
    labs(x="log2(OR)", y="-log10(fisher p)", title="eCLIP enrichment (INFO)") +
    theme_minimal()
  ggsave_both(file.path(outdir, "A4_eCLIP_volcano.pdf"), p_eclip, 8, 5)
} else {
  eclip_long <- tibble(var_id=character(), INFO=numeric(), rbp=character())
  message("No eclip_tot column, skip eCLIP.")
}

# ============================================================
# B) Hits per variant: counts table + 0/1/≥2 bars
# ============================================================
hits_motif <- motif_long %>% count(var_id, name = "n_motif")
hits_mxs   <- rv_pairs   %>% count(var_id, name = "n_mxs")
hits_eclip <- eclip_long %>% count(var_id, name = "n_eclip")

hit_sum <- df %>%
  distinct(var_id, INFO) %>%
  left_join(hits_motif, by="var_id") %>%
  left_join(hits_mxs,   by="var_id") %>%
  left_join(hits_eclip, by="var_id") %>%
  mutate(across(c(n_motif, n_mxs, n_eclip), ~replace_na(., 0L)),
         any_hits   = (n_motif>0) | (n_mxs>0) | (n_eclip>0),
         total_hits = n_motif + n_mxs + n_eclip,
         n_hit_types= (n_motif>0) + (n_mxs>0) + (n_eclip>0))

write_csv(hit_sum, file.path(outdir, "B_hits_per_variant.csv"))

# —— 稳健版：返回图对象，避免白图
plot_one_vs_multi <- function(df, count_col, title_prefix, out_name) {
  df2 <- df %>%
    transmute(INFO_f = lab_info(INFO),
              bucket = case_when(
                .data[[count_col]] == 0 ~ "0",
                .data[[count_col]] == 1 ~ "1",
                .data[[count_col]] >= 2 ~ "≥2")) %>%
    filter(!is.na(bucket)) %>%
    count(INFO_f, bucket, name="n") %>%
    group_by(INFO_f) %>% mutate(prop = n/sum(n)) %>% ungroup()
  
  # 保护：数据不足直接跳过
  if (nrow(df2) == 0 || length(unique(df2$INFO_f)) < 2) {
    message("Skip ", title_prefix, ": not enough data for both groups.")
    return(invisible(NULL))
  }
  
  # 固定桶顺序
  df2$bucket <- factor(df2$bucket, levels = c("0","1","≥2"))
  
  # 全局检验（用 xtabs 避免 rowname 冲突）
  mat <- tryCatch(xtabs(n ~ INFO_f + bucket, data = df2), error = function(e) NULL)
  test_res <- NULL
  if (!is.null(mat) && all(dim(mat) >= 2)) {
    test_res <- tryCatch({
      expv <- suppressWarnings(chisq.test(mat)$expected)
      if (all(expv >= 5)) chisq.test(mat) else fisher.test(mat)
    }, error = function(e) NULL)
  }
  
  p <- ggplot(df2, aes(INFO_f, prop, fill = bucket)) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    labs(x="Class", y="Proportion", fill="Hits",
         title=paste0(title_prefix, " hits per variant",
                      if (!is.null(test_res)) paste0(" (P=", signif(test_res$p.value,3), ")") else "")) +
    theme_minimal()
  
  ggsave_both(file.path(outdir, out_name), p, 7, 5)
  return(p)
}

p_b1 <- plot_one_vs_multi(hit_sum, "n_motif",    "Motif",       "B1_hits_motif_one_vs_multi.pdf")
p_b2 <- plot_one_vs_multi(hit_sum, "n_mxs",      "Motif×State", "B2_hits_mxs_one_vs_multi.pdf")
p_b3 <- plot_one_vs_multi(hit_sum, "n_eclip",    "eCLIP",       "B3_hits_eclip_one_vs_multi.pdf")
p_b4 <- plot_one_vs_multi(hit_sum, "n_hit_types","ANY types",   "B4_hits_any_types_one_vs_multi.pdf")
p_b5 <- plot_one_vs_multi(hit_sum, "total_hits", "ANY total",   "B5_hits_any_total_one_vs_multi.pdf")

# ============================================================
# C) Hit Count vs Pathogenicity: density
# ============================================================
for (col in c("n_motif","n_mxs","n_eclip","total_hits")) {
  dfp <- hit_sum %>%
    mutate(INFO_f = lab_info(INFO)) %>%
    filter(!is.na(.data[[col]]))
  
  p <- ggplot(dfp, aes(.data[[col]], fill = INFO_f)) +
    geom_density(alpha = 0.4) +
    theme_minimal() +
    labs(title = paste0("Distribution of ", col, " by Group"),
         x = paste0(col, " (hit count)"), y = "Density", fill = "Group") +
    theme(legend.position = "top")
  ggsave_both(file.path(outdir, paste0("C_density_", col, ".pdf")), p, 7, 5)
  
  if (col == "n_eclip") p_c_eclip <<- p  # 供 Figure1 使用
}

# ============================================================
# D) Scatter matrix (features inter-relationship)
#     Includes pred_eqtl / pred_gwas if present in CSV
# ============================================================
scatter_df <- hit_sum %>%
  select(var_id, INFO, n_motif, n_mxs, n_eclip, total_hits)

extra_cols <- intersect(c("pred_eqtl","pred_gwas"), names(df))
if (length(extra_cols) > 0) {
  df_num <- df %>% select(var_id, INFO, all_of(extra_cols))
  scatter_df <- scatter_df %>% left_join(df_num, by=c("var_id","INFO"))
}

scatter_df <- scatter_df %>% mutate(INFO_f = lab_info(INFO))
scatter_cols <- c("n_motif", "n_mxs", "n_eclip", "total_hits", extra_cols)
scatter_cols <- scatter_cols[scatter_cols %in% names(scatter_df)]
if (length(scatter_cols) >= 2) {
  p_matrix <- GGally::ggpairs(
    scatter_df,
    columns = scatter_cols,
    mapping = aes(color = INFO_f, alpha = 0.6),
    upper = list(continuous = wrap("cor", size = 3)),
    lower = list(continuous = wrap("points", alpha = 0.3, size = 0.7)),
    diag  = list(continuous = wrap("densityDiag", alpha = 0.4))
  ) + theme_bw(base_size = 9) + theme(legend.position = "bottom")
  ggsave_both(file.path(outdir, "D_scatter_matrix.pdf"), p_matrix, 10, 10)
}

# ============================================================
# F) pred_eqtl / pred_gwas 的组间分析（自动识别二元/连续）
# ============================================================
analyze_pred_col <- function(df, col, out_prefix) {
  if (!col %in% names(df)) { message("Skip ", col, ": not found."); return(invisible(NULL)) }
  dd <- df %>% select(var_id, INFO, !!sym(col)) %>% filter(!is.na(.data[[col]]), !is.na(INFO))
  vals <- unique(na.omit(dd[[col]]))
  is_binary <- length(setdiff(vals, c(0,1))) == 0
  
  if (is_binary) {
    tab <- dd %>%
      mutate(pred = factor(.data[[col]], levels=c(0,1), labels=c("0","1")),
             INFO_f = lab_info(INFO)) %>%
      count(INFO_f, pred, name="n") %>%
      group_by(INFO_f) %>% mutate(prop = n/sum(n)) %>% ungroup()
    write_csv(tab, file.path(outdir, paste0(out_prefix, "_count_prop.csv")))
    
    tab1 <- tab %>% filter(pred=="1")
    p1 <- ggplot(tab1, aes(x=INFO_f, y=prop, fill=INFO_f)) +
      geom_col(width=0.6) + scale_y_continuous(labels=scales::percent) +
      labs(title=paste0(col, " = 1 proportion by group"), x="Class", y="Proportion") +
      theme_minimal() + theme(legend.position="none")
    ggsave_both(file.path(outdir, paste0(out_prefix, "_prop_pred1_bar.pdf")), p1, 6, 4)
    
    p2 <- ggplot(tab, aes(x=INFO_f, y=prop, fill=pred)) +
      geom_col(position="fill") + scale_y_continuous(labels=scales::percent) +
      labs(title=paste0(col, " distribution by group"), x="Class", y="Proportion", fill="Pred") +
      theme_minimal()
    ggsave_both(file.path(outdir, paste0(out_prefix, "_stacked_prop.pdf")), p2, 6, 4)
    
    return(list(type="binary", prop1=p1, stacked=p2))
    
  } else {
    p_box <- ggplot(dd, aes(x=lab_info(INFO), y=.data[[col]], fill=factor(INFO))) +
      geom_boxplot(alpha=0.7, outlier.shape=NA) +
      geom_jitter(width=0.2, alpha=0.3, size=0.8) +
      theme_minimal() +
      labs(title=paste0(col, " by Group"), x="Class", y=col) +
      theme(legend.position="none")
    ggsave_both(file.path(outdir, paste0(out_prefix, "_box.pdf")), p_box, 7, 5)
    
    p_den <- ggplot(dd, aes(x=.data[[col]], fill=lab_info(INFO))) +
      geom_density(alpha=0.4) + theme_minimal() +
      labs(title=paste0("Density of ", col, " by Group"), x=col, y="Density", fill="Group")
    ggsave_both(file.path(outdir, paste0(out_prefix, "_density.pdf")), p_den, 7, 5)
    
    return(list(type="continuous", box=p_box, density=p_den))
  }
}
res_eqtl <- analyze_pred_col(df, "pred_eqtl", "F1_pred_eqtl")
res_gwas <- analyze_pred_col(df, "pred_gwas", "F2_pred_gwas")

# ============================================================
# G) motif_RBPs / Rep_info 的“每变异条目个数”统计（0/1/2/≥3）
# ============================================================
count_items <- function(x, sep_guess = NULL) {
  if (is.na(x) || x=="") return(0L)
  s <- as.character(x)
  seps <- if (is.null(sep_guess)) c("__",";","\\|",",","_","\\s+") else sep_guess
  sep <- NULL
  if (is.null(sep_guess)) { for (sp in seps) { if (grepl(sp, s, perl=TRUE)) { sep <- sp; break } } }
  else sep <- sep_guess
  if (is.null(sep)) return(1L)
  parts <- unlist(strsplit(s, sep, perl=TRUE))
  sum(nzchar(trimws(parts)))
}

motif_items <- if ("motif_RBPs" %in% names(df)) {
  df %>% transmute(var_id, INFO,
                   n_motif_items = ifelse(is.na(motif_RBPs) | motif_RBPs=="",
                                          0L, sapply(motif_RBPs, count_items, sep_guess = "_")))
} else tibble(var_id=df$var_id, INFO=df$INFO, n_motif_items=NA_integer_)

rep_items <- if ("Rep_info" %in% names(df)) {
  df %>% transmute(var_id, INFO,
                   n_rep_items = ifelse(is.na(Rep_info) | Rep_info=="",
                                        0L, sapply(Rep_info, count_items, sep_guess = NULL)))
} else tibble(var_id=df$var_id, INFO=df$INFO, n_rep_items=NA_integer_)

items_sum <- motif_items %>% left_join(rep_items, by = c("var_id","INFO"))
write_csv(items_sum, file.path(outdir, "G_items_per_variant.csv"))

plot_items_bar <- function(df_items, count_col, title_prefix, out_name_base) {
  d1 <- df_items %>%
    transmute(bucket = case_when(
      .data[[count_col]] %in% 0 ~ "0",
      .data[[count_col]] %in% 1 ~ "1",
      .data[[count_col]] %in% 2 ~ "2",
      .data[[count_col]] >= 3 ~ "≥3",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(bucket)) %>% count(bucket, name="n") %>%
    mutate(prop = n / sum(n))
  p_all <- ggplot(d1, aes(bucket, prop, fill=bucket)) +
    geom_col() + scale_y_continuous(labels=scales::percent) +
    labs(title=paste0(title_prefix, " items per variant (overall)"),
         x="Items", y="Proportion", fill="Items") + theme_minimal()
  ggsave_both(file.path(outdir, paste0(out_name_base, "_overall.pdf")), p_all, 6, 4)
  
  d2 <- df_items %>%
    mutate(INFO_f = lab_info(INFO),
           bucket = case_when(
             .data[[count_col]] %in% 0 ~ "0",
             .data[[count_col]] %in% 1 ~ "1",
             .data[[count_col]] %in% 2 ~ "2",
             .data[[count_col]] >= 3 ~ "≥3",
             TRUE ~ NA_character_
           )) %>%
    filter(!is.na(bucket)) %>%
    count(INFO_f, bucket, name="n") %>%
    group_by(INFO_f) %>% mutate(prop = n/sum(n)) %>% ungroup()
  p_grp <- ggplot(d2, aes(INFO_f, prop, fill=bucket)) +
    geom_col(position="fill") +
    scale_y_continuous(labels=scales::percent) +
    labs(title=paste0(title_prefix, " items per variant (by Group)"),
         x="Class", y="Proportion", fill="Items") + theme_minimal()
  ggsave_both(file.path(outdir, paste0(out_name_base, "_byGroup.pdf")), p_grp, 7, 5)
  write_csv(d2, file.path(outdir, paste0(out_name_base, "_byGroup_table.csv")))
}

if (!all(is.na(items_sum$n_motif_items))) {
  plot_items_bar(items_sum, "n_motif_items", "motif_RBPs", "G1_motif_items")
} else message("motif_RBPs not found, skip motif items plots.")
if (!all(is.na(items_sum$n_rep_items))) {
  plot_items_bar(items_sum, "n_rep_items", "Rep_info", "G2_rep_items")
} else message("Rep_info not found/empty, skip rep items plots.")


