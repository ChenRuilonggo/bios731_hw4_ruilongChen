suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# =========================
# 1. 路径设置
# =========================
raw_dir <- "/projects/YangLabData/Ruilong/results/raw"
summary_dir <- "/projects/YangLabData/Ruilong/731/results/summary"
fig_dir <- "/projects/YangLabData/Ruilong/731/figures"

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 2. 读取所有 raw 结果
# =========================
raw_files <- list.files(
  path = raw_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(raw_files) == 0) {
  stop("No raw result files found in results/raw.")
}

all_results <- map_dfr(raw_files, read_csv, show_col_types = FALSE)

# 检查关键列
required_cols <- c(
  "method", "n", "sim_id", "component",
  "mu_true", "mu_hat", "ci_lower", "ci_upper",
  "covered", "elapsed_time"
)

missing_cols <- setdiff(required_cols, names(all_results))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# =========================
# 3. 去重检查（防止重复提交导致重复文件）
# =========================
dup_check <- all_results %>%
  count(method, n, sim_id, component) %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  warning("Duplicate rows detected. Keeping the first occurrence for each method-n-sim_id-component combination.")
  
  all_results <- all_results %>%
    distinct(method, n, sim_id, component, .keep_all = TRUE)
}

# =========================
# 4. 计算 bias 和 coverage 汇总
# =========================
summary_metrics <- all_results %>%
  mutate(
    bias_i = mu_hat - mu_true
  ) %>%
  group_by(method, n, component, mu_true) %>%
  summarise(
    nsim = n(),
    mean_mu_hat = mean(mu_hat),
    bias = mean(bias_i),
    bias_mcse = sd(mu_hat) / sqrt(nsim),
    coverage = mean(covered),
    coverage_mcse = sqrt(coverage * (1 - coverage) / nsim),
    mean_elapsed_time = mean(elapsed_time),
    median_elapsed_time = median(elapsed_time),
    .groups = "drop"
  )

# =========================
# 5. 计算时间汇总
# =========================
time_summary <- all_results %>%
  group_by(method, n, sim_id) %>%
  summarise(
    elapsed_time = first(elapsed_time),
    .groups = "drop"
  ) %>%
  group_by(method, n) %>%
  summarise(
    nsim = n(),
    mean_elapsed_time = mean(elapsed_time),
    sd_elapsed_time = sd(elapsed_time),
    median_elapsed_time = median(elapsed_time),
    min_elapsed_time = min(elapsed_time),
    max_elapsed_time = max(elapsed_time),
    .groups = "drop"
  )

# =========================
# 6. 保存汇总结果
# =========================
write_csv(all_results, file.path(summary_dir, "all_results_combined.csv"))
write_csv(summary_metrics, file.path(summary_dir, "summary_metrics.csv"))
write_csv(time_summary, file.path(summary_dir, "time_summary.csv"))

# =========================
# 7. Bias plot with MCSE error bars
# =========================
bias_plot <- ggplot(
  summary_metrics,
  aes(
    x = factor(mu_true),
    y = bias,
    color = method,
    group = method
  )
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(
    position = position_dodge(width = 0.35),
    size = 2
  ) +
  geom_errorbar(
    aes(
      ymin = bias - 1.96 * bias_mcse,
      ymax = bias + 1.96 * bias_mcse
    ),
    width = 0.15,
    position = position_dodge(width = 0.35)
  ) +
  facet_wrap(~ n, scales = "free_y") +
  labs(
    title = "Bias of mu-hat with Monte Carlo standard error bars",
    x = expression(True~mu[k]),
    y = "Bias"
  ) +
  theme_bw()

ggsave(
  filename = file.path(fig_dir, "bias_plot.png"),
  plot = bias_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# =========================
# 8. Coverage plot with MCSE error bars
# =========================
coverage_plot <- ggplot(
  summary_metrics,
  aes(
    x = factor(mu_true),
    y = coverage,
    color = method,
    group = method
  )
) +
  geom_hline(yintercept = 0.95, linetype = 2) +
  geom_point(
    position = position_dodge(width = 0.35),
    size = 2
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(0, coverage - 1.96 * coverage_mcse),
      ymax = pmin(1, coverage + 1.96 * coverage_mcse)
    ),
    width = 0.15,
    position = position_dodge(width = 0.35)
  ) +
  facet_wrap(~ n) +
  labs(
    title = "Coverage of 95% intervals with Monte Carlo standard error bars",
    x = expression(True~mu[k]),
    y = "Coverage"
  ) +
  theme_bw()

ggsave(
  filename = file.path(fig_dir, "coverage_plot.png"),
  plot = coverage_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# =========================
# 9. 终端打印简要结果
# =========================
cat("\nSaved summary files to:\n")
cat(summary_dir, "\n")

cat("\nSaved figures to:\n")
cat(fig_dir, "\n")

cat("\nTime summary:\n")
print(time_summary)

cat("\nSummary metrics (first 12 rows):\n")
print(head(summary_metrics, 12))