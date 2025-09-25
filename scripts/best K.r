# Load libraries
library(ggplot2)
library(dplyr)
library(patchwork)

# -------------------------
# 1. Read space-delimited input file
# -------------------------
data_file <- "ofavNgsAdmixLogfile_formatted"  # replace with your file
df <- read.table(data_file, header = FALSE, sep = "", 
                 col.names = c("K", "LogLikelihood"), stringsAsFactors = FALSE)

# -------------------------
# 2. Clean columns and convert to numeric
# -------------------------
df$K <- as.numeric(trimws(df$K))
df$LogLikelihood <- as.numeric(trimws(df$LogLikelihood))
df <- df[!is.na(df$K) & !is.na(df$LogLikelihood), ]

if(nrow(df) == 0) stop("No valid data found in the file. Check file formatting.")

# -------------------------
# 3. Summarize mean and SD for each K
# -------------------------
summary_df <- df %>%
  group_by(K) %>%
  summarise(MeanLnP = mean(LogLikelihood, na.rm = TRUE),
            SD = sd(LogLikelihood, na.rm = TRUE)) %>%
  arrange(K)

# -------------------------
# 4. Calculate ΔK robustly
# -------------------------
L <- summary_df$MeanLnP
n <- nrow(summary_df)
deltaK <- rep(NA, n)

if(n >= 3) {
  for(i in 2:(n-1)) {
    Li_minus1 <- L[i-1]; Li <- L[i]; Li_plus1 <- L[i+1]; SDi <- summary_df$SD[i]
    if(all(is.finite(c(Li_minus1, Li, Li_plus1, SDi))) && SDi != 0) {
      deltaK[i] <- abs(Li_plus1 - 2*Li + Li_minus1) / SDi
    }
  }
}

summary_df$deltaK <- deltaK

# -------------------------
# 5. Identify all ΔK peaks
# -------------------------
is_peak <- rep(FALSE, n)
if(n >= 3 && any(!is.na(deltaK))) {
  for(i in 2:(n-1)) {
    if(!is.na(deltaK[i])) {
      left <- deltaK[i-1]; right <- deltaK[i+1]
      if((is.na(left) || deltaK[i] > left) && (is.na(right) || deltaK[i] > right)) {
        is_peak[i] <- TRUE
      }
    }
  }
}
peaks_df <- summary_df[is_peak, ]

# Identify global best K (largest ΔK)
global_best <- if(nrow(peaks_df) > 0) peaks_df[which.max(peaks_df$deltaK), ] else NULL

# -------------------------
# 6. Create plots
# -------------------------

# ΔK plot
if(!all(is.na(summary_df$deltaK))) {
  p_deltaK <- ggplot(summary_df, aes(x = K, y = deltaK)) +
    geom_line(color = "steelblue") +
    geom_point() +
    theme_minimal() +
    labs(title = "ΔK Plot (Evanno Method)",
         x = "K",
         y = expression(Delta*K)) +
    coord_cartesian(clip = "off")  # prevent label clipping
  
  if(nrow(peaks_df) > 0) {
    p_deltaK <- p_deltaK +
      geom_point(data = peaks_df, aes(x = K, y = deltaK), color = "red", size = 3) +
      geom_text(data = peaks_df, aes(x = K, y = deltaK, label = paste("K =", K)),
                vjust = -0.5, color = "red", fontface = "bold", clip = "off")
  }
  
  if(!is.null(global_best)) {
    p_deltaK <- p_deltaK +
      geom_point(data = global_best, aes(x = K, y = deltaK), color = "darkred", size = 5) +
      geom_text(data = global_best, aes(x = K, y = deltaK, label = paste("Best K =", K)),
                vjust = -1, color = "darkred", fontface = "bold", clip = "off")
  }
} else {
  p_deltaK <- NULL
}

# Mean LogLikelihood plot
p_mean <- ggplot(summary_df, aes(x = K, y = MeanLnP)) +
  geom_line(color = "darkgreen") +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Mean Log Likelihood vs K",
       x = "K",
       y = "Mean Log Likelihood") +
  coord_cartesian(clip = "off")

# Combine plots vertically
if(!is.null(p_deltaK)) {
  p_combined <- p_deltaK / p_mean + plot_layout(heights = c(1,1))
} else {
  p_combined <- p_mean
}

# -------------------------
# 7. Save combined plot
# -------------------------
output_file <- "DeltaK_and_MeanLnP_plot.png"
ggsave(output_file, plot = p_combined, width = 6, height = 8, dpi = 300)

# -------------------------
# 8. Print peak info
# -------------------------
cat("Plot saved as", output_file, "\n")
if(!is.null(global_best)) {
  cat("ΔK peaks at K =", paste(peaks_df$K, collapse = ", "), "\n")
  cat("ΔK values:", paste(round(peaks_df$deltaK, 2), collapse = ", "), "\n")
  cat("Global best K =", global_best$K, "with ΔK =", round(global_best$deltaK, 2), "\n")
} else if(nrow(summary_df) > 0) {
  cat("ΔK could not be calculated; showing Mean Log Likelihood plot.\n")
} else {
  cat("No valid data to plot.\n")
}
