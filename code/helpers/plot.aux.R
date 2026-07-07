library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(MASS)  # for kde2d
library(viridis)  # for nice color scales

# Function to compute posterior probability that pos > neg at each horizon
compute_prob_pos_gt_neg <- function(pos_shock, neg_shock) {
  H <- ncol(pos_shock)
  sapply(1:H, function(h) mean(median(pos_shock[, h]) < (neg_shock[, h])))
}

plot_irf_with_prob_labels_top <- function(pos_shock, neg_shock) {
  H <- ncol(pos_shock)
  horizon <- 0:(H - 1)
  
  # Step 1: Summarize IRFs
  get_summary <- function(shock_mat, label) {
    data.frame(
      horizon = horizon,
      t(apply(shock_mat, 2, quantile, probs = c(0.16, 0.5, 0.84))),
      shock = label
    ) |> 
      rename(p16 = X16., median = X50., p84 = X84.)
  }
  
  pos_df <- get_summary(pos_shock, "positive")
  neg_df <- get_summary(neg_shock, "negative")
  plot_df <- bind_rows(pos_df, neg_df)
  
  # Step 2: Compute posterior probabilities
  probs <- compute_prob_pos_gt_neg(pos_shock, neg_shock)
  
  # Step 3: Position label and line
  pmax_ci <- sapply(1:H, function(h) max(pos_df$p84[h], neg_df$p84[h]))
  label_y <- pmax_ci + 0.05 * abs(pmax_ci)
  
  label_df <- data.frame(
    horizon = horizon,
    label = sprintf("%.2f", probs),
    y = label_y
  )
  
  line_df <- data.frame(
    horizon = horizon,
    y_start = label_y,
    y_end = 0  # always down to x-axis
  )
  
  # Step 4: Plot
  ggplot(plot_df, aes(x = horizon)) +
    # IRF ribbons (drawn first, below everything)
    geom_ribbon(aes(ymin = p16, ymax = p84, fill = shock), alpha = 0.2, color = NA) +
    # IRF lines
    geom_line(aes(y = median, color = shock), linewidth = 1) +
    # Dashed vertical lines — drawn after IRFs to be visible
    geom_segment(data = line_df,
                 aes(x = horizon, xend = horizon, y = y_start, yend = y_end),
                 inherit.aes = FALSE, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    # Probability labels on top
    geom_text(data = label_df, aes(x = horizon, y = y, label = label),
              inherit.aes = FALSE, size = 3, vjust = 0) +
    # Styling
    scale_color_manual(values = c("positive" = "blue", "negative" = "red")) +
    scale_fill_manual(values = c("positive" = "blue", "negative" = "red")) +
    labs(x = "Horizon", y = "Impulse Response") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top")
}
test_irf_difference <- function(pos_shock, neg_shock, alpha = 0.05, method = c("wilcoxon", "ks", "mean")) {
  method <- match.arg(method)
  H <- ncol(pos_shock)
  sig_vector <- integer(H)
  
  for (h in 1:H) {
    x <- pos_shock[, h]
    y <- neg_shock[, h]
    
    if (method == "wilcoxon") {
      test_result <- wilcox.test(x, y, exact = FALSE)
      sig_vector[h] <- as.integer(test_result$p.value < alpha)
      
    } else if (method == "ks") {
      test_result <- ks.test(x, y)
      sig_vector[h] <- as.integer(test_result$p.value < alpha)
      
    } else if (method == "mean") {
      diff_samples <- x - y
      ci <- quantile(diff_samples, probs = c(alpha / 2, 1 - alpha / 2))
      sig_vector[h] <- as.integer(ci[1] > 0 | ci[2] < 0)
    }
  }
  
  return(sig_vector)
}

plot_irf_with_asterisks <- function(pos_shock, neg_shock, diff.stat) {
  H <- ncol(pos_shock)
  horizon <- 0:(H - 1)
  
  # Compute 16th, 50th, and 84th percentiles
  get_summary <- function(shock_mat, label) {
    as.data.frame(t(apply(shock_mat, 2, quantile, probs = c(0.16, 0.5, 0.84)))) %>%
      mutate(horizon = horizon, shock = label) %>%
      rename(p16 = `16%`, median = `50%`, p84 = `84%`)
  }
  
  pos_df <- get_summary(pos_shock, "positive")
  neg_df <- get_summary(neg_shock, "negative")
  plot_df <- bind_rows(pos_df, neg_df)
  
  # Data frame for significance asterisks
  star_df <- plot_df %>%
    filter(shock == "positive") %>%  # only need one line for horizon
    mutate(significant = diff.stat,
           label = ifelse(significant == 1, "*", NA),
           y = min(p16) - 0.05 * max(abs(p84)))  # place slightly below IRF area
  
  # Plot
  ggplot(plot_df, aes(x = horizon, y = median, color = shock, fill = shock)) +
    geom_ribbon(aes(ymin = p16, ymax = p84), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    geom_text(data = star_df, aes(x = horizon, y = y, label = label),
              inherit.aes = FALSE, size = 4, vjust = 1) +
    theme_minimal(base_size = 13) +
    labs(x = "Horizon", y = "Impulse Response") +
    scale_color_manual(values = c("positive" = "blue", "negative" = "red")) +
    scale_fill_manual(values = c("positive" = "blue", "negative" = "red")) +
    theme(legend.position = "top")
}

identify_factor_loadings <- function(loadings_posterior) {
  N <- dim(loadings_posterior)[1]  # draws
  M <- dim(loadings_posterior)[2]  # series
  Q <- dim(loadings_posterior)[3]  # factors
  
  aligned_loadings <- array(NA, dim = c(N, M, Q))
  
  # Step 1: Posterior mean as reference
  reference <- apply(loadings_posterior, c(2, 3), mean)  # M x Q
  reference <- scale(reference, center = FALSE, scale = apply(reference, 2, function(x) sqrt(sum(x^2))))
  
  for (i in 1:N) {
    current <- loadings_posterior[i, , ]  # M x Q
    current_norm <- scale(current, center = FALSE, scale = apply(current, 2, function(x) sqrt(sum(x^2))))
    
    # Compute absolute correlation matrix
    corr_mat <- t(reference) %*% current_norm  # Q x Q
    abs_corr <- abs(corr_mat)
    
    # Make a cost matrix with non-negative entries
    cost_mat <- max(abs_corr) - abs_corr  # ensure non-negativity
    
    # Hungarian algorithm for optimal column permutation
    perm <- clue::solve_LSAP(cost_mat)
    
    aligned <- matrix(NA, M, Q)
    for (j in 1:Q) {
      col_idx <- perm[j]
      sign_flip <- sign(sum(reference[, j] * current_norm[, col_idx]))
      aligned[, j] <- sign_flip * current[, col_idx]
    }
    
    aligned_loadings[i, , ] <- aligned
  }
  
  return(aligned_loadings)
}
