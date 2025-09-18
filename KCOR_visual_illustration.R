# KCOR (Kirsch Cumulative Outcomes Ratio) Visual Illustration
# This script creates geometric and visual illustrations of the KCOR method
# based on the method described in README.md

# Load required libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(cowplot)

# Set theme for all plots
theme_set(theme_minimal() + 
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5, size = 12),
                axis.title = element_text(size = 11),
                legend.position = "bottom"))

# Function to generate synthetic mortality data with realistic patterns
generate_mortality_data <- function(weeks = 120, enrollment_week = 0) {
  
  # Base mortality rates (weekly) - realistic values
  # Unvaccinated cohort: baseline exponential trend
  base_mortality_unvax <- 0.002  # 0.2% weekly baseline mortality
  slope_unvax <- 0.001  # slight upward trend
  
  # Vaccinated cohort: initially higher (immediate vaccine effect), then similar slope
  base_mortality_vax <- 0.002
  slope_vax <- 0.001
  vaccine_effect <- 0.0008  # additional mortality from vaccination
  
  # Create time series
  time <- 0:weeks
  
  # Raw mortality rates with some noise
  set.seed(123)
  mr_unvax_raw <- base_mortality_unvax * exp(slope_unvax * time/52) + 
                  rnorm(length(time), 0, 0.0003)
  mr_vax_raw <- (base_mortality_unvax + vaccine_effect) * exp(slope_vax * time/52) + 
                rnorm(length(time), 0, 0.0003)
  
  # Ensure positive rates
  mr_unvax_raw <- pmax(mr_unvax_raw, 0.0001)
  mr_vax_raw <- pmax(mr_vax_raw, 0.0001)
  
  # Apply 8-week centered moving average (smoothing)
  smooth_ma <- function(x, n = 8) {
    stats::filter(x, rep(1/n, n), sides = 2)
  }
  
  mr_unvax_smooth <- as.numeric(smooth_ma(mr_unvax_raw))
  mr_vax_smooth <- as.numeric(smooth_ma(mr_vax_raw))
  
  # Fill NAs at edges with original values
  mr_unvax_smooth[is.na(mr_unvax_smooth)] <- mr_unvax_raw[is.na(mr_unvax_smooth)]
  mr_vax_smooth[is.na(mr_vax_smooth)] <- mr_vax_raw[is.na(mr_vax_smooth)]
  
  # Calculate slopes using anchor points (weeks 53 and 114 as per methodology)
  anchor1 <- 53; anchor2 <- 114
  if (anchor2 <= length(time)) {
    # Geometric mean around anchor points (±2 weeks window)
    window <- 2
    
    # Anchor 1 windows
    idx1_unvax <- max(1, anchor1-window):min(length(mr_unvax_smooth), anchor1+window)
    idx1_vax <- max(1, anchor1-window):min(length(mr_vax_smooth), anchor1+window)
    gm1_unvax <- exp(mean(log(mr_unvax_smooth[idx1_unvax])))
    gm1_vax <- exp(mean(log(mr_vax_smooth[idx1_vax])))
    
    # Anchor 2 windows
    idx2_unvax <- max(1, anchor2-window):min(length(mr_unvax_smooth), anchor2+window)
    idx2_vax <- max(1, anchor2-window):min(length(mr_vax_smooth), anchor2+window)
    gm2_unvax <- exp(mean(log(mr_unvax_smooth[idx2_unvax])))
    gm2_vax <- exp(mean(log(mr_vax_smooth[idx2_vax])))
    
    # Calculate slopes
    delta_t <- anchor2 - anchor1
    slope_calc_unvax <- (1/delta_t) * log(gm2_unvax / gm1_unvax)
    slope_calc_vax <- (1/delta_t) * log(gm2_vax / gm1_vax)
  } else {
    # Use simpler slope calculation if anchors are outside range
    slope_calc_unvax <- slope_unvax
    slope_calc_vax <- slope_vax
  }
  
  return(list(
    time = time,
    mr_unvax_raw = mr_unvax_raw,
    mr_vax_raw = mr_vax_raw,
    mr_unvax_smooth = mr_unvax_smooth,
    mr_vax_smooth = mr_vax_smooth,
    slope_unvax = slope_calc_unvax,
    slope_vax = slope_calc_vax,
    anchor1 = anchor1,
    anchor2 = anchor2
  ))
}

# Generate synthetic data
data <- generate_mortality_data(weeks = 120)

# Step 1: Raw vs Smoothed Mortality Rates
create_raw_smoothed_plot <- function(data) {
  df <- data.frame(
    Week = rep(data$time, 4),
    MortalityRate = c(data$mr_unvax_raw, data$mr_vax_raw, 
                      data$mr_unvax_smooth, data$mr_vax_smooth),
    Cohort = rep(c("Unvaccinated", "Vaccinated"), each = length(data$time), times = 2),
    Type = rep(c("Raw", "Raw", "Smoothed (8-wk MA)", "Smoothed (8-wk MA)"), 
               each = length(data$time))
  )
  
  p <- ggplot(df, aes(x = Week, y = MortalityRate, color = Cohort, linetype = Type)) +
    geom_line(size = 0.8) +
    scale_color_manual(values = c("Unvaccinated" = "#2E86AB", "Vaccinated" = "#A23B72")) +
    scale_linetype_manual(values = c("Raw" = "dotted", "Smoothed (8-wk MA)" = "solid")) +
    labs(title = "Step 1: Raw vs Smoothed Mortality Rates",
         subtitle = "8-week centered moving average reduces noise",
         x = "Week from Enrollment",
         y = "Weekly Mortality Rate",
         color = "Cohort", linetype = "Data Type") +
    geom_vline(xintercept = c(data$anchor1, data$anchor2), 
               color = "red", linetype = "dashed", alpha = 0.7) +
    annotate("text", x = data$anchor1, y = max(df$MortalityRate) * 0.9, 
             label = "Anchor 1", angle = 90, vjust = -0.5, color = "red") +
    annotate("text", x = data$anchor2, y = max(df$MortalityRate) * 0.9, 
             label = "Anchor 2", angle = 90, vjust = -0.5, color = "red")
  
  return(p)
}

# Step 2: Slope Normalization
create_slope_normalization_plot <- function(data) {
  # Apply slope correction: MR_adj(t) = MR(t) * exp(-r(t - t0))
  # where t0 = 0 (enrollment week)
  
  mr_unvax_adj <- data$mr_unvax_smooth * exp(-data$slope_unvax * (data$time - 0))
  mr_vax_adj <- data$mr_vax_smooth * exp(-data$slope_vax * (data$time - 0))
  
  df <- data.frame(
    Week = rep(data$time, 4),
    MortalityRate = c(data$mr_unvax_smooth, data$mr_vax_smooth,
                      mr_unvax_adj, mr_vax_adj),
    Cohort = rep(c("Unvaccinated", "Vaccinated"), each = length(data$time), times = 2),
    Type = rep(c("Original", "Original", "Slope Adjusted", "Slope Adjusted"), 
               each = length(data$time))
  )
  
  p <- ggplot(df, aes(x = Week, y = MortalityRate, color = Cohort, linetype = Type)) +
    geom_line(size = 0.8) +
    scale_color_manual(values = c("Unvaccinated" = "#2E86AB", "Vaccinated" = "#A23B72")) +
    scale_linetype_manual(values = c("Original" = "dashed", "Slope Adjusted" = "solid")) +
    labs(title = "Step 2: Slope Normalization",
         subtitle = "MR_adj(t) = MR(t) × exp(-r(t - t₀))",
         x = "Week from Enrollment",
         y = "Weekly Mortality Rate",
         color = "Cohort", linetype = "Adjustment") +
    geom_vline(xintercept = 0, color = "green", linetype = "dotted", alpha = 0.7) +
    annotate("text", x = 5, y = max(df$MortalityRate) * 0.9, 
             label = "Enrollment\n(t₀ = 0)", color = "green")
  
  return(p)
}

# Step 3: Discrete Hazard Function Transform
create_hazard_transform_plot <- function(data) {
  # Calculate adjusted mortality rates
  mr_unvax_adj <- data$mr_unvax_smooth * exp(-data$slope_unvax * (data$time - 0))
  mr_vax_adj <- data$mr_vax_smooth * exp(-data$slope_vax * (data$time - 0))
  
  # Clip to prevent log(0): MR_adj clipped to 0.999
  mr_unvax_adj <- pmin(mr_unvax_adj, 0.999)
  mr_vax_adj <- pmin(mr_vax_adj, 0.999)
  
  # Apply discrete hazard transform: hazard(t) = -ln(1 - MR_adj(t))
  hazard_unvax <- -log(1 - mr_unvax_adj)
  hazard_vax <- -log(1 - mr_vax_adj)
  
  # Create comparison plot
  df <- data.frame(
    Week = rep(data$time, 4),
    Value = c(mr_unvax_adj, mr_vax_adj, hazard_unvax, hazard_vax),
    Cohort = rep(c("Unvaccinated", "Vaccinated"), each = length(data$time), times = 2),
    Transform = rep(c("Adjusted MR", "Adjusted MR", "Discrete Hazard", "Discrete Hazard"), 
                    each = length(data$time))
  )
  
  p <- ggplot(df, aes(x = Week, y = Value, color = Cohort)) +
    geom_line(size = 0.8) +
    facet_wrap(~Transform, scales = "free_y", ncol = 1,
               labeller = labeller(Transform = c("Adjusted MR" = "Adjusted Mortality Rate",
                                                "Discrete Hazard" = "Discrete Hazard: -ln(1 - MR_adj)"))) +
    scale_color_manual(values = c("Unvaccinated" = "#2E86AB", "Vaccinated" = "#A23B72")) +
    labs(title = "Step 3: Discrete Hazard Function Transform",
         subtitle = "Mathematical exactness through hazard transformation",
         x = "Week from Enrollment",
         y = "Value",
         color = "Cohort")
  
  return(p)
}

# Step 4: Cumulative Hazards and KCOR
create_kcor_plot <- function(data) {
  # Calculate adjusted mortality rates and hazards
  mr_unvax_adj <- data$mr_unvax_smooth * exp(-data$slope_unvax * (data$time - 0))
  mr_vax_adj <- data$mr_vax_smooth * exp(-data$slope_vax * (data$time - 0))
  
  mr_unvax_adj <- pmin(mr_unvax_adj, 0.999)
  mr_vax_adj <- pmin(mr_vax_adj, 0.999)
  
  hazard_unvax <- -log(1 - mr_unvax_adj)
  hazard_vax <- -log(1 - mr_vax_adj)
  
  # Calculate cumulative hazards: CH(t) = sum(hazard(i)) for i=0 to t
  ch_unvax <- cumsum(hazard_unvax)
  ch_vax <- cumsum(hazard_vax)
  
  # Calculate KCOR with baseline normalization at week 4
  baseline_week <- 5  # Week 4 (0-indexed, so position 5)
  if (length(ch_vax) >= baseline_week && length(ch_unvax) >= baseline_week) {
    baseline_ratio <- ch_vax[baseline_week] / ch_unvax[baseline_week]
    kcor <- (ch_vax / ch_unvax) / baseline_ratio
  } else {
    kcor <- ch_vax / ch_unvax
  }
  
  # Create plots
  # Subplot 1: Cumulative Hazards
  df_ch <- data.frame(
    Week = rep(data$time, 2),
    CumulativeHazard = c(ch_unvax, ch_vax),
    Cohort = rep(c("Unvaccinated", "Vaccinated"), each = length(data$time))
  )
  
  p1 <- ggplot(df_ch, aes(x = Week, y = CumulativeHazard, color = Cohort)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Unvaccinated" = "#2E86AB", "Vaccinated" = "#A23B72")) +
    labs(title = "Cumulative Hazards: CH(t) = Σhazard(i)",
         x = "Week from Enrollment",
         y = "Cumulative Hazard",
         color = "Cohort") +
    geom_vline(xintercept = 4, color = "orange", linetype = "dashed", alpha = 0.7) +
    annotate("text", x = 4, y = max(df_ch$CumulativeHazard) * 0.9, 
             label = "Baseline\n(Week 4)", angle = 90, vjust = -0.5, color = "orange")
  
  # Subplot 2: KCOR Values
  df_kcor <- data.frame(
    Week = data$time,
    KCOR = kcor
  )
  
  p2 <- ggplot(df_kcor, aes(x = Week, y = KCOR)) +
    geom_line(size = 1, color = "#E63946") +
    geom_hline(yintercept = 1, linetype = "solid", color = "black", alpha = 0.5) +
    labs(title = "KCOR(t) = [CH_vax(t)/CH_unvax(t)] / [CH_vax(t₀)/CH_unvax(t₀)]",
         x = "Week from Enrollment",
         y = "KCOR Value") +
    geom_vline(xintercept = 4, color = "orange", linetype = "dashed", alpha = 0.7) +
    annotate("text", x = 4, y = max(df_kcor$KCOR) * 0.9, 
             label = "Baseline\n(Week 4)", angle = 90, vjust = -0.5, color = "orange") +
    annotate("text", x = max(data$time) * 0.7, y = 1.05, 
             label = "KCOR = 1.0\n(No difference)", color = "black", size = 3)
  
  # Combine plots
  combined <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 1))
  
  return(combined)
}

# Create KCOR method flowchart
create_methodology_flowchart <- function() {
  # Create a simple flowchart showing the 4 main steps
  library(grid)
  
  # Create a new plot
  p <- ggplot() + 
    xlim(0, 10) + ylim(0, 10) +
    theme_void() +
    labs(title = "KCOR Flowchart",
         subtitle = "Four-Step Process for Mortality Risk Assessment")
  
  # Add boxes and text for each step
  # Step 1
  p <- p + 
    geom_rect(aes(xmin = 1, xmax = 3, ymin = 8, ymax = 9), 
              fill = "#E3F2FD", color = "#1976D2", size = 1) +
    annotate("text", x = 2, y = 8.5, label = "Step 1:\nSmoothing\n(8-week MA)", 
             hjust = 0.5, vjust = 0.5, size = 3)
  
  # Step 2
  p <- p + 
    geom_rect(aes(xmin = 1, xmax = 3, ymin = 6, ymax = 7), 
              fill = "#F3E5F5", color = "#7B1FA2", size = 1) +
    annotate("text", x = 2, y = 6.5, label = "Step 2:\nSlope Adjustment\nMR_adj(t)", 
             hjust = 0.5, vjust = 0.5, size = 3)
  
  # Step 3
  p <- p + 
    geom_rect(aes(xmin = 1, xmax = 3, ymin = 4, ymax = 5), 
              fill = "#E8F5E8", color = "#388E3C", size = 1) +
    annotate("text", x = 2, y = 4.5, label = "Step 3:\nHazard Transform\n-ln(1 - MR_adj)", 
             hjust = 0.5, vjust = 0.5, size = 3)
  
  # Step 4
  p <- p + 
    geom_rect(aes(xmin = 1, xmax = 3, ymin = 2, ymax = 3), 
              fill = "#FFEBEE", color = "#D32F2F", size = 1) +
    annotate("text", x = 2, y = 2.5, label = "Step 4:\nKCOR Ratio\nCH_vax/CH_unvax", 
             hjust = 0.5, vjust = 0.5, size = 3)
  
  # Add arrows
  p <- p + 
    geom_segment(aes(x = 2, y = 7.9, xend = 2, yend = 7.1), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 1, color = "black") +
    geom_segment(aes(x = 2, y = 5.9, xend = 2, yend = 5.1), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 1, color = "black") +
    geom_segment(aes(x = 2, y = 3.9, xend = 2, yend = 3.1), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 1, color = "black")
  
  # Add key formulas on the right
  p <- p + 
    annotate("text", x = 6.5, y = 8.5, 
             label = "Smoothing:\n8-week centered moving average", 
             hjust = 0, size = 3, color = "#1976D2") +
    annotate("text", x = 6.5, y = 6.5, 
             label = "Slope Adjustment:\nMR_adj(t) = MR(t) × exp(-r(t - t₀))", 
             hjust = 0, size = 3, color = "#7B1FA2") +
    annotate("text", x = 6.5, y = 4.5, 
             label = "Discrete Hazard:\nhazard(t) = -ln(1 - MR_adj(t))", 
             hjust = 0, size = 3, color = "#388E3C") +
    annotate("text", x = 6.5, y = 2.5, 
             label = "KCOR Formula:\nKCOR(t) = [CH_vax(t)/CH_unvax(t)] / baseline", 
             hjust = 0, size = 3, color = "#D32F2F")
  
  return(p)
}

# Generate all plots
cat("Generating KCOR visual illustrations...\n")

# Create individual plots
plot1 <- create_raw_smoothed_plot(data)
plot2 <- create_slope_normalization_plot(data)
plot3 <- create_hazard_transform_plot(data)
plot4 <- create_kcor_plot(data)
plot5 <- create_methodology_flowchart()

# Save individual plots
ggsave("KCOR_step1_smoothing.png", plot1, width = 12, height = 8, dpi = 300)
ggsave("KCOR_step2_slope_normalization.png", plot2, width = 12, height = 8, dpi = 300)
ggsave("KCOR_step3_hazard_transform.png", plot3, width = 12, height = 10, dpi = 300)
ggsave("KCOR_step4_cumulative_hazard_kcor.png", plot4, width = 12, height = 10, dpi = 300)
ggsave("KCOR_methodology_flowchart.png", plot5, width = 12, height = 10, dpi = 300)

# Create comprehensive summary plot
summary_plot <- plot_grid(
  plot5,
  plot_grid(plot1, plot2, ncol = 1),
  plot_grid(plot3, plot4, ncol = 1),
  ncol = 3,
  rel_widths = c(1, 1, 1),
  labels = c("A", "B", "C"),
  label_size = 16
)

ggsave("KCOR_complete_methodology_illustration.png", summary_plot, 
       width = 20, height = 12, dpi = 300)

cat("Visualizations saved successfully!\n")
cat("Files created:\n")
cat("- KCOR_step1_smoothing.png\n")
cat("- KCOR_step2_slope_normalization.png\n") 
cat("- KCOR_step3_hazard_transform.png\n")
cat("- KCOR_step4_cumulative_hazard_kcor.png\n")
cat("- KCOR_methodology_flowchart.png\n")
cat("- KCOR_complete_methodology_illustration.png\n")

# Print key insights
cat("\n=== KEY KCOR METHOD INSIGHTS ===\n")
cat("1. SMOOTHING: 8-week centered moving average reduces noise in raw mortality data\n")
cat("2. SLOPE NORMALIZATION: Adjusts for baseline mortality trends: MR_adj(t) = MR(t) × exp(-r(t-t₀))\n")
cat("3. HAZARD TRANSFORM: Converts to discrete hazard functions for mathematical exactness\n")
cat("4. KCOR CALCULATION: Ratio of cumulative hazards, normalized to baseline (week 4)\n")
cat("\nKCOR > 1.0: Higher mortality in vaccinated group\n")
cat("KCOR < 1.0: Lower mortality in vaccinated group\n")
cat("KCOR = 1.0: No difference between groups\n")

# Display final KCOR values from the simulation
mr_unvax_adj <- data$mr_unvax_smooth * exp(-data$slope_unvax * (data$time - 0))
mr_vax_adj <- data$mr_vax_smooth * exp(-data$slope_vax * (data$time - 0))
mr_unvax_adj <- pmin(mr_unvax_adj, 0.999)
mr_vax_adj <- pmin(mr_vax_adj, 0.999)
hazard_unvax <- -log(1 - mr_unvax_adj)
hazard_vax <- -log(1 - mr_vax_adj)
ch_unvax <- cumsum(hazard_unvax)
ch_vax <- cumsum(hazard_vax)
baseline_week <- 5
if (length(ch_vax) >= baseline_week && length(ch_unvax) >= baseline_week) {
  baseline_ratio <- ch_vax[baseline_week] / ch_unvax[baseline_week]
  kcor <- (ch_vax / ch_unvax) / baseline_ratio
  final_kcor <- tail(kcor, 1)
  cat(sprintf("\nSimulated Final KCOR Value: %.4f\n", final_kcor))
  cat(sprintf("This represents a %.1f%% %s in mortality risk for the vaccinated group\n", 
              abs(final_kcor - 1) * 100,
              ifelse(final_kcor > 1, "increase", "decrease")))
}