
library(tidyverse)
library(ggplot2)
library(quadprog)

install.packages("tidysynth")
data(smoking, package = "tidysynth")
smoking
view(smoking)

state_ids <- unique(smoking$state)
smoking$state_id <- as.numeric(factor(smoking$state, levels = state_ids))


california_id <- smoking$state_id[smoking$state == "California"][1]


outcome_wide <- smoking %>%
  select(state_id, year, cigsale) %>%
  pivot_wider(names_from = year, values_from = cigsale) %>%
  column_to_rownames("state_id")


treated_id <- california_id
control_ids <- setdiff(1:nrow(outcome_wide), treated_id)


pre_years <- 1970:1987
post_years <- 1989:2000


Y_pre_treated <- as.numeric(outcome_wide[treated_id, as.character(pre_years)])
Y_pre_controls <- as.matrix(outcome_wide[control_ids, as.character(pre_years)])

n_controls <- nrow(Y_pre_controls)


Dmat <- Y_pre_controls %*% t(Y_pre_controls)
ridge <- 1e-5 * diag(n_controls)
Dmat <- Dmat + ridge

dvec <- Y_pre_controls %*% Y_pre_treated

Amat <- cbind(rep(1, n_controls), diag(n_controls))
bvec <- c(1, rep(0, n_controls))
meq <- 1

qp_solution <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
w_opt <- qp_solution$solution
names(w_opt) <- rownames(outcome_wide)[control_ids]

cat("Sum of weights:", sum(w_opt), "\n")
print(round(sort(w_opt[w_opt > 1e-4], decreasing = TRUE), 4))


Y_full_controls <- as.matrix(outcome_wide[control_ids, ])
synthetic <- colSums(w_opt * Y_full_controls)
actual <- as.numeric(outcome_wide[treated_id, ])

years <- as.numeric(colnames(outcome_wide))
plot_df <- data.frame(
  year = years,
  actual = actual,
  synthetic = synthetic,
  gap = actual - synthetic
)

p1 <- ggplot(plot_df, aes(x = year)) +
  geom_line(aes(y = actual, color = "Actual California")) +
  geom_line(aes(y = synthetic, color = "Synthetic California")) +
  geom_vline(xintercept = 1988, linetype = "dashed", color = "red") +
  labs(title = "California Smoking Ban: Actual vs Synthetic Cigarette Sales",
       x = "Year", y = "Cigarette Sales (packs per capita)", color = "Series") +
  theme_minimal()
print(p1)
ggsave("california_actual_vs_synthetic.png", width = 8, height = 5)

p2 <- ggplot(plot_df, aes(x = year, y = gap)) +
  geom_line(color = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 1988, linetype = "dashed", color = "red") +
  labs(title = "Estimated Effect of Proposition 99 (Gap: Actual - Synthetic)",
       x = "Year", y = "Reduction in Cigarette Sales") +
  theme_minimal()
print(p2)
ggsave("california_treatment_gap.png", width = 8, height = 5)


compute_synthetic_for_state <- function(treated, controls, outcome_wide, pre_years) {
  Y_pre_treated <- as.numeric(outcome_wide[treated, as.character(pre_years)])
  Y_pre_controls <- as.matrix(outcome_wide[controls, as.character(pre_years)])
  n_ctl <- nrow(Y_pre_controls)
  
  Dmat <- Y_pre_controls %*% t(Y_pre_controls) + 1e-5 * diag(n_ctl)
  dvec <- Y_pre_controls %*% Y_pre_treated
  Amat <- cbind(rep(1, n_ctl), diag(n_ctl))
  bvec <- c(1, rep(0, n_ctl))
  meq <- 1
  
  qp_sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  w <- qp_sol$solution
  
  Y_full_controls <- as.matrix(outcome_wide[controls, ])
  synthetic <- colSums(w * Y_full_controls)
  gap <- as.numeric(outcome_wide[treated, ]) - synthetic
  return(gap)
}

pre_years <- 1970:1987
all_states <- 1:nrow(outcome_wide)
treated_id <- california_id
control_states <- setdiff(all_states, treated_id)

cal_gap <- compute_synthetic_for_state(treated_id, control_states, outcome_wide, pre_years)

placebo_gaps <- matrix(NA, nrow = length(control_states), ncol = length(years))
for (i in seq_along(control_states)) {
  other_controls <- setdiff(control_states, control_states[i])
  placebo_gaps[i, ] <- compute_synthetic_for_state(control_states[i], other_controls, outcome_wide, pre_years)
  if (i %% 10 == 0) cat("Processed", i, "of", length(control_states), "\n")
}



post_indices <- which(years >= 1989)
cal_post_avg <- mean(cal_gap[post_indices])

placebo_post_avg <- apply(placebo_gaps[, post_indices, drop = FALSE], 1, mean)

p_value <- mean(abs(placebo_post_avg) >= abs(cal_post_avg))

cat("\n Results \n")
cat("California average post‑treatment gap (reduction):", round(cal_post_avg, 2), "packs per capita\n")
cat("p‑value (two‑sided):", round(p_value, 4), "\n")


placebo_df <- as.data.frame(t(placebo_gaps))
colnames(placebo_df) <- paste0("Placebo_", control_states)
placebo_df$year <- years

placebo_long <- pivot_longer(placebo_df, cols = starts_with("Placebo_"),
                             names_to = "state", values_to = "gap")

p3 <- ggplot() +
  geom_line(data = placebo_long, aes(x = year, y = gap, group = state),
            color = "gray", alpha = 0.5, linewidth = 0.3) +
  geom_line(data = plot_df, aes(x = year, y = gap),
            color = "red", linewidth = 1.2) +
  geom_vline(xintercept = 1988, linetype = "dashed") +
  labs(title = "Placebo Test: California vs Control States",
       x = "Year", y = "Gap (Actual - Synthetic)") +
  theme_minimal()
print(p3)
ggsave("placebo_test_all_states.png", width = 8, height = 5)

