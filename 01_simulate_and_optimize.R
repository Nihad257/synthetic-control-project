library(tidyverse)
library(ggplot2)
library(quadprog)

set.seed(123)

n_units <- 21          
n_time <- 30          
treatment_time <- 16  
true_effect <- 5       

unit_intercept <- rnorm(n_units, mean = 10, sd = 2)
unit_slope <- rnorm(n_units, mean = 0.5, sd = 0.1)
noise_sd <- 0.5        

sim_data <- expand.grid(unit = 1:n_units, time = 1:n_time) %>%
  mutate(
    unit_intercept = unit_intercept[unit],
    unit_slope = unit_slope[unit],
    noise = rnorm(n(), 0, noise_sd),
    outcome = unit_intercept + unit_slope * time + noise
  )

sim_data <- sim_data %>%
  mutate(
    treated = ifelse(unit == 1 & time >= treatment_time, 1, 0),
    outcome = outcome + true_effect * treated
  )

pre_time <- 1:(treatment_time - 1)

outcome_matrix <- sim_data %>%
  select(unit, time, outcome) %>%
  pivot_wider(names_from = time, values_from = outcome) %>%
  select(-unit) %>%
  as.matrix()

Y_pre_treated <- outcome_matrix[1, pre_time]
Y_pre_controls <- outcome_matrix[2:n_units, pre_time]

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
names(w_opt) <- paste0("control_", 2:n_units)

cat("Sum of weights:", sum(w_opt), "\n")
print(round(w_opt, 4))

Y_full_controls <- outcome_matrix[2:n_units, ]
synthetic_outcome <- colSums(w_opt * Y_full_controls)
actual_outcome <- outcome_matrix[1, ]

plot_df <- data.frame(
  time = 1:n_time,
  actual = actual_outcome,
  synthetic = synthetic_outcome,
  gap = actual_outcome - synthetic_outcome
)

ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = actual, color = "Actual")) +
  geom_line(aes(y = synthetic, color = "Synthetic")) +
  geom_vline(xintercept = treatment_time, linetype = "dashed", color = "red") +
  labs(title = "Synthetic Control: Actual vs Synthetic (Ridge Penalty)",
       x = "Time", y = "Outcome", color = "Series") +
  theme_minimal()

ggplot(plot_df, aes(x = time, y = gap)) +
  geom_line(color = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = treatment_time, linetype = "dashed", color = "red") +
  labs(title = "Estimated Treatment Effect (Gap)",
       x = "Time", y = "Effect (Actual - Synthetic)") +
  theme_minimal()

post_indices <- treatment_time:n_time
avg_effect <- mean(plot_df$gap[post_indices])
pre_indices <- 1:(treatment_time - 1)
pre_mse <- mean(plot_df$gap[pre_indices]^2)

cat("True effect:", true_effect, "\n")
cat("Estimated average post-treatment effect:", round(avg_effect, 3), "\n")
cat("Pre-treatment MSE (should be < 0.1):", round(pre_mse, 4), "\n")
