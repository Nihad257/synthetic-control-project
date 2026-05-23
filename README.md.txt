# Synthetic Control from Scratch – California Smoking Ban

## Project Overview
Estimates the causal effect of California's Proposition 99 (1988 tobacco tax and control program) on cigarette sales using the **Synthetic Control Method** implemented entirely from scratch in R.

## Methods
- Quadratic programming with non‑negativity constraints and sum‑to‑one weights
- Ridge stabilisation to handle near‑singular matrices
- Permutation‑based inference (in‑space placebo tests)

## Key Results
- **Estimated reduction** : 19.7 packs per capita per year (95% credible via placebo distribution)
- **p‑value** : 0.0263 (two‑sided)
- **Interpretation** : Proposition 99 caused a statistically significant decline in cigarette sales.

## Visualisations
### Actual vs Synthetic California
![Actual vs Synthetic](california_actual_vs_synthetic.png)

### Treatment Effect (Gap)
![Gap](california_treatment_gap.png)

### Placebo Test (all control states in grey, California in red)
![Placebo](placebo_test_all_states.png)

## Repository Structure
- `01_simulate_and_optimize.R` – simulated data validation


## How to Run
1. Open `02_real_data_placebo.R` in RStudio.
2. Run line by line.
3. Plots and summary will be saved automatically.

## References
Abadie, Diamond, Hainmueller (2010). "Synthetic Control Methods for Comparative Case Studies".