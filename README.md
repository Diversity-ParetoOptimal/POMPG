# PO-MPG

PO-MPG extends the Pareto-optimal predictor weighting framework of the [`ParetoR`](https://github.com/Diversity-ParetoOptimal/ParetoR) package (Normal Boundary Intersection / NBI) from a single adverse-impact protected subgroup to multiple protected subgroups simultaneously.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Diversity-ParetoOptimal/ParetoR")
remotes::install_github("Diversity-ParetoOptimal/POMPG")
```

## Usage

```r
library(POMPG)

# Predictor inter-correlation matrix
Rx <- matrix(c(1,  .37, .51,
               .37,   1, .03,
               .51, .03,   1), 3, 3)

# Predictor-criterion validities
Rxy <- c(.32, .52, .22)

# Standardized subgroup mean differences (majority - minority),
# one row per protected subgroup
subd <- rbind(c(.39, .72, -.09),   # subgroup 1
              c(.20, .30, -.05))   # subgroup 2

# Proportion of each subgroup in the applicant pool
prop <- c(0.15, 0.10)

out <- POMA_Focal(
  sr = 0.15, prop = prop, Rx = Rx, Rxy = Rxy, subd = subd,
  Spac = 20, graph = FALSE, display_solution = FALSE
)

out$Pareto_Fmat  # composite validity (Ry) and adverse impact ratio (AIR1) per solution
out$Pareto_Xmat  # predictor weights per solution
```

## Functions

- `POMA_Focal()`: multi-subgroup Pareto-optimal predictor weighting via NBI.
- `get_weights()`: convenience wrapper dispatching to `ParetoR::ParetoR()` (single subgroup) or `POMA_Focal()` (multiple subgroups).
- `get_sample()`: draw a stratified sample from subgroup-specific multivariate normal distributions, for simulation studies.
- `hire()`: simulate hiring under a selection ratio and predictor weights, returning composite validity and subgroup adverse impact ratios.
- `calc_biserial()` / `calc_d()`: convert between standardized mean differences and biserial correlations.

## License

MIT
