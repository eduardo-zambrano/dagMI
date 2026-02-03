# dagMI: Multilinear Inequality Tests for DAG Structures

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

The `dagMI` package implements multilinear inequality tests for directed acyclic graph (DAG) structures based on Carbery's generalization of the Cauchy-Schwarz inequality. This methodology provides a principled framework for:

- **Testing DAG compatibility**: Determine whether observed data is consistent with a hypothesized causal structure
- **Distinguishing Markov-equivalent DAGs**: Unlike conditional independence tests, multilinear inequalities can discriminate between DAGs that encode the same conditional independence relations
- **Confirmatory causal analysis**: Obtain formal p-values for DAG hypotheses (not just rankings or scores)

## Associated Paper

This package accompanies the paper:

> Zambrano, E. (2026). "Testing Causal DAG Structures via Multilinear Inequalities."

The paper establishes the theoretical foundation for using density-corrected Carbery inequalities to test causal DAG structures, including:

- The Q_n^G functional that encodes DAG-specific distributional constraints
- The density correction factor p_i(x)^{1/(n+1)} required for probability measures
- Constrained bootstrap inference respecting DAG structure
- Power characterization and comparison with LiNGAM, ANM, and NOTEARS

## Installation

```r
# Install from GitHub (development version)
# install.packages("devtools")
devtools::install_github("ezambrano/dagMI")
```

## Quick Start

```r
library(dagMI)

# Define a chain DAG: X1 -> X2 -> X3
A <- matrix(c(0, 1, 0,
              0, 0, 1,
              0, 0, 0), 3, 3, byrow = TRUE)
chain <- dag(A, nodes = c("X1", "X2", "X3"))

# Generate data consistent with the chain
set.seed(123)
n <- 500
X1 <- rnorm(n)
X2 <- 0.7 * X1 + rnorm(n, sd = sqrt(1 - 0.7^2))
X3 <- 0.7 * X2 + rnorm(n, sd = sqrt(1 - 0.7^2))
data <- cbind(X1, X2, X3)

# Test DAG compatibility
result <- mi_test(data, chain, B = 500)
print(result)
```

## Key Features

### Core Functions

| Function | Description |
|----------|-------------|
| `mi_test()` | Main interface for multilinear inequality testing |
| `compute_qn()` | Compute the Q_n^G functional for a DAG |
| `compute_test_stat()` | Compute the test statistic T_h^G |
| `compare_dags()` | Compare two competing DAG structures |
| `subgraph_test()` | Test subgraphs for scalability |

### Test Function Selection

| Function | Description |
|----------|-------------|
| `select_test_function()` | Adaptive selection via data splitting (Algorithm 4) |
| `h_squared()`, `h_poly()`, `h_exp()`, `h_indicator()` | Built-in test functions |
| `compute_test_stat_multiple()` | Evaluate multiple test functions |

### Inference and Diagnostics

| Function | Description |
|----------|-------------|
| `constrained_bootstrap()` | Bootstrap inference respecting DAG structure |
| `power_mi_test()` | Power computation for sample size planning |
| `diagnose_confounding()` | Diagnostic strategies for suspected confounding |
| `fisher_combine()`, `bonferroni_orderings()` | Multiple testing corrections |

### Visualization

| Function | Description |
|----------|-------------|
| `plot_bootstrap()` | Bootstrap distribution plots |
| `plot_dag_comparison()` | Visual comparison of DAG tests |
| `plot_diagnostics()` | Diagnostic plots for test results |

## Methodology

### The Multilinear Inequality

For a DAG G with n nodes, the density-corrected Carbery inequality states:

$$E\left[\prod_{j=1}^n h_j(X_j) \cdot p_j(X_j)^{1/(n+1)}\right] \leq Q_n^{\mathcal{G}} \prod_{j=1}^n \left(E[h_j(X_j)^{n+1}]\right)^{1/(n+1)}$$

where:
- $h_j$ are non-negative test functions
- $p_j(x_j)$ are marginal densities
- $Q_n^{\mathcal{G}}$ is a DAG-specific functional encoding the causal structure

### Why It Works

The Q_n functional captures **more than conditional independence**:

1. **Marginal density structure**: Q_n involves marginal densities $p_j(x_j)$ raised to the power $1/(n+1)$
2. **Bivariate dependencies**: Q_n includes consecutive bivariate density products $\int p_{j-1}(x)p_j(x)dx$
3. **DAG-specific geometry**: Different DAGs impose different constraints on these quantities

This allows discrimination between Markov-equivalent DAGs that conditional independence tests cannot achieve.

### Test Function Guidance

| Distribution Type | Recommended Test Function |
|-------------------|---------------------------|
| Near-Gaussian | `h_squared` (x²) |
| Heavy-tailed | `h_poly` with k=2 (x⁴) |
| Skewed | `h_exp` with t=0.5 |
| Bounded support | `h_indicator` |

Use `select_test_function()` for automatic adaptive selection.

## Computational Limits

| DAG Size | Full Bootstrap | Point Estimate |
|----------|----------------|----------------|
| n ≤ 5 | Quadrature-based (fast) | Fast |
| 6 ≤ n ≤ 8 | Monte Carlo (moderate) | Moderate |
| 9 ≤ n ≤ 15 | Limited | Monte Carlo |
| n > 15 | Use subgraph testing | Markov blanket strategy |

For large DAGs, use `subgraph_test()` on Markov blankets or ancestral subgraphs.

## Documentation

- **Vignette**: `vignette("dagMI-intro")` - Comprehensive worked examples
- **Function help**: `?mi_test`, `?compute_qn`, etc.
- **Paper**: Full theoretical development and simulation studies

## Dependencies

**Required:**
- R (≥ 4.0.0)
- Rcpp, RcppArmadillo (C++ integration)
- stats, methods (base R)
- ggplot2 (visualization)

**Suggested:**
- future, furrr (parallel processing)
- igraph (DAG visualization)
- testthat (testing)

## Citation

If you use this package, please cite:

```bibtex
@article{zambrano2026dag,
  title={Testing Causal DAG Structures via Multilinear Inequalities},
  author={Zambrano, Eduardo},
  journal={Journal of the American Statistical Association},
  year={2026},
  note={Under revision}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contributing

Issues and pull requests welcome at [GitHub](https://github.com/ezambrano/dagMI).
