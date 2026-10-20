# Spectral and Subspace Analysis

Null spaces, complete solution families, complex eigenpairs, and rank diagnostics.

## Run

Requires R 4.4+ and packages MASS and jsonlite.

```r
install.packages(c("MASS", "jsonlite"))
```

```sh
Rscript scripts/run.R
Rscript tests/checks.R
```

Generated metrics and figures are saved in `results/`. Read [input provenance](data/README.md) and [method changes](docs/provenance.md).

## Executed results

Numerical checks verify residuals rather than relying on exact floating-point equality. Full measured diagnostics are in [metrics.json](results/metrics.json).
The underdetermined system is represented by the complete family `x = x_p + N z`. Complex eigenvectors and eigenvalues retain their phase; eigenvector rank and residuals matter, not merely the number of returned values.

## Reproducibility and scope

Seeds are fixed in source. Validation uses numerical invariants or input/split checks. This is a maintained reconstruction of the supplied work; methodological corrections and newly generated data mean some results intentionally differ from the original reports. Original PDFs are omitted to avoid publishing student identifiers.
