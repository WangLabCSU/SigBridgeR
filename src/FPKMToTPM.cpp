#include <RcppArmadillo.h>
#include <cmath>
#include <string>
#include <vector>

#include "SigBridgeR.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::NumericMatrix FPKMToTPM_impl(const Rcpp::NumericMatrix &fpkm,
                                   bool na_as_zero = true,
                                   bool verbose = true) {
  // copy once to avoid modifying the original matrix from R
  Rcpp::NumericMatrix out = Rcpp::clone(fpkm);

  const int nr = out.nrow();
  const int nc = out.ncol();

  std::vector<double> sample_sums(static_cast<std::size_t>(nc), 0.0);

  R_xlen_t na_count = 0;

  /*
   * Armadillo view of the external memory:
   * - does not copy out's data
   * - both R and Armadillo use column-major ordering
   */
  if (nr > 0 && nc > 0) {
    arma::mat x(out.begin(), static_cast<arma::uword>(nr),
                static_cast<arma::uword>(nc),
                false, // copy_aux_mem = false
                true   // strict = true
    );

    // first pass: handle NA and compute the FPKM total for each sample
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      double *col = x.colptr(j);
      double total = 0.0;

      for (arma::uword i = 0; i < x.n_rows; ++i) {
        const double value = col[i];

        // std::isnan detects both NA_REAL and NaN
        if (std::isnan(value)) {
          ++na_count;

          if (!na_as_zero) {
            Rcpp::stop("The input matrix contains missing values. "
                       "Set na_as_zero = TRUE to replace them with zero.");
          }

          col[i] = 0.0;
        } else {
          total += value;
        }
      }

      sample_sums[static_cast<std::size_t>(j)] = total;
    }

    if (verbose && na_count > 0) {
      cli_emit("info",
               "Replacing " +
                   std::to_string(static_cast<unsigned long long>(na_count)) +
                   " missing value(s) with zero");
    }

    // second pass: perform TPM normalization in place
    int zero_sum_samples = 0;

    for (arma::uword j = 0; j < x.n_cols; ++j) {
      const double total = sample_sums[static_cast<std::size_t>(j)];

      if (total > 0.0) {
        double *col = x.colptr(j);
        const double scale = 1e6 / total;

        for (arma::uword i = 0; i < x.n_rows; ++i) {
          col[i] *= scale;
        }
      } else if (total == 0.0) {
        ++zero_sum_samples;
      }
      // keep the original value when total < 0;
      // negative inputs are checked on the R side
    }

    if (zero_sum_samples > 0) {
      Rcpp::warning("%d sample(s) have a total FPKM of zero; "
                    "their TPM values will be set to zero.",
                    zero_sum_samples);
    }
  } else {
    // for a zero-row matrix, every sample total is 0
    if (nc > 0) {
      Rcpp::warning("%d sample(s) have a total FPKM of zero; "
                    "their TPM values will be set to zero.",
                    nc);
    }
  }

  if (verbose) {
    cli_emit("success", "TPM conversion completed: " + std::to_string(nr) +
                            " genes × " + std::to_string(nc) + " samples");
  }

  return out;
}