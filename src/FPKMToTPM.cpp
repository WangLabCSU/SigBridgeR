#include <Rcpp.h>
#include "Rtatami.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "SigBridgeR.h"

// [[Rcpp::depends(beachmat, assorthead)]]
// [[Rcpp::plugins(cpp17)]]

using namespace Rcpp;

/*
 * FPKM to TPM conversion.
 *
 * The initialized_fpkm passed from R must be the external pointer produced by
 *
 *     beachmat::initializeCpp(fpkm)
 *
 * The beachmat backend is only used for reading, so the output is always a
 * new dense NumericMatrix.
 */
// [[Rcpp::export]]
NumericMatrix fpkm_to_tpm_cpp(SEXP initialized_fpkm, bool na_as_zero = true,
                              bool verbose = true) {
  /*
   * Parse the beachmat/tatami external pointer.
   */
  Rtatami::BoundNumericPointer parsed(initialized_fpkm);
  auto matrix = parsed->ptr;

  const std::size_t nr = static_cast<std::size_t>(matrix->nrow());

  const std::size_t nc = static_cast<std::size_t>(matrix->ncol());

  if (nr > static_cast<std::size_t>(INT_MAX) ||
      nc > static_cast<std::size_t>(INT_MAX)) {
    stop("Matrix dimensions exceed R integer limits.");
  }

  /*
   * Output matrix.
   *
   * beachmat handles reading the input matrix; the output is created by Rcpp.
   */
  NumericMatrix output(static_cast<int>(nr), static_cast<int>(nc));

  /*
   * FPKM total for each sample.
   */
  std::vector<double> sample_sums(nc, 0.0);

  R_xlen_t na_count = 0;

  /*
   * Access the matrix column by column.
   *
   * For sparse matrices, unstored positions are returned as 0;
   * for DelayedMatrix, fetch() reads or computes the corresponding column.
   */
  auto accessor = matrix->dense_column();

  std::vector<double> buffer(nr);

  /*
   * First pass:
   *
   * - check for NA/NaN;
   * - decide whether to treat them as 0 according to na_as_zero;
   * - compute the FPKM total for each sample.
   *
   * Because the beachmat backend is usually read-only, values cannot be
   * modified in place, so this pass only records the totals and the second
   * pass writes into output.
   */
  for (std::size_t col = 0; col < nc; ++col) {
    auto values = accessor->fetch(col, buffer.data());

    double total = 0.0;

    for (std::size_t row = 0; row < nr; ++row) {
      const double value = values[row];

      /*
       * std::isnan detects both NA_real_ and NaN.
       */
      if (std::isnan(value)) {
        ++na_count;

        if (!na_as_zero) {
          stop("The input matrix contains missing values. "
               "Set na_as_zero = TRUE to replace them with zero.");
        }

        /*
         * When na_as_zero = TRUE, NA/NaN is not added to total.
         */
        continue;
      }

      total += value;
    }

    sample_sums[col] = total;
  }

  /*
   * Emit the NA replacement message.
   *
   * cli_emit() must be called on the R main thread.
   */
  if (verbose && na_count > 0) {
    cli_emit("info",
             "Replacing " +
                 std::to_string(static_cast<unsigned long long>(na_count)) +
                 " missing value(s) with zero");
  }

  /*
   * Second pass:
   *
   * TPM_i = FPKM_i / sum(FPKM) * 1e6
   */
  int zero_sum_samples = 0;

  for (std::size_t col = 0; col < nc; ++col) {
    const double total = sample_sums[col];

    if (!std::isfinite(total)) {
      stop("Input contains a non-finite sample total.");
    }

    auto values = accessor->fetch(col, buffer.data());

    if (total > 0.0) {
      const double scale = 1e6 / total;

      for (std::size_t row = 0; row < nr; ++row) {
        const double value = values[row];

        if (std::isnan(value)) {
          /*
           * When na_as_zero = TRUE, NA/NaN is written out as 0.
           * When na_as_zero = FALSE, the first pass has already stopped.
           */
          output(static_cast<int>(row), static_cast<int>(col)) = 0.0;
        } else {
          output(static_cast<int>(row), static_cast<int>(col)) = value * scale;
        }
      }
    } else if (total == 0.0) {
      ++zero_sum_samples;

      /*
       * For a sample whose total is 0, all of its TPM values are 0.
       */
      for (std::size_t row = 0; row < nr; ++row) {
        output(static_cast<int>(row), static_cast<int>(col)) = 0.0;
      }
    } else {
      /*
       * Keep consistent with the original implementation:
       * no normalization is performed when total < 0.
       *
       * FPKM values are normally non-negative, so checking this on the R
       * side beforehand is recommended.
       */
      for (std::size_t row = 0; row < nr; ++row) {
        const double value = values[row];

        output(static_cast<int>(row), static_cast<int>(col)) =
            std::isnan(value) ? 0.0 : value;
      }
    }
  }

  if (zero_sum_samples > 0) {
    warning("%d sample(s) have a total FPKM of zero; "
            "their TPM values will be set to zero.",
            zero_sum_samples);
  }

  if (verbose) {
    cli_emit("success", "TPM conversion completed: " + std::to_string(nr) +
                            " genes x " + std::to_string(nc) + " samples");
  }

  return output;
}