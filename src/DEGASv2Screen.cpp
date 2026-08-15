#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::CharacterVector get_bulk_hvg(
    const arma::mat& patdata,
    const Rcpp::CharacterVector& common_genes,
    const int n_hvg
)
{
    const arma::uword n_genes = patdata.n_rows;
    if (common_genes.size() != static_cast<R_xlen_t>(n_genes))
    {
        Rcpp::stop("The length of common_genes must equal the number of rows of patdata");
    }
    if (n_hvg <= 0 || n_genes == 0)
    {
        return Rcpp::CharacterVector(0);
    }
    const arma::uword k =
        std::min<arma::uword>(static_cast<arma::uword>(n_hvg), n_genes);
    arma::vec gene_sd(n_genes);
    gene_sd.fill(NA_REAL);
    for (arma::uword i = 0; i < n_genes; ++i)
    {
        // convert to a column vector to avoid rowvec layout errors
        arma::vec x = patdata.row(i).t();
        arma::uvec valid = arma::find_finite(x);
        // R::sd(x, na.rm=TRUE): returns NA when fewer than 2 valid values
        if (valid.n_elem >= 2)
        {
            arma::vec x_valid = x.elem(valid);
            gene_sd[i] = arma::stddev(x_valid, 0);
            // norm_type = 0 means sample standard deviation, denominator is n - 1
        }
    }
    // sort only finite standard deviations, NA last
    arma::uvec finite_idx = arma::find_finite(gene_sd);
    arma::uvec na_idx = arma::find_nonfinite(gene_sd);
    arma::uvec sorted_finite;
    if (finite_idx.n_elem > 0)
    {
        sorted_finite = finite_idx(
                            arma::sort_index(gene_sd.elem(finite_idx), "descend")
                        );
    }
    arma::uvec order_idx(sorted_finite.n_elem + na_idx.n_elem);
    if (sorted_finite.n_elem > 0)
    {
        order_idx.head(sorted_finite.n_elem) = sorted_finite;
    }
    if (na_idx.n_elem > 0)
    {
        order_idx.tail(na_idx.n_elem) = na_idx;
    }
    Rcpp::CharacterVector result(k);
    for (arma::uword i = 0; i < k; ++i)
    {
        result[i] = common_genes[order_idx[i]];
    }
    return result;
}


namespace
{

    constexpr double EPS = 1e-3;

// Both R's NA_real_ and NaN can be detected with std::isnan()
    inline bool is_missing(const double x)
    {
        return std::isnan(x);
    }


    /*
     * Normalizes a vector:
     *
     *     (x - mean(x, na.rm = TRUE)) /
     *     (sd(x, na.rm = TRUE) + 1e-3)
     *
     * Input:
     *   x       : input vector pointer
     *   n       : vector length
     *   out     : output vector pointer
     *
     * Notes:
     *   - does not create a temporary arma::vec
     *   - reads and writes directly through pointers
     *   - preserves NA/NaN
     */
    inline void norm_inplace_ptr(
        const double *x,
        double *out,
        const arma::uword n
    )
    {
        double sum = 0.0;
        arma::uword n_valid = 0;
        // first pass: compute the mean
        for (arma::uword i = 0; i < n; ++i)
        {
            const double value = x[i];
            if (!is_missing(value))
            {
                sum += value;
                ++n_valid;
            }
        }
        // when there are no valid values, R's mean(..., na.rm = TRUE) is NaN
        if (n_valid == 0)
        {
            for (arma::uword i = 0; i < n; ++i)
            {
                out[i] = NA_REAL;
            }
            return;
        }
        const double mean_value =
            sum / static_cast<double>(n_valid);
        // second pass: compute the sample variance
        double ss = 0.0;
        for (arma::uword i = 0; i < n; ++i)
        {
            const double value = x[i];
            if (!is_missing(value))
            {
                const double d = value - mean_value;
                ss += d * d;
            }
        }
        /*
         * R's sd():
         * - sample standard deviation when n_valid >= 2
         * - NA when n_valid < 2
         */
        double denominator;
        if (n_valid < 2)
        {
            denominator = NA_REAL;
        }
        else
        {
            const double sd_value =
                std::sqrt(ss / static_cast<double>(n_valid - 1));
            denominator = sd_value + EPS;
        }
        // third pass: write out the standardized result
        for (arma::uword i = 0; i < n; ++i)
        {
            const double value = x[i];
            if (is_missing(value))
            {
                out[i] = value;
            }
            else
            {
                out[i] = (value - mean_value) / denominator;
            }
        }
    }


    /*
     * Scales an already-normalized vector:
     *
     *     (x - min(x, na.rm = TRUE)) /
     *     (max(x, na.rm = TRUE) - min(x, na.rm = TRUE) + 1e-3)
     *
     * Input:
     *   x       : input vector pointer
     *   n       : vector length
     *   out     : output vector pointer
     *
     * Note:
     *   does not modify x, writes the result directly into out.
     */
    inline void scale_ptr(
        const double *x,
        double *out,
        const arma::uword n
    )
    {
        double xmin = std::numeric_limits<double>::infinity();
        double xmax = -std::numeric_limits<double>::infinity();
        arma::uword n_valid = 0;
        // first pass: compute the minimum and maximum
        for (arma::uword i = 0; i < n; ++i)
        {
            const double value = x[i];
            if (!is_missing(value))
            {
                if (value < xmin)
                {
                    xmin = value;
                }
                if (value > xmax)
                {
                    xmax = value;
                }
                ++n_valid;
            }
        }
        // when all values are NA/NaN, the result is all NA
        if (n_valid == 0)
        {
            for (arma::uword i = 0; i < n; ++i)
            {
                out[i] = x[i];
            }
            return;
        }
        const double denominator =
            xmax - xmin + EPS;
        // second pass: perform 0-1 scaling
        for (arma::uword i = 0; i < n; ++i)
        {
            const double value = x[i];
            if (is_missing(value))
            {
                out[i] = value;
            }
            else
            {
                out[i] = (value - xmin) / denominator;
            }
        }
    }

} // anonymous namespace


// [[Rcpp::export]]
arma::vec normFunc_ptr(const arma::vec& x)
{
    arma::vec out(x.n_elem);
    norm_inplace_ptr(
        x.memptr(),
        out.memptr(),
        x.n_elem
    );
    return out;
}


// [[Rcpp::export]]
arma::vec scaleFunc_ptr(const arma::vec& x)
{
    arma::vec out(x.n_elem);
    scale_ptr(
        x.memptr(),
        out.memptr(),
        x.n_elem
    );
    return out;
}


/*
 * Identical to the original R normalizeScale:
 *
 * normalizeScale <- function(X) {
 *   t(
 *     apply(
 *       t(apply(as.matrix(t(X)), 1, normFunc)),
 *       1,
 *       scaleFunc
 *     )
 *   )
 * }
 *
 * For each column of X:
 *
 *   1. normFunc
 *   2. scaleFunc
 *
 * Final output dimension:
 *
 *   ncol(X) x nrow(X)
 */
// [[Rcpp::export]]
arma::mat normalizeScale_ptr(const arma::mat& X)
{
    const arma::uword nrow = X.n_rows;
    const arma::uword ncol = X.n_cols;
    // the original R code finally returns ncol x nrow
    arma::mat result(ncol, nrow);
    // temporary buffer, allocated only once
    arma::vec normalized(nrow);
    const double *x_ptr = X.memptr();
    double *result_ptr = result.memptr();
    /*
     * Armadillo stores matrices in column-major order by default.
     *
     * Starting address of column j of X:
     *     x_ptr + j * nrow
     *
     * Row j of result:
     *     result(j, i)
     *
     * Since result is stored column-major, the elements of row j of result
     * are not contiguous in memory, so a pointer stride of ncol is used.
     */
    for (arma::uword j = 0; j < ncol; ++j)
    {
        const double *column_ptr =
            x_ptr + j * nrow;
        // step 1: column-wise normalization
        norm_inplace_ptr(
            column_ptr,
            normalized.memptr(),
            nrow
        );
        // step 2: scale the normalized column
        double xmin = std::numeric_limits<double>::infinity();
        double xmax = -std::numeric_limits<double>::infinity();
        arma::uword n_valid = 0;
        const double *norm_ptr = normalized.memptr();
        for (arma::uword i = 0; i < nrow; ++i)
        {
            const double value = norm_ptr[i];
            if (!is_missing(value))
            {
                if (value < xmin)
                {
                    xmin = value;
                }
                if (value > xmax)
                {
                    xmax = value;
                }
                ++n_valid;
            }
        }
        double *output_row_start =
            result_ptr + j;
        if (n_valid == 0)
        {
            for (arma::uword i = 0; i < nrow; ++i)
            {
                output_row_start[i * ncol] = norm_ptr[i];
            }
        }
        else
        {
            const double denominator =
                xmax - xmin + EPS;
            for (arma::uword i = 0; i < nrow; ++i)
            {
                const double value = norm_ptr[i];
                if (is_missing(value))
                {
                    output_row_start[i * ncol] = value;
                }
                else
                {
                    output_row_start[i * ncol] =
                        (value - xmin) / denominator;
                }
            }
        }
    }
    return result;
}


/*
 * Optimized version of preprocessCounts.
 *
 * Original operation:
 *
 *   normalizeScale(1.5^log2(X + 1))
 *
 * Instead of constructing a transformed matrix,
 * each element is computed directly on the fly:
 *
 *   transformed_ij = 1.5^log2(X_ij + 1)
 *
 * then normFunc and scaleFunc are applied directly.
 */
// [[Rcpp::export]]
arma::mat preprocessCounts_ptr(const arma::mat& X)
{
    const arma::uword nrow = X.n_rows;
    const arma::uword ncol = X.n_cols;
    // the original R function finally outputs ncol x nrow
    arma::mat result(ncol, nrow);
    // buffer allocated only once
    arma::vec transformed(nrow);
    arma::vec normalized(nrow);
    const double *x_ptr = X.memptr();
    double *result_ptr = result.memptr();
    for (arma::uword j = 0; j < ncol; ++j)
    {
        const double *column_ptr =
            x_ptr + j * nrow;
        /*
         * For the current column:
         *
         *   1.5^log2(x + 1)
         *
         * The result is written into transformed
         */
        for (arma::uword i = 0; i < nrow; ++i)
        {
            const double value = column_ptr[i];
            if (is_missing(value))
            {
                transformed[i] = value;
            }
            else
            {
                transformed[i] =
                    std::pow(1.5, std::log2(value + 1.0));
            }
        }
        /*
         * Phase 1: normFunc
         */
        norm_inplace_ptr(
            transformed.memptr(),
            normalized.memptr(),
            nrow
        );
        /*
         * Phase 2: scaleFunc
         */
        double xmin = std::numeric_limits<double>::infinity();
        double xmax = -std::numeric_limits<double>::infinity();
        arma::uword n_valid = 0;
        const double *norm_ptr = normalized.memptr();
        for (arma::uword i = 0; i < nrow; ++i)
        {
            const double value = norm_ptr[i];
            if (!is_missing(value))
            {
                if (value < xmin)
                {
                    xmin = value;
                }
                if (value > xmax)
                {
                    xmax = value;
                }
                ++n_valid;
            }
        }
        double *output_row_start =
            result_ptr + j;
        if (n_valid == 0)
        {
            for (arma::uword i = 0; i < nrow; ++i)
            {
                output_row_start[i * ncol] = norm_ptr[i];
            }
        }
        else
        {
            const double denominator =
                xmax - xmin + EPS;
            for (arma::uword i = 0; i < nrow; ++i)
            {
                const double value = norm_ptr[i];
                if (is_missing(value))
                {
                    output_row_start[i * ncol] = value;
                }
                else
                {
                    output_row_start[i * ncol] =
                        (value - xmin) / denominator;
                }
            }
        }
    }
    return result;
}
