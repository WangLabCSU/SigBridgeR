#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
bool IsCountsMatrixImpl(const arma::mat& x,
                        bool verbose = false,
                        double integer_tol = 1e-8,
                        double min_integer_fraction = 0.95)
{
    const arma::uword n = x.n_elem;
    if (n == 0)
    {
        if (verbose)
            Rcpp::Rcout << "Result: FALSE; matrix is empty.\n";
        return false;
    }
    // Check NA, NaN, and Inf
    if (!x.is_finite())
    {
        if (verbose)
            Rcpp::Rcout
                    << "Result: FALSE; matrix contains NA, NaN, or Inf.\n";
        return false;
    }
    // Check non-negativity
    double nonnegative_fraction =
        static_cast<double>(arma::accu(x >= 0)) / n;
    if (nonnegative_fraction < 1.0)
    {
        if (verbose)
            Rcpp::Rcout
                    << "Result: FALSE; matrix contains negative values.\n";
        return false;
    }
    // Check whether values are integers
    arma::mat rounded = arma::round(x);
    double integer_fraction =
        static_cast<double>(
            arma::accu(arma::abs(x - rounded) <= integer_tol)
        ) / n;
    if (integer_fraction < min_integer_fraction)
    {
        if (verbose)
        {
            Rcpp::Rcout
                    << "Integer fraction: " << integer_fraction << "\n"
                    << "Result: FALSE; insufficient integer-like values.\n";
        }
        return false;
    }
    // Check library sizes for all samples
    arma::rowvec library_sizes = arma::sum(x, 0);
    if (arma::any(library_sizes <= 0))
    {
        if (verbose)
            Rcpp::Rcout
                    << "Result: FALSE; at least one sample has a library size <= 0.\n";
        return false;
    }
    if (verbose)
    {
        double zero_fraction =
            static_cast<double>(arma::accu(x == 0)) / n;
        Rcpp::Rcout
                << "Non-negative fraction: " << nonnegative_fraction << "\n"
                << "Integer fraction: " << integer_fraction << "\n"
                << "Zero fraction: " << zero_fraction << "\n"
                << "Library sizes: ";
        for (arma::uword j = 0; j < library_sizes.n_elem && j < 10; ++j)
        {
            Rcpp::Rcout << library_sizes[j];
            if (j + 1 < library_sizes.n_elem)
                Rcpp::Rcout << ", ";
            if (j == 9 && j + 1 < library_sizes.n_elem)
                Rcpp::Rcout << " ...";
        }
        Rcpp::Rcout
                << "\nResult: TRUE; matrix is likely a counts matrix.\n";
    }
    return true;
}