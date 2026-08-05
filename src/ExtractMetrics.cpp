// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace Rcpp;

struct RowStats
{
    arma::vec mean;
    arma::vec var;
    int nrow;
    int ncol;
};

struct DropoutStats
{
    arma::vec dropout;
    int nrow;
    int ncol;
};

inline bool is_na_or_nan(double x)
{
    return ISNAN(x);
}

// ------------------------------------------------------------
// Quantile, type = 7, same as stats::quantile default
// ------------------------------------------------------------
double quantile_type7(const arma::vec& x, double prob)
{
    std::vector<double> vals;
    vals.reserve(x.n_elem);
    for (arma::uword i = 0; i < x.n_elem; ++i)
    {
        double v = x[i];
        if (!is_na_or_nan(v))
        {
            vals.push_back(v);
        }
    }
    if (vals.empty()) return NA_REAL;
    std::sort(vals.begin(), vals.end());
    if (prob <= 0.0) return vals.front();
    if (prob >= 1.0) return vals.back();
    double h = (vals.size() - 1) * prob;
    std::size_t lo = static_cast<std::size_t>(std::floor(h));
    std::size_t hi = static_cast<std::size_t>(std::ceil(h));
    double gamma = h - lo;
    return vals[lo] * (1.0 - gamma) + vals[hi] * gamma;
}

// ------------------------------------------------------------
// Pearson correlation with complete finite observations
// ------------------------------------------------------------
double pearson_cor(const std::vector<double> &x,
                   const std::vector<double> &y)
{
    std::size_t n = x.size();
    if (n != y.size() || n < 2) return NA_REAL;
    double sx = 0.0, sy = 0.0;
    std::size_t m = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (R_FINITE(x[i]) && R_FINITE(y[i]))
        {
            sx += x[i];
            sy += y[i];
            ++m;
        }
    }
    if (m < 2) return NA_REAL;
    double mx = sx / m;
    double my = sy / m;
    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (R_FINITE(x[i]) && R_FINITE(y[i]))
        {
            double dx = x[i] - mx;
            double dy = y[i] - my;
            sxx += dx * dx;
            syy += dy * dy;
            sxy += dx * dy;
        }
    }
    if (sxx <= 0.0 || syy <= 0.0) return NA_REAL;
    return sxy / std::sqrt(sxx * syy);
}

// ------------------------------------------------------------
// Average rank for Spearman correlation
// ------------------------------------------------------------
std::vector<double> average_rank(const std::vector<double> &x)
{
    std::size_t n = x.size();
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b)
    {
        return x[a] < x[b];
    });
    std::vector<double> r(n);
    std::size_t start = 0;
    while (start < n)
    {
        std::size_t end = start + 1;
        while (end < n && x[idx[end]] == x[idx[start]])
        {
            ++end;
        }
        double avg_rank = 0.5 * (static_cast<double>(start + 1) +
                                 static_cast<double>(end));
        for (std::size_t k = start; k < end; ++k)
        {
            r[idx[k]] = avg_rank;
        }
        start = end;
    }
    return r;
}

// ------------------------------------------------------------
// Spearman correlation with complete finite observations
// ------------------------------------------------------------
double spearman_cor(const arma::vec& x, const arma::vec& y)
{
    if (x.n_elem != y.n_elem) return NA_REAL;
    std::vector<double> xx;
    std::vector<double> yy;
    xx.reserve(x.n_elem);
    yy.reserve(y.n_elem);
    for (arma::uword i = 0; i < x.n_elem; ++i)
    {
        double xi = x[i];
        double yi = y[i];
        if (!is_na_or_nan(xi) && !is_na_or_nan(yi) &&
                R_FINITE(xi) && R_FINITE(yi))
        {
            xx.push_back(xi);
            yy.push_back(yi);
        }
    }
    if (xx.size() < 2) return NA_REAL;
    std::vector<double> rx = average_rank(xx);
    std::vector<double> ry = average_rank(yy);
    return pearson_cor(rx, ry);
}

// ------------------------------------------------------------
// Row means and row variances for dgCMatrix
// ------------------------------------------------------------
RowStats row_stats_dgCMatrix(SEXP mat)
{
    S4 m(mat);
    IntegerVector dim = m.slot("Dim");
    IntegerVector p = m.slot("p");
    IntegerVector i = m.slot("i");
    NumericVector x = m.slot("x");
    int nr = dim[0];
    int nc = dim[1];
    arma::vec sum(nr, arma::fill::zeros);
    arma::vec sumsq(nr, arma::fill::zeros);
    std::vector<unsigned char> has_na(nr, 0);
    for (int col = 0; col < nc; ++col)
    {
        for (int k = p[col]; k < p[col + 1]; ++k)
        {
            int row = i[k];
            double val = x[k];
            if (is_na_or_nan(val))
            {
                has_na[row] = 1;
            }
            else
            {
                sum[row] += val;
                sumsq[row] += val * val;
            }
        }
    }
    arma::vec mean(nr);
    arma::vec var(nr);
    for (int r = 0; r < nr; ++r)
    {
        if (has_na[r] || nc == 0)
        {
            mean[r] = NA_REAL;
            var[r] = NA_REAL;
            continue;
        }
        mean[r] = sum[r] / static_cast<double>(nc);
        if (nc <= 1)
        {
            var[r] = NA_REAL;
        }
        else
        {
            double v = (sumsq[r] - sum[r] * sum[r] / nc) /
                       static_cast<double>(nc - 1);
            if (v < 0.0 && v > -1e-12) v = 0.0;
            var[r] = v;
        }
    }
    return {mean, var, nr, nc};
}

// ------------------------------------------------------------
// Row means and variances for dense matrix
// ------------------------------------------------------------
RowStats row_stats_dense(SEXP mat, bool s4_dense)
{
    IntegerVector dim;
    NumericVector x;
    if (s4_dense)
    {
        S4 m(mat);
        dim = m.slot("Dim");
        x = m.slot("x");
    }
    else
    {
        dim = Rf_getAttrib(mat, R_DimSymbol);
        x = as<NumericVector>(mat);
    }
    int nr = dim[0];
    int nc = dim[1];
    arma::vec mean(nr);
    arma::vec var(nr);
    for (int r = 0; r < nr; ++r)
    {
        bool has_na = false;
        double s = 0.0;
        double ss = 0.0;
        for (int c = 0; c < nc; ++c)
        {
            double val = x[r + c * nr];
            if (is_na_or_nan(val))
            {
                has_na = true;
                break;
            }
            s += val;
            ss += val * val;
        }
        if (has_na || nc == 0)
        {
            mean[r] = NA_REAL;
            var[r] = NA_REAL;
            continue;
        }
        mean[r] = s / static_cast<double>(nc);
        if (nc <= 1)
        {
            var[r] = NA_REAL;
        }
        else
        {
            double v = (ss - s * s / nc) / static_cast<double>(nc - 1);
            if (v < 0.0 && v > -1e-12) v = 0.0;
            var[r] = v;
        }
    }
    return {mean, var, nr, nc};
}

// ------------------------------------------------------------
// Dispatcher for row stats
// ------------------------------------------------------------
RowStats row_stats(SEXP mat)
{
    if (Rf_isS4(mat) && Rf_inherits(mat, "dgCMatrix"))
    {
        return row_stats_dgCMatrix(mat);
    }
    if (Rf_isS4(mat) && Rf_inherits(mat, "dgeMatrix"))
    {
        return row_stats_dense(mat, true);
    }
    if (Rf_isMatrix(mat))
    {
        return row_stats_dense(mat, false);
    }
    stop("Unsupported assay_data matrix type. Supported: dgCMatrix, dgeMatrix, base matrix.");
}

// ------------------------------------------------------------
// Dropout rate for dgCMatrix
// rowMeans(counts == 0)
// ------------------------------------------------------------
DropoutStats dropout_dgCMatrix(SEXP mat)
{
    S4 m(mat);
    IntegerVector dim = m.slot("Dim");
    IntegerVector p = m.slot("p");
    IntegerVector i = m.slot("i");
    NumericVector x = m.slot("x");
    int nr = dim[0];
    int nc = dim[1];
    std::vector<int> detected(nr, 0);
    std::vector<unsigned char> has_na(nr, 0);
    for (int col = 0; col < nc; ++col)
    {
        for (int k = p[col]; k < p[col + 1]; ++k)
        {
            int row = i[k];
            double val = x[k];
            if (is_na_or_nan(val))
            {
                has_na[row] = 1;
            }
            else if (val != 0.0)
            {
                detected[row]++;
            }
        }
    }
    arma::vec dropout(nr);
    for (int r = 0; r < nr; ++r)
    {
        if (has_na[r] || nc == 0)
        {
            dropout[r] = NA_REAL;
        }
        else
        {
            dropout[r] = 1.0 - static_cast<double>(detected[r]) /
                         static_cast<double>(nc);
        }
    }
    return {dropout, nr, nc};
}

// ------------------------------------------------------------
// Dropout for dense matrix
// ------------------------------------------------------------
DropoutStats dropout_dense(SEXP mat, bool s4_dense)
{
    IntegerVector dim;
    NumericVector x;
    if (s4_dense)
    {
        S4 m(mat);
        dim = m.slot("Dim");
        x = m.slot("x");
    }
    else
    {
        dim = Rf_getAttrib(mat, R_DimSymbol);
        x = as<NumericVector>(mat);
    }
    int nr = dim[0];
    int nc = dim[1];
    arma::vec dropout(nr);
    for (int r = 0; r < nr; ++r)
    {
        bool has_na = false;
        int zeros = 0;
        for (int c = 0; c < nc; ++c)
        {
            double val = x[r + c * nr];
            if (is_na_or_nan(val))
            {
                has_na = true;
                break;
            }
            if (val == 0.0) zeros++;
        }
        if (has_na || nc == 0)
        {
            dropout[r] = NA_REAL;
        }
        else
        {
            dropout[r] = static_cast<double>(zeros) / static_cast<double>(nc);
        }
    }
    return {dropout, nr, nc};
}

// ------------------------------------------------------------
// Dispatcher for dropout
// ------------------------------------------------------------
DropoutStats dropout_stats(SEXP mat)
{
    if (Rf_isS4(mat) && Rf_inherits(mat, "dgCMatrix"))
    {
        return dropout_dgCMatrix(mat);
    }
    if (Rf_isS4(mat) && Rf_inherits(mat, "dgeMatrix"))
    {
        return dropout_dense(mat, true);
    }
    if (Rf_isMatrix(mat))
    {
        return dropout_dense(mat, false);
    }
    stop("Unsupported counts matrix type. Supported: dgCMatrix, dgeMatrix, base matrix.");
}

// ------------------------------------------------------------
// Main exported function
// ------------------------------------------------------------

// [[Rcpp::export]]
List ExtractMetricsCpp(SEXP assay_data,
                       SEXP counts,
                       double low_expressed_thresh = 0.2)
{
    RowStats rs = row_stats(assay_data);
    DropoutStats ds = dropout_stats(counts);
    if (rs.nrow != ds.nrow)
    {
        stop("assay_data and counts must have the same number of rows.");
    }
    if (rs.ncol != ds.ncol)
    {
        stop("assay_data and counts must have the same number of columns.");
    }
    // ----------------------------------------------------------
    // 1. Variance-mean correlation
    // ----------------------------------------------------------
    double q = quantile_type7(rs.mean, low_expressed_thresh);
    double variance_mean_cor = NA_REAL;
    if (!is_na_or_nan(q))
    {
        std::vector<double> lx;
        std::vector<double> ly;
        int n_valid_basic = 0;
        for (int g = 0; g < rs.nrow; ++g)
        {
            double mu = rs.mean[g];
            double va = rs.var[g];
            if (!is_na_or_nan(mu) &&
                    !is_na_or_nan(va) &&
                    mu > q &&
                    va > 0.0)
            {
                n_valid_basic++;
                double x = std::log10(mu + 1.0);
                double y = std::log10(va + 1e-6);
                lx.push_back(x);
                ly.push_back(y);
            }
        }
        if (n_valid_basic > 100)
        {
            variance_mean_cor = pearson_cor(lx, ly);
        }
    }
    // ----------------------------------------------------------
    // 3. Dropout residual analysis
    // ----------------------------------------------------------
    double dropout_residual = spearman_cor(ds.dropout, rs.mean);
    double mean_dropout_residual = NA_REAL;
    if (!is_na_or_nan(dropout_residual))
    {
        mean_dropout_residual = std::fabs(dropout_residual);
    }
    return List::create(
               _["variance_mean_cor"] = variance_mean_cor,
               _["mean_dropout_residual"] = mean_dropout_residual,
               _["n_cells"] = rs.ncol,
               _["n_genes"] = rs.nrow
           );
}