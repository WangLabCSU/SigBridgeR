#include <Rcpp.h>
#include <limits>

using namespace Rcpp;

// Min-max normalization
// [[Rcpp::export]]
NumericVector normalize_metric(NumericVector x, bool invert = false)
{
    int n = x.size();
    double sum = 0.0;
    double sumsq = 0.0;
    int count = 0;
    double xmin = R_PosInf;
    double xmax = R_NegInf;
    // single pass: compute sum, sumsq, min, max, ignoring NA/NaN
    for (int i = 0; i < n; ++i)
    {
        double xi = x[i];
        if (NumericVector::is_na(xi) || R_IsNaN(xi))
        {
            continue;
        }
        sum += xi;
        sumsq += xi * xi;
        count++;
        if (xi < xmin) xmin = xi;
        if (xi > xmax) xmax = xi;
    }
    NumericVector out(n);
    // corresponds to stats::sd(x, na.rm = TRUE)
    // when the number of non-NA values <= 1, sd is NA; treat as non-normalizable -> 0.5
    if (count <= 1)
    {
        std::fill(out.begin(), out.end(), 0.5);
        return out;
    }
    double mean = sum / count;
    double variance = (sumsq - count * mean * mean) / (count - 1);
    // prevent tiny negative values caused by floating-point error
    if (variance < 0.0 && variance > -1e-15)
    {
        variance = 0.0;
    }
    double sd1 = std::sqrt(variance);
    if (sd1 <  std::numeric_limits<double>::epsilon())
    {
        std::fill(out.begin(), out.end(), 0.5);
        return out;
    }
    double range = xmax - xmin;
    for (int i = 0; i < n; ++i)
    {
        double xi = x[i];
        if (NumericVector::is_na(xi) || R_IsNaN(xi))
        {
            out[i] = NA_REAL;
        }
        else
        {
            double val = (xi - xmin) / range;
            out[i] = invert ? 1.0 - val : val;
        }
    }
    return out;
}

// -----------------------------------------------------------------------------
