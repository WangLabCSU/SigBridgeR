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
    // 一次遍历：计算 sum、sumsq、min、max，忽略 NA/NaN
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
    // 对应 stats::sd(x, na.rm = TRUE)
    // 当非 NA 数量 <= 1 时，sd 为 NA；这里按无法归一化处理为 0.5
    if (count <= 1)
    {
        std::fill(out.begin(), out.end(), 0.5);
        return out;
    }
    double mean = sum / count;
    double variance = (sumsq - count * mean * mean) / (count - 1);
    // 防止浮点误差导致极小负数
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
