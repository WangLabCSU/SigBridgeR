#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

enum AggMethod
{
    AGG_MAX,
    AGG_SUM,
    AGG_MEAN,
    AGG_MEDIAN,
    AGG_FIRST
};

static AggMethod parse_method(const std::string& method)
{
    if (method == "max")    return AGG_MAX;
    if (method == "sum")    return AGG_SUM;
    if (method == "mean")   return AGG_MEAN;
    if (method == "median") return AGG_MEDIAN;
    if (method == "first")  return AGG_FIRST;
    stop("Unsupported method: " + method);
    return AGG_MAX;
}

static std::string name_key(SEXP s)
{
    if (s == NA_STRING)
    {
        return std::string("\001NA");
    }
    return std::string("\002") + Rf_translateCharUTF8(s);
}

struct Groups
{
    std::vector< std::vector<int>> index;
    std::vector<int> first_pos;
};

static Groups make_groups(CharacterVector names)
{
    std::unordered_map<std::string, int> seen;
    Groups groups;
    const int n = names.size();
    for (int i = 0; i < n; ++i)
    {
        SEXP s = STRING_ELT(names, i);
        std::string key = name_key(s);
        auto it = seen.find(key);
        if (it == seen.end())
        {
            int gid = groups.index.size();
            seen[key] = gid;
            groups.index.push_back(std::vector<int>());
            groups.first_pos.push_back(i);
            groups.index[gid].push_back(i);
        }
        else
        {
            groups.index[it->second].push_back(i);
        }
    }
    return groups;
}

static double median_no_na(std::vector<double> &vals)
{
    const int n = vals.size();
    if (n == 0)
    {
        return NA_REAL;
    }
    std::sort(vals.begin(), vals.end());
    if (n % 2 == 1)
    {
        return vals[n / 2];
    }
    else
    {
        return (vals[n / 2 - 1] + vals[n / 2]) / 2.0;
    }
}

// [[Rcpp::export()]]
NumericMatrix aggregate_dup_cols_cpp(
    NumericMatrix x,
    CharacterVector col_names,
    std::string method
)
{
    const int nr = x.nrow();
    const int nc = x.ncol();
    if (col_names.size() != nc)
    {
        stop("Length of col_names must equal ncol(x).");
    }
    AggMethod m = parse_method(method);
    Groups groups = make_groups(col_names);
    const int ng = groups.index.size();
    NumericMatrix out(nr, ng);
    arma::mat X(REAL(x), nr, nc, false, true);
    arma::mat Y(REAL(out), nr, ng, false, true);
    for (int g = 0; g < ng; ++g)
    {
        const std::vector<int> &idx = groups.index[g];
        if (idx.size() == 1 || m == AGG_FIRST)
        {
            Y.col(g) = X.col(idx[0]);
            continue;
        }
        for (int i = 0; i < nr; ++i)
        {
            if (m == AGG_SUM)
            {
                double s = 0.0;
                for (int j : idx)
                {
                    double v = X(i, j);
                    if (!ISNAN(v))
                    {
                        s += v;
                    }
                }
                Y(i, g) = s;
            }
            else if (m == AGG_MEAN)
            {
                double s = 0.0;
                int n = 0;
                for (int j : idx)
                {
                    double v = X(i, j);
                    if (!ISNAN(v))
                    {
                        s += v;
                        ++n;
                    }
                }
                Y(i, g) = n == 0 ? R_NaN : s / n;
            }
            else if (m == AGG_MAX)
            {
                double best = R_NegInf;
                bool has_value = false;
                for (int j : idx)
                {
                    double v = X(i, j);
                    if (!ISNAN(v))
                    {
                        if (!has_value || v > best)
                        {
                            best = v;
                        }
                        has_value = true;
                    }
                }
                Y(i, g) = has_value ? best : NA_REAL;
            }
            else if (m == AGG_MEDIAN)
            {
                std::vector<double> vals;
                vals.reserve(idx.size());
                for (int j : idx)
                {
                    double v = X(i, j);
                    if (!ISNAN(v))
                    {
                        vals.push_back(v);
                    }
                }
                Y(i, g) = median_no_na(vals);
            }
        }
    }
    return out;
}

// [[Rcpp::export()]]
NumericMatrix aggregate_dup_rows_cpp(
    NumericMatrix x,
    CharacterVector row_names,
    std::string method
)
{
    const int nr = x.nrow();
    const int nc = x.ncol();
    if (row_names.size() != nr)
    {
        stop("Length of row_names must equal nrow(x).");
    }
    AggMethod m = parse_method(method);
    Groups groups = make_groups(row_names);
    const int ng = groups.index.size();
    NumericMatrix out(ng, nc);
    arma::mat X(REAL(x), nr, nc, false, true);
    arma::mat Y(REAL(out), ng, nc, false, true);
    for (int g = 0; g < ng; ++g)
    {
        const std::vector<int> &idx = groups.index[g];
        if (idx.size() == 1 || m == AGG_FIRST)
        {
            Y.row(g) = X.row(idx[0]);
            continue;
        }
        for (int j = 0; j < nc; ++j)
        {
            if (m == AGG_SUM)
            {
                double s = 0.0;
                for (int i : idx)
                {
                    double v = X(i, j);
                    if (!ISNAN(v))
                    {
                        s += v;
                    }
                }
                Y(g, j) = s;
            }
            else if (m == AGG_MEAN)
            {
                double s = 0.0;
                int n = 0;
                for (int i : idx)
                {
                    double v = X(i, j);
                    if (!ISNAN(v))
                    {
                        s += v;
                        ++n;
                    }
                }
                Y(g, j) = n == 0 ? R_NaN : s / n;
            }
            else if (m == AGG_MAX)
            {
                double best = R_NegInf;
                bool has_value = false;
                for (int i : idx)
                {
                    double v = X(i, j);
                    if (!ISNAN(v))
                    {
                        if (!has_value || v > best)
                        {
                            best = v;
                        }
                        has_value = true;
                    }
                }
                Y(g, j) = has_value ? best : NA_REAL;
            }
            else if (m == AGG_MEDIAN)
            {
                std::vector<double> vals;
                vals.reserve(idx.size());
                for (int i : idx)
                {
                    double v = X(i, j);
                    if (!ISNAN(v))
                    {
                        vals.push_back(v);
                    }
                }
                Y(g, j) = median_no_na(vals);
            }
        }
    }
    return out;
}