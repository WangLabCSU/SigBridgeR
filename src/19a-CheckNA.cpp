#include <Rcpp.h>
#include <R_ext/Arith.h>
#include <climits>

using namespace Rcpp;

inline bool is_one_na(SEXP x)
{
    if (Rf_isNull(x) || XLENGTH(x) != 1) return false;
    switch (TYPEOF(x))
    {
        case REALSXP:
            return ISNAN(REAL(x)[0]);
        case INTSXP:
            return INTEGER(x)[0] == NA_INTEGER;
        case LGLSXP:
            return LOGICAL(x)[0] == NA_LOGICAL;
        case STRSXP:
            return STRING_ELT(x, 0) == NA_STRING;
        case CPLXSXP:
            {
                Rcomplex z = COMPLEX(x)[0];
                return ISNAN(z.r) || ISNAN(z.i);
            }
        default:
            return false;
    }
}

template <typename Add>
void scan_atomic(SEXP x, Add add)
{
    R_xlen_t n = XLENGTH(x);
    switch (TYPEOF(x))
    {
        case REALSXP:
            {
                const double *p = REAL(x);
                for (R_xlen_t k = 0; k < n; ++k)
                {
                    if (ISNAN(p[k])) add(k);
                }
                break;
            }
        case INTSXP:
            {
                const int *p = INTEGER(x);
                for (R_xlen_t k = 0; k < n; ++k)
                {
                    if (p[k] == NA_INTEGER) add(k);
                }
                break;
            }
        case LGLSXP:
            {
                const int *p = LOGICAL(x);
                for (R_xlen_t k = 0; k < n; ++k)
                {
                    if (p[k] == NA_LOGICAL) add(k);
                }
                break;
            }
        case STRSXP:
            {
                for (R_xlen_t k = 0; k < n; ++k)
                {
                    if (STRING_ELT(x, k) == NA_STRING) add(k);
                }
                break;
            }
        case CPLXSXP:
            {
                const Rcomplex* p = COMPLEX(x);
                for (R_xlen_t k = 0; k < n; ++k)
                {
                    if (ISNAN(p[k].r) || ISNAN(p[k].i)) add(k);
                }
                break;
            }
        case RAWSXP:
            break;
        case VECSXP:
            {
                for (R_xlen_t k = 0; k < n; ++k)
                {
                    if (is_one_na(VECTOR_ELT(x, k))) add(k);
                }
                break;
            }
        default:
            stop("Unsupported type for NA scan: %s", Rf_type2char(TYPEOF(x)));
    }
}

inline bool is_na_at(SEXP x, R_xlen_t k)
{
    switch (TYPEOF(x))
    {
        case REALSXP:
            return ISNAN(REAL(x)[k]);
        case INTSXP:
            return INTEGER(x)[k] == NA_INTEGER;
        case LGLSXP:
            return LOGICAL(x)[k] == NA_LOGICAL;
        case STRSXP:
            return STRING_ELT(x, k) == NA_STRING;
        case CPLXSXP:
            {
                Rcomplex z = COMPLEX(x)[k];
                return ISNAN(z.r) || ISNAN(z.i);
            }
        default:
            stop("Unsupported sparse x slot type: %s", Rf_type2char(TYPEOF(x)));
    }
}

// [[Rcpp::export]]
List check_na_vector_cpp(SEXP x)
{
    std::vector<int> pos;
    double count = 0;
    scan_atomic(x, [&](R_xlen_t k)
    {
        if (k >= INT_MAX) stop("Long vector index exceeds integer range.");
        ++count;
        pos.push_back(static_cast<int>(k + 1));
    });
    return List::create(
               _["count"] = count,
               _["positions"] = wrap(pos)
           );
}

// [[Rcpp::export]]
List check_na_dense2d_cpp(SEXP x, IntegerVector dim)
{
    if (dim.size() < 2) stop("dim must have length >= 2.");
    int nr = dim[0];
    int nc = dim[1];
    R_xlen_t expected = static_cast<R_xlen_t>(nr) * static_cast<R_xlen_t>(nc);
    if (XLENGTH(x) != expected)
    {
        stop("Data length does not match supplied dim.");
    }
    std::vector<int> rows;
    std::vector<int> cols;
    double count = 0;
    scan_atomic(x, [&](R_xlen_t k)
    {
        ++count;
        rows.push_back(static_cast<int>(k % nr) + 1);
        cols.push_back(static_cast<int>(k / nr) + 1);
    });
    return List::create(
               _["count"] = count,
               _["row"] = wrap(rows),
               _["col"] = wrap(cols)
           );
}

// [[Rcpp::export]]
List check_na_dataframe_cpp(List df, int nr)
{
    int nc = df.size();
    std::vector<int> rows;
    std::vector<int> cols;
    double count = 0;
    for (int j = 0; j < nc; ++j)
    {
        SEXP col = df[j];
        if (XLENGTH(col) != nr)
        {
            stop("All data.frame columns must have length nrow(data).");
        }
        scan_atomic(col, [&](R_xlen_t r)
        {
            ++count;
            rows.push_back(static_cast<int>(r + 1));
            cols.push_back(j + 1);
        });
    }
    return List::create(
               _["count"] = count,
               _["row"] = wrap(rows),
               _["col"] = wrap(cols)
           );
}

// For dgCMatrix / lgCMatrix / sparse compressed-column matrices
// [[Rcpp::export]]
List check_na_sparse_csc_cpp(SEXP x, IntegerVector i, IntegerVector p,
                             IntegerVector dim)
{
    if (Rf_isNull(x))
    {
        return List::create(
                   _["count"] = 0.0,
                   _["row"] = IntegerVector(0),
                   _["col"] = IntegerVector(0)
               );
    }
    int nc = dim[1];
    std::vector<int> rows;
    std::vector<int> cols;
    double count = 0;
    for (int col = 0; col < nc; ++col)
    {
        for (int k = p[col]; k < p[col + 1]; ++k)
        {
            if (is_na_at(x, k))
            {
                ++count;
                rows.push_back(i[k] + 1);
                cols.push_back(col + 1);
            }
        }
    }
    return List::create(
               _["count"] = count,
               _["row"] = wrap(rows),
               _["col"] = wrap(cols)
           );
}

// For dgTMatrix / triplet sparse matrices
// [[Rcpp::export]]
List check_na_sparse_triplet_cpp(SEXP x, IntegerVector i, IntegerVector j,
                                 IntegerVector dim)
{
    if (Rf_isNull(x))
    {
        return List::create(
                   _["count"] = 0.0,
                   _["row"] = IntegerVector(0),
                   _["col"] = IntegerVector(0)
               );
    }
    R_xlen_t n = XLENGTH(x);
    std::vector<int> rows;
    std::vector<int> cols;
    double count = 0;
    for (R_xlen_t k = 0; k < n; ++k)
    {
        if (is_na_at(x, k))
        {
            ++count;
            rows.push_back(i[k] + 1);
            cols.push_back(j[k] + 1);
        }
    }
    return List::create(
               _["count"] = count,
               _["row"] = wrap(rows),
               _["col"] = wrap(cols)
           );
}