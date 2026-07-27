#include <RcppArmadillo.h>
#include <unordered_map>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

static CharacterVector get_rownames(SEXP counts)
{
    if (Rf_isS4(counts))
    {
        S4 obj(counts);
        if (!obj.hasSlot("Dimnames"))
        {
            return CharacterVector();
        }
        SEXP dnS = obj.slot("Dimnames");
        if (Rf_isNull(dnS))
        {
            return CharacterVector();
        }
        List dn(dnS);
        if (dn.size() < 1)
        {
            return CharacterVector();
        }
        SEXP rnS = dn[0];
        if (Rf_isNull(rnS))
        {
            return CharacterVector();
        }
        return as<CharacterVector>(rnS);
    }
    else
    {
        SEXP dnS = Rf_getAttrib(counts, R_DimNamesSymbol);
        if (Rf_isNull(dnS))
        {
            return CharacterVector();
        }
        List dn(dnS);
        if (dn.size() < 1)
        {
            return CharacterVector();
        }
        SEXP rnS = dn[0];
        if (Rf_isNull(rnS))
        {
            return CharacterVector();
        }
        return as<CharacterVector>(rnS);
    }
}

static arma::vec make_inv_gene_length(
    CharacterVector rownames,
    NumericVector gene_length,
    int nrow
)
{
    if (rownames.size() != nrow)
    {
        stop("counts must have rownames, and length(rownames(counts)) must equal nrow(counts).");
    }
    SEXP namesS = gene_length.attr("names");
    if (Rf_isNull(namesS))
    {
        stop("gene_length must be a named numeric vector.");
    }
    CharacterVector gl_names(namesS);
    if (gl_names.size() != gene_length.size())
    {
        stop("gene_length names are invalid.");
    }
    std::unordered_map<std::string, double> len_map;
    len_map.reserve(static_cast<size_t>(gene_length.size()));
    for (R_xlen_t i = 0; i < gene_length.size(); ++i)
    {
        SEXP nmS = STRING_ELT(gl_names, i);
        if (nmS == NA_STRING)
        {
            stop("gene_length contains NA names.");
        }
        std::string nm(CHAR(nmS));
        if (nm.empty())
        {
            stop("gene_length contains empty names.");
        }
        double len = gene_length[i];
        if (!std::isfinite(len) || len <= 0.0)
        {
            stop("gene_length must contain positive finite gene lengths in bp.");
        }
        auto inserted = len_map.emplace(nm, len);
        if (!inserted.second)
        {
            stop("gene_length contains duplicated gene name: %s", nm);
        }
    }
    arma::vec inv_len(nrow);
    for (int i = 0; i < nrow; ++i)
    {
        SEXP rnS = STRING_ELT(rownames, i);
        if (rnS == NA_STRING)
        {
            stop("counts rownames contain NA.");
        }
        std::string gene(CHAR(rnS));
        auto it = len_map.find(gene);
        if (it == len_map.end())
        {
            stop("gene_length misses gene: %s", gene);
        }
        inv_len[i] = 1.0 / it->second;
    }
    return inv_len;
}

static NumericVector compute_tpm_dense(
    double *x_ptr,
    int nrow,
    int ncol,
    const arma::vec& inv_len
)
{
    arma::mat X(
        x_ptr,
        static_cast<arma::uword>(nrow),
        static_cast<arma::uword>(ncol),
        false,
        true
    );
    if (!X.is_finite())
    {
        stop("counts contains NA, NaN, or Inf.");
    }
    if (X.min() < 0.0)
    {
        stop("counts must be non-negative.");
    }
    R_xlen_t n_elem =
        static_cast<R_xlen_t>(nrow) * static_cast<R_xlen_t>(ncol);
    NumericVector out(n_elem);
    arma::mat Y(
        out.begin(),
        static_cast<arma::uword>(nrow),
        static_cast<arma::uword>(ncol),
        false,
        true
    );
    arma::rowvec denom = inv_len.t() * X;
    for (int j = 0; j < ncol; ++j)
    {
        double d = denom[j];
        if (!std::isfinite(d) || d < 0.0)
        {
            stop("invalid TPM denominator.");
        }
        if (d == 0.0)
        {
            Y.col(j).zeros();
        }
        else
        {
            Y.col(j) = (X.col(j) % inv_len) * (1e6 / d);
        }
    }
    return out;
}

static NumericVector compute_tpm_dgC(
    IntegerVector i_slot,
    IntegerVector p_slot,
    NumericVector x_slot,
    int nrow,
    int ncol,
    const arma::vec& inv_len
)
{
    if (p_slot.size() != ncol + 1)
    {
        stop("invalid dgCMatrix: length(p) must equal ncol + 1.");
    }
    if (i_slot.size() != x_slot.size())
    {
        stop("invalid dgCMatrix: length(i) must equal length(x).");
    }
    NumericVector out = clone(x_slot);
    std::vector<double> denom(ncol, 0.0);
    for (int col = 0; col < ncol; ++col)
    {
        int start = p_slot[col];
        int end   = p_slot[col + 1];
        for (int k = start; k < end; ++k)
        {
            int row = i_slot[k];
            if (row < 0 || row >= nrow)
            {
                stop("invalid dgCMatrix row index.");
            }
            double val = x_slot[k];
            if (!std::isfinite(val))
            {
                stop("counts contains NA, NaN, or Inf.");
            }
            if (val < 0.0)
            {
                stop("counts must be non-negative.");
            }
            denom[col] += val * inv_len[row];
        }
    }
    for (int col = 0; col < ncol; ++col)
    {
        double d = denom[col];
        if (!std::isfinite(d) || d < 0.0)
        {
            stop("invalid TPM denominator.");
        }
        double factor = d == 0.0 ? 0.0 : 1e6 / d;
        int start = p_slot[col];
        int end   = p_slot[col + 1];
        for (int k = start; k < end; ++k)
        {
            int row = i_slot[k];
            out[k] = x_slot[k] * inv_len[row] * factor;
        }
    }
    return out;
}

// [[Rcpp::export]]
SEXP CountsToTPM_impl(SEXP counts, NumericVector gene_length)
{
    if (Rf_isS4(counts))
    {
        S4 obj(counts);
        if (!obj.hasSlot("Dim"))
        {
            stop("S4 Matrix object must have slot 'Dim'.");
        }
        IntegerVector dim = obj.slot("Dim");
        if (dim.size() != 2)
        {
            stop("counts must be a 2D matrix.");
        }
        int nrow = dim[0];
        int ncol = dim[1];
        if (nrow <= 0 || ncol <= 0)
        {
            stop("counts must have positive nrow and ncol.");
        }
        CharacterVector rn = get_rownames(counts);
        arma::vec inv_len = make_inv_gene_length(rn, gene_length, nrow);
        if (Rf_inherits(counts, "dgeMatrix"))
        {
            NumericVector x_slot = obj.slot("x");
            R_xlen_t expected =
                static_cast<R_xlen_t>(nrow) * static_cast<R_xlen_t>(ncol);
            if (x_slot.size() != expected)
            {
                stop("invalid dgeMatrix: length(x) != nrow * ncol.");
            }
            NumericVector out_x =
                compute_tpm_dense(x_slot.begin(), nrow, ncol, inv_len);
            S4 ans = clone(obj);
            ans.slot("x") = out_x;
            if (ans.hasSlot("factors"))
            {
                ans.slot("factors") = List::create();
            }
            return ans;
        }
        if (Rf_inherits(counts, "dgCMatrix"))
        {
            IntegerVector i_slot = obj.slot("i");
            IntegerVector p_slot = obj.slot("p");
            NumericVector x_slot = obj.slot("x");
            NumericVector out_x =
                compute_tpm_dgC(i_slot, p_slot, x_slot, nrow, ncol, inv_len);
            S4 ans = clone(obj);
            ans.slot("x") = out_x;
            if (ans.hasSlot("factors"))
            {
                ans.slot("factors") = List::create();
            }
            return ans;
        }
        stop("S4 counts must be dgeMatrix or dgCMatrix after R-side coercion.");
    }
    if (!Rf_isMatrix(counts))
    {
        stop("counts must be a base matrix or an S4 Matrix object.");
    }
    NumericMatrix X = as<NumericMatrix>(counts);
    int nrow = X.nrow();
    int ncol = X.ncol();
    if (nrow <= 0 || ncol <= 0)
    {
        stop("counts must have positive nrow and ncol.");
    }
    CharacterVector rn = get_rownames(counts);
    arma::vec inv_len = make_inv_gene_length(rn, gene_length, nrow);
    NumericVector out = compute_tpm_dense(X.begin(), nrow, ncol, inv_len);
    out.attr("dim") = IntegerVector::create(nrow, ncol);
    SEXP dnS = Rf_getAttrib(counts, R_DimNamesSymbol);
    if (!Rf_isNull(dnS))
    {
        out.attr("dimnames") = dnS;
    }
    return out;
}