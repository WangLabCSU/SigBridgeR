#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

static inline double safe_div(double num, double den) {
  return den > 0.0 ? num / den : NA_REAL;
}
//' Calculate Binary Classification Metrics for Screening Labels
//'
//' Treats `"Positive"` as the positive class and computes accuracy metrics
//' for each column of predicted labels.
//'
//' @param metrics A data frame of predicted-label columns, one per screening
//'   method.
//' @param bench_label Character vector of ground-truth labels, one per row of
//'   `metrics`.
//' @param n_threads `integer`. Number of OpenMP threads. Default: `1L`.
//'
//' @return A named list with one element per column of `metrics`. Each
//'   element is a 1-row numeric matrix with columns `TPR`, `FPR`, `TNR`,
//'   `FNR`, `F1` and `Accuracy` and the method name as row name.
//'
//' @keywords internal
// [[Rcpp::export]]
List calc_binary_metrics(DataFrame metrics, CharacterVector bench_label,
                         int n_threads = 1) {
  const R_xlen_t n = bench_label.size();
  const int p = metrics.size();

  if (p == 0)
    stop("`metrics` has no columns.");
  if (metrics.nrow() != n)
    stop("nrow(metrics) (%d) must equal length(bench_label) (%d).",
         (int)metrics.nrow(), (int)n);

  CharacterVector keys =
      CharacterVector::create("Positive", "Negative", "Neutral", "Other");
  const SEXP K0 = STRING_ELT(keys, 0), K1 = STRING_ELT(keys, 1),
             K2 = STRING_ELT(keys, 2), K3 = STRING_ELT(keys, 3);

  const SEXP *bp = STRING_PTR_RO(bench_label);
  std::vector<unsigned char> truth(n);
  double P = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    const SEXP s = bp[i];
    if (s == K0) {
      truth[i] = 1;
      P += 1;
    } else if (s == K1 || s == K2 || s == K3)
      truth[i] = 0;
    else
      stop("`bench_label` contains invalid value or NA at index %d.",
           (int)(i + 1));
  }

  CharacterVector col_names = metrics.names();
  std::vector<const SEXP *> cols(p);
  for (int j = 0; j < p; ++j) {
    SEXP col = metrics[j];
    if (TYPEOF(col) != STRSXP)
      stop("Column '%s' of `metrics` must be a character vector (not factor).",
           CHAR(STRING_ELT(col_names, j)));
    cols[j] = STRING_PTR_RO(col);
  }

  std::vector<double> TPv(p), PPv(p);
  std::vector<R_xlen_t> bad(p, -1);
  const unsigned char *t = truth.data();

#pragma omp parallel for num_threads(n_threads)                                \
    schedule(static) if (n_threads > 1 && p > 1)
  for (int j = 0; j < p; ++j) {
    const SEXP *x = cols[j];
    R_xlen_t TP = 0, PP = 0;
    for (R_xlen_t i = 0; i < n; ++i) {
      const SEXP s = x[i];
      const bool isP = (s == K0);
      if (!(isP | (s == K1) | (s == K2) | (s == K3))) {
        bad[j] = i;
        break;
      }
      PP += isP;
      TP += (isP & (t[i] != 0));
    }
    TPv[j] = (double)TP;
    PPv[j] = (double)PP;
  }

  for (int j = 0; j < p; ++j)
    if (bad[j] >= 0)
      stop("Column '%s' contains invalid value or NA at row %d.",
           CHAR(STRING_ELT(col_names, j)), (int)(bad[j] + 1));

  CharacterVector stat_names =
      CharacterVector::create("TPR", "FPR", "TNR", "FNR", "F1", "Accuracy");
  const double N = (double)n;

  List out(p);
  for (int j = 0; j < p; ++j) {
    const double TP = TPv[j];
    const double FN = P - TP;
    const double FP = PPv[j] - TP;
    const double TN = N - P - FP;

    NumericMatrix m(1, 6);
    m(0, 0) = safe_div(TP, TP + FN);
    m(0, 1) = safe_div(FP, FP + TN);
    m(0, 2) = safe_div(TN, TN + FP);
    m(0, 3) = safe_div(FN, FN + TP);
    m(0, 4) = safe_div(2.0 * TP, 2.0 * TP + FP + FN);
    m(0, 5) = safe_div(TP + TN, N);
    colnames(m) = stat_names;
    rownames(m) = CharacterVector::create(col_names[j]);
    out[j] = m;
  }
  out.attr("names") = col_names;
  return out;
}