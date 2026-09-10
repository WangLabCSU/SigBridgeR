#include "Rtatami.h"
#include <R_ext/Arith.h>
#include <Rcpp.h>

#include <climits>
#include <string>
#include <vector>

// [[Rcpp::depends(beachmat, assorthead)]]
// [[Rcpp::plugins(cpp17)]]

using namespace Rcpp;

/*
 * Test whether a length-1 R object is NA.
 *
 * Used for vectors and data.frame columns.
 * Not used by the beachmat matrix scan.
 */
static inline bool is_one_na(SEXP x) {
  if (Rf_isNull(x) || XLENGTH(x) != 1) {
    return false;
  }

  switch (TYPEOF(x)) {
  case REALSXP:
    return ISNAN(REAL(x)[0]);

  case INTSXP:
    return INTEGER(x)[0] == NA_INTEGER;

  case LGLSXP:
    return LOGICAL(x)[0] == NA_LOGICAL;

  case STRSXP:
    return STRING_ELT(x, 0) == NA_STRING;

  case CPLXSXP: {
    const Rcomplex z = COMPLEX(x)[0];
    return ISNAN(z.r) || ISNAN(z.i);
  }

  default:
    return false;
  }
}

/*
 * Scan a plain R atomic vector, list column, or similar object.
 *
 * This function is not used for matrix-like objects.
 */
template <typename Add> static void scan_atomic(SEXP x, Add add) {
  const R_xlen_t n = XLENGTH(x);

  switch (TYPEOF(x)) {
  case REALSXP: {
    const double *p = REAL(x);

    for (R_xlen_t k = 0; k < n; ++k) {
      if (ISNAN(p[k])) {
        add(k);
      }
    }

    break;
  }

  case INTSXP: {
    const int *p = INTEGER(x);

    for (R_xlen_t k = 0; k < n; ++k) {
      if (p[k] == NA_INTEGER) {
        add(k);
      }
    }

    break;
  }

  case LGLSXP: {
    const int *p = LOGICAL(x);

    for (R_xlen_t k = 0; k < n; ++k) {
      if (p[k] == NA_LOGICAL) {
        add(k);
      }
    }

    break;
  }

  case STRSXP:
    for (R_xlen_t k = 0; k < n; ++k) {
      if (STRING_ELT(x, k) == NA_STRING) {
        add(k);
      }
    }
    break;

  case CPLXSXP: {
    const Rcomplex *p = COMPLEX(x);

    for (R_xlen_t k = 0; k < n; ++k) {
      if (ISNAN(p[k].r) || ISNAN(p[k].i)) {
        add(k);
      }
    }

    break;
  }

  case RAWSXP:
    // raw has no NA concept
    break;

  case VECSXP:
    for (R_xlen_t k = 0; k < n; ++k) {
      if (is_one_na(VECTOR_ELT(x, k))) {
        add(k);
      }
    }
    break;

  default:
    stop("Unsupported type for NA scan: ", Rf_type2char(TYPEOF(x)));
  }
}

/*
 * Append a position to a one-dimensional position result.
 *
 * R's IntegerVector uses int indices, so INT_MAX must be checked.
 */
static inline void append_vector_position(std::vector<int> &positions,
                                          R_xlen_t index) {
  if (index >= INT_MAX) {
    stop("Long vector index exceeds integer range.");
  }

  if (positions.size() >= static_cast<size_t>(INT_MAX)) {
    stop("Too many NA positions for an integer result.");
  }

  positions.push_back(static_cast<int>(index + 1));
}

/*
 * Append a row/column to a two-dimensional position result.
 */
static inline void append_matrix_position(std::vector<int> &rows,
                                          std::vector<int> &cols,
                                          std::size_t row, std::size_t col) {
  /*
   * row and col are 0-based.
   * The largest valid 1-based R index is INT_MAX.
   */
  if (row >= static_cast<std::size_t>(INT_MAX) ||
      col >= static_cast<std::size_t>(INT_MAX)) {
    stop("Matrix dimension exceeds integer index range.");
  }

  if (rows.size() >= static_cast<size_t>(INT_MAX)) {
    stop("Too many NA positions for an integer result.");
  }

  rows.push_back(static_cast<int>(row + 1));
  cols.push_back(static_cast<int>(col + 1));
}

/*
 * Scan a beachmat/tatami matrix.
 *
 * initialized_matrix must be the external pointer produced on the R side by
 *
 *     beachmat::initializeCpp(x)
 *
 * Data is fetched column by column via dense_column():
 *
 *     auto values = accessor->fetch(column, buffer.data());
 *
 * For sparse matrices, unstored positions are returned as the numeric 0;
 * stored NAs are returned as NA/NaN.
 */
// [[Rcpp::export]]
List check_na_beachmat_matrix_cpp(SEXP initialized_matrix) {
  Rtatami::BoundNumericPointer parsed(initialized_matrix);
  const auto &matrix = parsed->ptr;

  const std::size_t nr = static_cast<std::size_t>(matrix->nrow());

  const std::size_t nc = static_cast<std::size_t>(matrix->ncol());

  std::vector<int> rows;
  std::vector<int> cols;

  /*
   * Access column by column.
   *
   * For a dense matrix, the buffer usually does not cause an extra data
   * copy; for sparse or DelayedMatrix backends, the buffer receives the
   * materialized column.
   */
  auto accessor = matrix->dense_column();
  std::vector<double> buffer(nr);

  for (std::size_t col = 0; col < nc; ++col) {
    auto values = accessor->fetch(col, buffer.data());

    for (std::size_t row = 0; row < nr; ++row) {
      /*
       * ISNAN detects both NA_real_ and NaN,
       * consistent with the REALSXP behavior in the original code.
       */
      if (ISNAN(values[row])) {
        append_matrix_position(rows, cols, row, col);
      }
    }
  }

  const double count = static_cast<double>(rows.size());

  return List::create(_["count"] = count, _["row"] = wrap(rows),
                      _["col"] = wrap(cols));
}

/*
 * Scan a plain vector.
 *
 * This keeps the original R atomic-vector implementation, because
 * beachmat is a matrix backend rather than a general vector backend.
 */
// [[Rcpp::export]]
List check_na_vector_cpp(SEXP x) {
  std::vector<int> positions;

  scan_atomic(
      x, [&](R_xlen_t index) { append_vector_position(positions, index); });

  const double count = static_cast<double>(positions.size());

  return List::create(_["count"] = count, _["positions"] = wrap(positions));
}

/*
 * Scan a data.frame.
 *
 * A data.frame is not a beachmat matrix, so columns are scanned
 * one by one with scan_atomic().
 */
// [[Rcpp::export]]
List check_na_dataframe_cpp(List df, int nr) {
  if (nr < 0) {
    stop("nr must be non-negative.");
  }

  const int nc = df.size();

  std::vector<int> rows;
  std::vector<int> cols;

  for (int j = 0; j < nc; ++j) {
    SEXP column = df[j];

    if (XLENGTH(column) != nr) {
      stop("All data.frame columns must have length nrow(data).");
    }

    scan_atomic(column, [&](R_xlen_t row) {
      if (row >= INT_MAX) {
        stop("Row index exceeds integer range.");
      }

      if (rows.size() >= static_cast<size_t>(INT_MAX)) {
        stop("Too many NA positions for an integer result.");
      }

      rows.push_back(static_cast<int>(row + 1));
      cols.push_back(j + 1);
    });
  }

  const double count = static_cast<double>(rows.size());

  return List::create(_["count"] = count, _["row"] = wrap(rows),
                      _["col"] = wrap(cols));
}