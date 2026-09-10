#include "Rtatami.h"
#include <Rcpp.h>

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

// [[Rcpp::depends(beachmat, assorthead)]]
// [[Rcpp::plugins(cpp17)]]

using namespace Rcpp;

enum AggMethod { AGG_MAX, AGG_SUM, AGG_MEAN, AGG_MEDIAN, AGG_FIRST };

static AggMethod parse_method(const std::string &method) {
  if (method == "max")
    return AGG_MAX;
  if (method == "sum")
    return AGG_SUM;
  if (method == "mean")
    return AGG_MEAN;
  if (method == "median")
    return AGG_MEDIAN;
  if (method == "first")
    return AGG_FIRST;
  stop("Unsupported method: " + method);
  return AGG_MAX;
}

static std::string name_key(SEXP s) {
  if (s == NA_STRING) {
    return std::string("\001NA");
  }
  return std::string("\002") + Rf_translateCharUTF8(s);
}

struct Groups {
  std::vector<std::vector<int>> index;
  std::vector<int> first_pos;
};

static Groups make_groups(CharacterVector names) {
  std::unordered_map<std::string, int> seen;
  Groups groups;

  const int n = names.size();

  for (int i = 0; i < n; ++i) {
    std::string key = name_key(STRING_ELT(names, i));

    auto it = seen.find(key);

    if (it == seen.end()) {
      int gid = groups.index.size();

      seen[key] = gid;
      groups.index.emplace_back();
      groups.first_pos.push_back(i);
      groups.index[gid].push_back(i);
    } else {
      groups.index[it->second].push_back(i);
    }
  }

  return groups;
}

static double median_no_na(std::vector<double> &values) {
  if (values.empty()) {
    return NA_REAL;
  }

  std::sort(values.begin(), values.end());

  const int n = values.size();

  if (n % 2 == 1) {
    return values[n / 2];
  }

  return (values[n / 2 - 1] + values[n / 2]) / 2.0;
}

static double aggregate_values(const std::vector<double> &values,
                               const std::vector<int> &indices,
                               AggMethod method) {
  if (indices.empty()) {
    return NA_REAL;
  }

  if (method == AGG_FIRST) {
    return values[indices[0]];
  }

  if (method == AGG_SUM) {
    double sum = 0.0;

    for (int i : indices) {
      double value = values[i];

      if (!ISNAN(value)) {
        sum += value;
      }
    }

    return sum;
  }

  if (method == AGG_MEAN) {
    double sum = 0.0;
    int count = 0;

    for (int i : indices) {
      double value = values[i];

      if (!ISNAN(value)) {
        sum += value;
        ++count;
      }
    }

    return count == 0 ? R_NaN : sum / count;
  }

  if (method == AGG_MAX) {
    double best = R_NegInf;
    bool found = false;

    for (int i : indices) {
      double value = values[i];

      if (!ISNAN(value)) {
        if (!found || value > best) {
          best = value;
        }
        found = true;
      }
    }

    return found ? best : NA_REAL;
  }

  // median
  std::vector<double> selected;
  selected.reserve(indices.size());

  for (int i : indices) {
    double value = values[i];

    if (!ISNAN(value)) {
      selected.push_back(value);
    }
  }

  return median_no_na(selected);
}

// [[Rcpp::export]]
NumericMatrix aggregate_dup_cols_cpp(SEXP initialized_matrix,
                                     CharacterVector col_names,
                                     std::string method) {
  Rtatami::BoundNumericPointer parsed(initialized_matrix);
  const auto &matrix = parsed->ptr;

  const int nr = matrix->nrow();
  const int nc = matrix->ncol();

  if (col_names.size() != nc) {
    stop("Length of col_names must equal ncol(x).");
  }

  AggMethod agg = parse_method(method);
  Groups groups = make_groups(col_names);

  const int ng = groups.index.size();
  NumericMatrix output(nr, ng);

  auto accessor = matrix->dense_column();

  std::vector<std::vector<double>> columns(nc, std::vector<double>(nr));

  std::vector<double> buffer(nr);

  for (int j = 0; j < nc; ++j) {
    auto values = accessor->fetch(j, buffer.data());

    for (int i = 0; i < nr; ++i) {
      columns[j][i] = values[i];
    }
  }

  std::vector<double> values(nc);

  for (int g = 0; g < ng; ++g) {
    const std::vector<int> &indices = groups.index[g];

    for (int i = 0; i < nr; ++i) {
      for (int j = 0; j < nc; ++j) {
        values[j] = columns[j][i];
      }

      output(i, g) = aggregate_values(values, indices, agg);
    }
  }

  return output;
}

// [[Rcpp::export]]
NumericMatrix aggregate_dup_rows_cpp(SEXP initialized_matrix,
                                     CharacterVector row_names,
                                     std::string method) {
  Rtatami::BoundNumericPointer parsed(initialized_matrix);
  const auto &matrix = parsed->ptr;

  const int nr = matrix->nrow();
  const int nc = matrix->ncol();

  if (row_names.size() != nr) {
    stop("Length of row_names must equal nrow(x).");
  }

  AggMethod agg = parse_method(method);
  Groups groups = make_groups(row_names);

  const int ng = groups.index.size();
  NumericMatrix output(ng, nc);

  auto accessor = matrix->dense_column();
  std::vector<double> buffer(nr);

  for (int j = 0; j < nc; ++j) {
    auto values = accessor->fetch(j, buffer.data());

    for (int g = 0; g < ng; ++g) {
      const std::vector<int> &indices = groups.index[g];

      std::vector<double> selected;
      selected.reserve(indices.size());

      for (int i : indices) {
        selected.push_back(values[i]);
      }

      if (indices.size() == 1 || agg == AGG_FIRST) {
        output(g, j) = values[indices[0]];
      } else {
        std::vector<int> local_indices(selected.size());

        for (size_t k = 0; k < selected.size(); ++k) {
          local_indices[k] = k;
        }

        output(g, j) = aggregate_values(selected, local_indices, agg);
      }
    }
  }

  return output;
}