#include <R_ext/Arith.h>
#include <Rcpp.h>
#include "Rtatami.h"

#include <climits>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

// [[Rcpp::depends(beachmat, assorthead)]]
// [[Rcpp::plugins(cpp17)]]

using namespace Rcpp;

/*
 * 根据 counts 的行名生成每个基因的 1 / gene_length。
 */
static std::vector<double> make_inv_gene_length(CharacterVector row_names,
                                                NumericVector gene_length,
                                                std::size_t nrow) {
  if (row_names.size() != static_cast<R_xlen_t>(nrow)) {
    stop("counts must have rownames, and length(rownames(counts)) "
         "must equal nrow(counts).");
  }

  SEXP names_sexp = gene_length.attr("names");

  if (Rf_isNull(names_sexp)) {
    stop("gene_length must be a named numeric vector.");
  }

  CharacterVector gene_length_names(names_sexp);

  if (gene_length_names.size() != gene_length.size()) {
    stop("gene_length names are invalid.");
  }

  std::unordered_map<std::string, double> length_map;

  length_map.reserve(static_cast<std::size_t>(gene_length.size()));

  for (R_xlen_t k = 0; k < gene_length.size(); ++k) {
    SEXP name_sexp = STRING_ELT(gene_length_names, k);

    if (name_sexp == NA_STRING) {
      stop("gene_length contains NA names.");
    }

    const char *name_chars = Rf_translateCharUTF8(name_sexp);

    std::string gene_name(name_chars);

    if (gene_name.empty()) {
      stop("gene_length contains empty names.");
    }

    const double length = gene_length[k];

    if (!std::isfinite(length) || length <= 0.0) {
      stop("gene_length must contain positive finite "
           "gene lengths in bp.");
    }

    const auto inserted = length_map.emplace(gene_name, length);

    if (!inserted.second) {
      stop("gene_length contains duplicated gene name: %s", gene_name);
    }
  }

  std::vector<double> inv_length(nrow);

  for (std::size_t i = 0; i < nrow; ++i) {
    SEXP row_name_sexp = STRING_ELT(row_names, static_cast<R_xlen_t>(i));

    if (row_name_sexp == NA_STRING) {
      stop("counts rownames contain NA.");
    }

    const char *row_name_chars = Rf_translateCharUTF8(row_name_sexp);

    std::string gene_name(row_name_chars);

    auto it = length_map.find(gene_name);

    if (it == length_map.end()) {
      stop("gene_length misses gene: %s", gene_name);
    }

    inv_length[i] = 1.0 / it->second;
  }

  return inv_length;
}

/*
 * 通过 beachmat/tatami 读取矩阵，并转换为 TPM。
 *
 * initialized_counts 必须是 R 端：
 *
 *     beachmat::initializeCpp(counts)
 *
 * 返回的 external pointer。
 *
 * row_names 由 R 端传入，因为 beachmat external pointer
 * 主要负责矩阵数据访问，并不负责 dimnames。
 */
// [[Rcpp::export]]
NumericMatrix counts_to_tpm_cpp(SEXP initialized_counts,
                                CharacterVector row_names,
                                NumericVector gene_length) {
  /*
   * 将 R 端的 external pointer 解析成 tatami 数值矩阵。
   */
  Rtatami::BoundNumericPointer parsed(initialized_counts);

  /*
   * 使用局部 shared pointer。
   */
  auto matrix = parsed->ptr;

  const std::size_t nrow = static_cast<std::size_t>(matrix->nrow());

  const std::size_t ncol = static_cast<std::size_t>(matrix->ncol());

  if (nrow == 0 || ncol == 0) {
    stop("counts must have positive nrow and ncol.");
  }

  /*
   * NumericMatrix 的维度使用 int。
   */
  if (nrow > static_cast<std::size_t>(INT_MAX) ||
      ncol > static_cast<std::size_t>(INT_MAX)) {
    stop("counts dimensions exceed R integer dimension limits.");
  }

  const std::vector<double> inv_length =
      make_inv_gene_length(row_names, gene_length, nrow);

  /*
   * 第一次遍历：计算每个样本的 TPM 分母。
   *
   * TPM 的计算公式：
   *
   *   RPK_i = count_i / gene_length_i
   *   TPM_i = RPK_i / sum(RPK) * 1e6
   *
   * denom[j] = sum_i(count[i,j] * inv_length[i])
   */
  std::vector<double> denominator(ncol, 0.0);

  std::vector<double> buffer(nrow);

  /*
   * 使用列访问器。
   *
   * 对于稀疏矩阵，未存储的位置会以 0 返回；
   * 对于 DelayedMatrix，fetch() 会请求对应列的数据。
   */
  auto accessor = matrix->dense_column();

  for (std::size_t col = 0; col < ncol; ++col) {
    auto values = accessor->fetch(col, buffer.data());

    double denom = 0.0;

    for (std::size_t row = 0; row < nrow; ++row) {
      const double value = values[row];

      if (!std::isfinite(value)) {
        stop("counts contains NA, NaN, or Inf.");
      }

      if (value < 0.0) {
        stop("counts must be non-negative.");
      }

      denom += value * inv_length[row];
    }

    denominator[col] = denom;
  }

  /*
   * 第二次遍历：生成输出矩阵。
   *
   * 输出为普通 dense NumericMatrix。
   */
  NumericMatrix output(static_cast<int>(nrow), static_cast<int>(ncol));

  for (std::size_t col = 0; col < ncol; ++col) {
    const double denom = denominator[col];

    if (!std::isfinite(denom) || denom < 0.0) {
      stop("invalid TPM denominator.");
    }

    auto values = accessor->fetch(col, buffer.data());

    /*
     * 如果某个样本的所有 counts 都为 0，
     * 则该样本的 TPM 全部设为 0。
     */
    if (denom == 0.0) {
      for (std::size_t row = 0; row < nrow; ++row) {
        output(static_cast<int>(row), static_cast<int>(col)) = 0.0;
      }

      continue;
    }

    const double factor = 1e6 / denom;

    for (std::size_t row = 0; row < nrow; ++row) {
      const double value = values[row];

      output(static_cast<int>(row), static_cast<int>(col)) =
          value * inv_length[row] * factor;
    }
  }

  return output;
}