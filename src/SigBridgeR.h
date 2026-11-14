// src/SigBridgeR.h
//
// Shared C++ helpers for the SigBridgeR package.
// Include this header from any .cpp that needs to emit terminal
// (verbose) progress messages while keeping the Rcpp R API available,
// or that needs the recursive list-modification helpers.
#ifndef SIGBRIDGER_SRC_SIGBRIDGER_H
#define SIGBRIDGER_SRC_SIGBRIDGER_H

#include <Rcpp.h>

#include <cstring>
#include <string>

// NOTE: 调用 R API，只能在主线程使用（不要在并行 worker 内调用）。

// 转义 glue/cli 模板中的花括号，避免 msg 被当作表达式求值
inline std::string cli_escape(const std::string &s) {
  std::string out;
  out.reserve(s.size() + 8);
  for (const char c : s) {
    out.push_back(c);
    if (c == '{' || c == '}') {
      out.push_back(c); // { -> {{, } -> }}
    }
  }
  return out;
}

inline void cli_emit(const std::string &kind, const std::string &msg) {
  if (kind != "info" && kind != "success" && kind != "warning" &&
      kind != "danger") {
    Rcpp::stop("cli_emit: unsupported kind '%s'.", kind);
  }
  Rcpp::Environment ns = Rcpp::Environment::namespace_env("cli");
  Rcpp::Function alert = ns["cli_alert_" + kind];
  alert(cli_escape(msg));
}

// ---------------------------------------------------------------------------
// Recursive list modification (implementation moved from modify_list.cpp).
//
// 以 inline 形式放在头文件：任何 .cpp 只要 #include "SigBridgeR.h" 就能调用
// sigbridger::modify_list(x, val, keep_null)，不依赖特定编译单元。
// src/modify_list.cpp 仅保留 // [[Rcpp::export]] 薄封装，因为
// Rcpp::compileAttributes() 只扫描 .cpp 文件，头文件中的导出属性不会生成
// R 端 wrapper，也不能在头文件里重复定义默认参数。
// ---------------------------------------------------------------------------
namespace sigbridger {

// Check whether x inherits the "data.frame" class.
inline bool is_data_frame(SEXP x) {
  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);

  if (cls == R_NilValue || TYPEOF(cls) != STRSXP) {
    return false;
  }

  for (R_xlen_t i = 0; i < XLENGTH(cls); ++i) {
    SEXP name = STRING_ELT(cls, i);

    if (name != NA_STRING && std::strcmp(CHAR(name), "data.frame") == 0) {
      return true;
    }
  }

  return false;
}

// Check whether x is a plain list.
inline bool is_plain_list(SEXP x) {
  return TYPEOF(x) == VECSXP && !is_data_frame(x);
}

// Check whether a name is neither NA nor empty.
inline bool is_valid_name(SEXP name) {
  return name != NA_STRING && CHAR(name)[0] != '\0';
}

// Find the first matching name.
inline R_xlen_t find_name(SEXP names, SEXP key) {
  if (names == R_NilValue || TYPEOF(names) != STRSXP || key == NA_STRING) {
    return -1;
  }

  const char *k = CHAR(key);

  for (R_xlen_t i = 0; i < XLENGTH(names); ++i) {
    SEXP name = STRING_ELT(names, i);

    if (name == NA_STRING) {
      continue;
    }

    if (std::strcmp(CHAR(name), k) == 0) {
      return i;
    }
  }

  return -1;
}

// Get a valid name at position i.
inline bool get_valid_name_at(SEXP names, R_xlen_t i, SEXP *key) {
  if (names == R_NilValue || TYPEOF(names) != STRSXP || i >= XLENGTH(names)) {
    return false;
  }

  SEXP name = STRING_ELT(names, i);

  if (!is_valid_name(name)) {
    return false;
  }

  *key = name;
  return true;
}

// Create a shallow copy of a plain list.
// Only names and class are retained, matching the original behavior.
inline Rcpp::List shallow_copy_list(SEXP x) {
  const R_xlen_t n = is_plain_list(x) ? XLENGTH(x) : 0;

  Rcpp::List out(n);

  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = VECTOR_ELT(x, i);
  }

  if (is_plain_list(x)) {
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);

    if (names != R_NilValue) {
      out.attr("names") = names;
    }
  }

  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);

  if (cls != R_NilValue) {
    out.attr("class") = cls;
  }

  return out;
}

// Append one value to a list.
inline Rcpp::List append_value(SEXP x, SEXP key, SEXP value, bool has_key) {
  const R_xlen_t n = XLENGTH(x);

  Rcpp::List out(n + 1);
  SEXP old_names = Rf_getAttrib(x, R_NamesSymbol);

  const bool make_names = has_key || old_names != R_NilValue;

  Rcpp::CharacterVector out_names;

  if (make_names) {
    out_names = Rcpp::CharacterVector(n + 1);

    for (R_xlen_t i = 0; i < n; ++i) {
      if (old_names != R_NilValue && i < XLENGTH(old_names)) {
        out_names[i] = STRING_ELT(old_names, i);
      } else {
        out_names[i] = "";
      }
    }

    out_names[n] = has_key ? key : Rf_mkChar("");
  }

  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = VECTOR_ELT(x, i);
  }

  out[n] = value;

  if (make_names) {
    out.attr("names") = out_names;
  }

  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);

  if (cls != R_NilValue) {
    out.attr("class") = cls;
  }

  return out;
}

// Remove the element at position idx.
inline Rcpp::List remove_at(SEXP x, R_xlen_t idx) {
  const R_xlen_t n = XLENGTH(x);

  if (idx < 0 || idx >= n) {
    return shallow_copy_list(x);
  }

  Rcpp::List out(n - 1);
  SEXP old_names = Rf_getAttrib(x, R_NamesSymbol);

  const bool has_names = old_names != R_NilValue;
  Rcpp::CharacterVector out_names;

  if (has_names) {
    out_names = Rcpp::CharacterVector(n - 1);
  }

  R_xlen_t j = 0;

  for (R_xlen_t i = 0; i < n; ++i) {
    if (i == idx) {
      continue;
    }

    out[j] = VECTOR_ELT(x, i);

    if (has_names) {
      if (i < XLENGTH(old_names)) {
        out_names[j] = STRING_ELT(old_names, i);
      } else {
        out_names[j] = "";
      }
    }

    ++j;
  }

  if (has_names) {
    out.attr("names") = out_names;
  }

  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);

  if (cls != R_NilValue) {
    out.attr("class") = cls;
  }

  return out;
}

// Recursive implementation.
inline Rcpp::List modify_list_rec(SEXP x, SEXP val, bool keep_null) {
  if (!is_plain_list(x)) {
    Rcpp::stop("`x` must be a list.");
  }

  if (!is_plain_list(val)) {
    Rcpp::stop("`val` must be a list.");
  }

  Rcpp::List out = shallow_copy_list(x);

  SEXP val_names = Rf_getAttrib(val, R_NamesSymbol);

  for (R_xlen_t i = 0; i < XLENGTH(val); ++i) {
    SEXP key = R_NilValue;

    const bool has_key = get_valid_name_at(val_names, i, &key);

    SEXP value = VECTOR_ELT(val, i);

    // Unnamed elements are appended without matching.
    if (!has_key) {
      if (value == R_NilValue && !keep_null) {
        continue;
      }

      out = append_value(out, key, value, false);
      continue;
    }

    SEXP out_names = Rf_getAttrib(out, R_NamesSymbol);
    const R_xlen_t pos = find_name(out_names, key);

    // NULL deletes the first matching element.
    if (value == R_NilValue && !keep_null) {
      if (pos >= 0) {
        out = remove_at(out, pos);
      }

      continue;
    }

    // Recursively merge when both elements are plain lists.
    if (pos >= 0 && is_plain_list(VECTOR_ELT(out, pos)) &&
        is_plain_list(value)) {

      Rcpp::List merged =
          modify_list_rec(VECTOR_ELT(out, pos), value, keep_null);

      out[pos] = merged;

    } else if (pos >= 0) {
      out[pos] = value;

    } else {
      out = append_value(out, key, value, true);
    }
  }

  return out;
}

// Public entry point: recursively merge `val` into `x`.
// keep_null = true 保留 NULL；false 时命名元素的 NULL 表示删除，
// 未命名元素的 NULL 直接忽略。
inline Rcpp::List modify_list(Rcpp::List x, Rcpp::List val,
                              bool keep_null = true) {
  if (!is_plain_list(x)) {
    Rcpp::stop("`x` must be a list.");
  }

  if (!is_plain_list(val)) {
    Rcpp::stop("`val` must be a list.");
  }

  return modify_list_rec(x, val, keep_null);
}

} // namespace sigbridger

#endif // SIGBRIDGER_SRC_SIGBRIDGER_H
