#include <Rcpp.h>
#include <cstring>

// [[Rcpp::plugins(cpp17)]]

using namespace Rcpp;

// Check whether x inherits the "data.frame" class.
static bool is_data_frame(SEXP x) {
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
static bool is_plain_list(SEXP x) {
  return TYPEOF(x) == VECSXP && !is_data_frame(x);
}

// Check whether a name is neither NA nor empty.
static bool is_valid_name(SEXP name) {
  return name != NA_STRING && CHAR(name)[0] != '\0';
}

// Find the first matching name.
static R_xlen_t find_name(SEXP names, SEXP key) {
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
static bool get_valid_name_at(SEXP names, R_xlen_t i, SEXP *key) {
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
static List shallow_copy_list(SEXP x) {
  const R_xlen_t n = is_plain_list(x) ? XLENGTH(x) : 0;

  List out(n);

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
static List append_value(SEXP x, SEXP key, SEXP value, bool has_key) {
  const R_xlen_t n = XLENGTH(x);

  List out(n + 1);
  SEXP old_names = Rf_getAttrib(x, R_NamesSymbol);

  const bool make_names = has_key || old_names != R_NilValue;

  CharacterVector out_names;

  if (make_names) {
    out_names = CharacterVector(n + 1);

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
static List remove_at(SEXP x, R_xlen_t idx) {
  const R_xlen_t n = XLENGTH(x);

  if (idx < 0 || idx >= n) {
    return shallow_copy_list(x);
  }

  List out(n - 1);
  SEXP old_names = Rf_getAttrib(x, R_NamesSymbol);

  const bool has_names = old_names != R_NilValue;
  CharacterVector out_names;

  if (has_names) {
    out_names = CharacterVector(n - 1);
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
static List modify_list_rec(SEXP x, SEXP val, bool keep_null) {
  if (!is_plain_list(x)) {
    stop("`x` must be a list.");
  }

  if (!is_plain_list(val)) {
    stop("`val` must be a list.");
  }

  List out = shallow_copy_list(x);

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

      List merged = modify_list_rec(VECTOR_ELT(out, pos), value, keep_null);

      out[pos] = merged;

    } else if (pos >= 0) {
      out[pos] = value;

    } else {
      out = append_value(out, key, value, true);
    }
  }

  return out;
}

// [[Rcpp::export]]
List modify_list(List x, List val, bool keep_null = true) {
  if (!is_plain_list(x)) {
    stop("`x` must be a list.");
  }

  if (!is_plain_list(val)) {
    stop("`val` must be a list.");
  }

  return modify_list_rec(x, val, keep_null);
}
