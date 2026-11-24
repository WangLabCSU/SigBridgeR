#include <Rcpp.h>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>

using namespace Rcpp;

void sexp_write(std::ostringstream &oss, SEXP x);

std::string type_or_class(SEXP x) {
  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);
  if (!Rf_isNull(cls) && Rf_length(cls) > 0) {
    CharacterVector cv(cls);
    std::string out = "object of class ";
    for (int i = 0; i < cv.size(); ++i) {
      if (i > 0)
        out += "/";
      if (cv[i] == NA_STRING) {
        out += "NA";
      } else {
        out += as<std::string>(cv[i]);
      }
    }
    return out;
  }
  switch (TYPEOF(x)) {
  case NILSXP:
    return "NULL";
  case SYMSXP:
    return "symbol";
  case LISTSXP:
    return "pairlist";
  case CLOSXP:
    return "function";
  case ENVSXP:
    return "environment";
  case PROMSXP:
    return "promise";
  case LANGSXP:
    return "call";
  case SPECIALSXP:
    return "special function";
  case BUILTINSXP:
    return "builtin function";
  case CHARSXP:
    return "CHARSXP";
  case LGLSXP:
    return "logical";
  case INTSXP:
    return "integer";
  case REALSXP:
    return "numeric";
  case CPLXSXP:
    return "complex";
  case STRSXP:
    return "character";
  case DOTSXP:
    return "...";
  case ANYSXP:
    return "any";
  case VECSXP:
    return "list";
  case EXPRSXP:
    return "expression";
  case BCODESXP:
    return "bytecode";
  case EXTPTRSXP:
    return "external pointer";
  case WEAKREFSXP:
    return "weak reference";
  case RAWSXP:
    return "raw";
  case S4SXP:
    return "S4 object";
  default:
    return "unknown type";
  }
}

void atomic_write(std::ostringstream &oss, SEXP x) {
  CharacterVector cv = as<CharacterVector>(x);
  for (int i = 0; i < cv.size(); ++i) {
    if (i > 0)
      oss << ", ";
    if (cv[i] == NA_STRING) {
      oss << "NA";
    } else {
      oss << as<std::string>(cv[i]);
    }
  }
}

void list_write(std::ostringstream &oss, List x) {
  SEXP nms_sexp = Rf_getAttrib(x, R_NamesSymbol);
  CharacterVector nms =
      Rf_isNull(nms_sexp) ? CharacterVector(0) : CharacterVector(nms_sexp);
  int n = x.size();
  for (int i = 0; i < n; ++i) {
    if (i > 0)
      oss << ", ";
    if (nms.size() == n && nms[i] != NA_STRING &&
        std::string(nms[i]).size() > 0) {
      oss << as<std::string>(nms[i]);
    } else {
      oss << i + 1;
    }
    oss << ": ";
    SEXP elem = x[i];
    if (Rf_isNewList(elem)) {
      oss << "(";
      list_write(oss, as<List>(elem));
      oss << ")";
    } else {
      sexp_write(oss, elem);
    }
  }
}

void sexp_write(std::ostringstream &oss, SEXP x) {
  if (Rf_isNull(x)) {
    oss << "NULL";
    return;
  }
  if (Rf_isNewList(x)) {
    oss << "(";
    list_write(oss, as<List>(x));
    oss << ")";
    return;
  }
  if (Rf_isVectorAtomic(x)) {
    atomic_write(oss, x);
    return;
  }
  oss << "<" << type_or_class(x) << ">";
}

//' Convert a List to a Compact Character Representation
//'
//' Converts an R `list` into a single character string using the format
//' `"key: value"`. Named elements use their names as keys, while unnamed
//' elements use their 1-based positions. Nested lists are converted
//' recursively and enclosed in parentheses.
//' This method is registered as `list_to_chr`, so it is called
//' automatically when `as.character()` is applied to a list.
//'
//' @param x A list to be converted.
//'
//' @return A character vector of length one.
//'
//' @examples
//' x <- list(key1 = "val1", key2 = "val2")
//' as.character(x)
//' # "key1: val1, key2: val2"
//' x <- list(key1 = list(key2 = list(key3 = "val3")))
//' as.character(x)
//' # "key1: (key2: (key3: val3))"
//' x <- list(
//'   a = 1,
//'   b = list(
//'     c = 2,
//'     d = list(e = "hello")
//'   ),
//'   f = c(1, 2, 3)
//' )
//' as.character(x)
//' # "a: 1, b: (c: 2, d: (e: hello)), f: 1, 2, 3"
//' x <- list("foo", bar = NA, baz = NULL)
//' as.character(x)
//' # "1: foo, bar: NA, baz: NULL"
//'
//' @export
// [[Rcpp::export]]
std::string list_to_chr(List x) {
  std::ostringstream oss;
  list_write(oss, x);
  return oss.str();
}

void list_write_formated(std::ostringstream &oss, List x, int depth,
                         bool &first) {
  SEXP nms_sexp = Rf_getAttrib(x, R_NamesSymbol);
  CharacterVector nms =
      Rf_isNull(nms_sexp) ? CharacterVector(0) : CharacterVector(nms_sexp);
  int n = x.size();
  std::string indent(static_cast<size_t>(depth) * 2, ' ');

  for (int i = 0; i < n; ++i) {
    if (!first)
      oss << "\n";
    first = false;

    oss << indent;
    if (nms.size() == n && nms[i] != NA_STRING &&
        std::string(nms[i]).size() > 0) {
      oss << as<std::string>(nms[i]);
    } else {
      oss << i + 1;
    }
    oss << ":";

    SEXP elem = x[i];
    if (Rf_isNull(elem)) {
      oss << " NULL";
    } else if (TYPEOF(elem) == VECSXP) {
      List sub(elem);
      if (sub.size() == 0) {
        oss << " ()";
      } else {
        list_write_formated(oss, sub, depth + 1, first);
      }
    } else {
      oss << " ";
      sexp_write(oss, elem);
    }
  }
}

//' Convert a List to an Indented Multi-line Character Representation
//'
//' Converts an R `list` into a single character string in which every
//' element occupies its own line, using the format `"key: value"`. Named
//' elements use their names as keys, while unnamed elements use their
//' 1-based positions. Nested lists are expanded recursively, and each
//' nesting level is indented by two additional spaces.
//'
//' @param x A list to be converted.
//'
//' @return A character vector of length one, with lines separated by `"\n"`
//'   (no trailing newline).
//'
//' @examples
//' x <- list(
//'   a = 1,
//'   b = list(
//'     c = 2,
//'     d = list(e = "hello")
//'   ),
//'   f = c(1, 2, 3)
//' )
//' cat(list_to_chr_formated(x), "\n")
//' # a: 1
//' # b:
//' #   c: 2
//' #   d:
//' #     e: hello
//' # f: 1, 2, 3
//'
//' @export
// [[Rcpp::export]]
std::string list_to_chr_formated(List x) {
  std::ostringstream oss;
  bool first = true;
  list_write_formated(oss, x, 0, first);
  return oss.str();
}

namespace {

// ANSI 前景色（与 cli 的 blue / grey 对应）
constexpr const char *ANSI_BLUE = "\033[34m";
constexpr const char *ANSI_GREY = "\033[90m";
constexpr const char *ANSI_RESET = "\033[39m";

inline void put_styled(std::ostringstream &oss, const char *code,
                       const std::string &text, bool color) {
  if (color) {
    oss << code << text << ANSI_RESET;
  } else {
    oss << text;
  }
}

// 模仿 cli 对向量的折叠：a / a and b / a, b, and c
std::string collapse_cli(const std::vector<std::string> &v) {
  const size_t n = v.size();
  if (n == 0)
    return "";
  if (n == 1)
    return v[0];
  if (n == 2)
    return v[0] + " and " + v[1];
  std::string out;
  for (size_t i = 0; i + 1 < n; ++i) {
    out += v[i];
    out += ", ";
  }
  out += "and ";
  out += v[n - 1];
  return out;
}

// class(x)[1L]
std::string first_class(SEXP x) {
  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);
  if (TYPEOF(cls) == STRSXP && Rf_length(cls) > 0) {
    return CHAR(STRING_ELT(cls, 0));
  }
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (!Rf_isNull(dim)) {
    return Rf_length(dim) == 2 ? "matrix" : "array";
  }
  switch (TYPEOF(x)) {
  case CLOSXP:
  case SPECIALSXP:
  case BUILTINSXP:
    return "function";
  case SYMSXP:
    return "name";
  default:
    return type_or_class(x); // 此时没有 class 属性，返回基础类型名
  }
}

// class(x) 的全部元素
std::vector<std::string> class_vector(SEXP x) {
  std::vector<std::string> out;
  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);
  if (TYPEOF(cls) == STRSXP && Rf_length(cls) > 0) {
    for (R_xlen_t i = 0; i < Rf_xlength(cls); ++i) {
      out.emplace_back(CHAR(STRING_ELT(cls, i)));
    }
  } else {
    out.push_back(first_class(x));
  }
  return out;
}

// rlang::fn_fmls_names()：闭包直接取 formals，原语先经 args()
std::vector<std::string> fn_formals(SEXP f) {
  Environment base = Environment::base_env();
  Function formals_fn = base["formals"];
  RObject fmls;
  if (TYPEOF(f) == CLOSXP) {
    fmls = formals_fn(f);
  } else {
    Function args_fn = base["args"];
    RObject a = args_fn(f);
    fmls = formals_fn(a);
  }
  std::vector<std::string> out;
  SEXP nms = Rf_getAttrib(fmls, R_NamesSymbol);
  if (TYPEOF(nms) == STRSXP) {
    for (R_xlen_t i = 0; i < Rf_xlength(nms); ++i) {
      out.emplace_back(CHAR(STRING_ELT(nms, i)));
    }
  }
  return out;
}

// 原子向量转字符串（最多 max_n 个）；带 class 的（factor / Date 等）走
// as.character
std::vector<std::string> atomic_items(SEXP x, R_xlen_t max_n) {
  CharacterVector cv;
  if (Rf_isObject(x)) {
    Function as_chr = Environment::base_env()["as.character"];
    cv = as_chr(x);
  } else {
    cv = as<CharacterVector>(x);
  }
  const R_xlen_t n = std::min<R_xlen_t>(cv.size(), max_n);
  std::vector<std::string> out;
  out.reserve(static_cast<size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) {
    out.emplace_back(cv[i] == NA_STRING ? "NA" : as<std::string>(cv[i]));
  }
  return out;
}

// S4 / S7 / R6 这类面向对象实例
inline bool is_oop_object(SEXP x) {
  const int t = TYPEOF(x);
  if (t == S4SXP || t == 25 /* OBJSXP, R >= 4.4 */)
    return true;
  return t == ENVSXP && Rf_isObject(x); // R6 等
}

inline bool is_2d_sexp(SEXP x) {
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  return !Rf_isNull(dim) && Rf_length(dim) == 2;
}

} // namespace

//' Print Method Backend for SigBridgeRBase
//'
//' Renders the class header and every S7 property of a `SigBridgeRBase`
//' object, one bullet per property, with ANSI colors when the terminal
//' supports them.
//'
//' @param cls_x Class name shown in the header.
//' @param nms_properties Property names.
//' @param properties A list of property values.
//'
//' @return `NULL`, invisibly. Called for its side effect of printing.
//' @noRd
// [[Rcpp::export]]
void print_sigbridger_base_cpp(std::string cls_x,
                               CharacterVector nms_properties,
                               List properties) {
  Environment cli = Environment::namespace_env("cli");

  Function num_ansi_colors = cli["num_ansi_colors"];
  Function is_utf8_output = cli["is_utf8_output"];

  const bool use_color = as<int>(num_ansi_colors()) > 1;
  const bool utf8 = as<bool>(is_utf8_output());
  const char *bullet = utf8 ? "\xE2\x80\xA2" : "*";

  std::ostringstream oss;
  oss << "<" << cls_x << ">\n";

  const R_xlen_t n = properties.size();
  const bool has_names = nms_properties.size() == n;

  for (R_xlen_t i = 0; i < n; ++i) {
    std::string name = (has_names && nms_properties[i] != NA_STRING)
                           ? as<std::string>(nms_properties[i])
                           : std::to_string(i + 1);
    SEXP key = properties[i];

    std::string inline_text;
    std::string block_text;

    if (TYPEOF(key) == CLOSXP || TYPEOF(key) == SPECIALSXP ||
        TYPEOF(key) == BUILTINSXP) {
      inline_text = collapse_cli(fn_formals(key));
    } else if (TYPEOF(key) == VECSXP) {
      std::ostringstream lo;
      bool first = true;
      list_write_formated(lo, List(key), 0, first);
      block_text = lo.str();
    } else if (is_oop_object(key)) {
      inline_text = collapse_cli(class_vector(key));
    } else if (is_2d_sexp(key)) {
      SEXP dim = Rf_getAttrib(key, R_DimSymbol);
      inline_text = "(" + std::to_string(INTEGER(dim)[0]) + " x " +
                    std::to_string(INTEGER(dim)[1]) + ")";
    } else if (Rf_isVectorAtomic(key)) {
      if (Rf_xlength(key) > 10) {
        std::vector<std::string> head = atomic_items(key, 5);
        for (size_t j = 0; j < head.size(); ++j) {
          inline_text += head[j];
          inline_text += ", ";
        }
        inline_text += "...";
      } else {
        inline_text = collapse_cli(atomic_items(key, 10));
      }
    } else if (!Rf_isNull(key)) {
      inline_text = "<" + type_or_class(key) + ">";
    }

    // • name [class]: value
    oss << bullet << " ";
    put_styled(oss, ANSI_BLUE, name, use_color);
    oss << " ";
    put_styled(oss, ANSI_GREY, "[" + first_class(key) + "]", use_color);
    oss << ":";

    if (!block_text.empty()) {
      std::istringstream iss(block_text);
      std::string line;
      while (std::getline(iss, line)) {
        oss << "\n    " << line;
      }
    } else if (!inline_text.empty()) {
      oss << " " << inline_text;
    }
    oss << "\n";
  }

  Rcout << oss.str();
}