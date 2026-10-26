// src/SigBridgeR.h
//
// Shared C++ helpers for the SigBridgeR package.
// Include this header from any .cpp that needs to emit terminal
// (verbose) progress messages while keeping the Rcpp R API available.
#ifndef SIGBRIDGER_SRC_SIGBRIDGER_H
#define SIGBRIDGER_SRC_SIGBRIDGER_H

#include <Rcpp.h>

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

#endif // SIGBRIDGER_SRC_SIGBRIDGER_H
