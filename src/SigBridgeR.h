// src/SigBridgeR.h
//
// Shared C++ helpers for the SigBridgeR package.
// Include this header from any .cpp that needs to emit terminal
// (verbose) progress messages while keeping the Rcpp R API available.
#ifndef SIGBRIDGER_SRC_SIGBRIDGER_H
#define SIGBRIDGER_SRC_SIGBRIDGER_H

#include <Rcpp.h>

#include <string>

// ============================================================
// cli-style colored terminal output (ANSI).
//
// Dispatches to the cli R package function cli::cli_alert_<kind>
// (e.g. kind = "info"  -> cli::cli_alert_info, a cyan "i" bullet;
//        kind = "success" -> cli::cli_alert_success, a green "v"
//                            bullet).
//
// The ANSI colors follow cli's num_colors auto-detection, so they
// are emitted only when the terminal supports them. Messages go to
// the message connection (like base message()).
//
// Callers should guard each call behind their own verbose flag, and
// kind must be one of the cli_alert_* functions exported by cli.
// ============================================================
inline void cli_emit(const char *kind, const std::string &msg) {
  Rcpp::Environment ns = Rcpp::Environment::namespace_env("cli");
  Rcpp::Function alert = ns[std::string("cli_alert_") + kind];
  alert(msg);
}

#endif // SIGBRIDGER_SRC_SIGBRIDGER_H
