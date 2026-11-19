// src/BenchMemory.cpp
#include <Rcpp.h>
#include <cstdlib>
#include <fstream>
#include <string>

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
#include <sys/sysctl.h>
#include <sys/types.h>
#else
#include <unistd.h>
#endif

namespace {
constexpr double kGiB = 1024.0 * 1024.0 * 1024.0;

#if defined(__linux__)
// Parse "<key>   <number> kB" lines; returns true and sets `kb` on success.
bool parse_kb_line(const std::string &line, const std::string &key,
                   unsigned long long &kb) {
  if (line.compare(0, key.size(), key) != 0) {
    return false;
  }
  const char *start = line.c_str() + key.size();
  char *end = nullptr;
  kb = std::strtoull(start, &end, 10);
  return end != start;
}
#endif
} // namespace

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <psapi.h>
#include <windows.h>
#else
#include <sys/resource.h>
#include <sys/time.h>
#endif

//' Peak resident memory of the current R process
//'
//' Returns the peak resident set size (RSS) reached by the current R
//' process since it started, in megabytes (MB, i.e. 1024^2 bytes).
//'
//' @details
//' The value is obtained directly from the operating system:
//' \itemize{
//'   \item **Windows**: `PeakWorkingSetSize` from `GetProcessMemoryInfo()`.
//'   \item **macOS**: `ru_maxrss` from `getrusage(RUSAGE_SELF)`, reported
//'     in bytes by the kernel.
//'   \item **Linux and other Unix-alikes**: `ru_maxrss` from
//'     `getrusage(RUSAGE_SELF)`, reported in kilobytes by the kernel
//'     (equivalent to `VmHWM` in `/proc/self/status`).
//' }
//' The peak is monotonically non-decreasing over the lifetime of the
//' process; it is not reset by garbage collection.
//'
//' @return A single `numeric` value: the peak RSS in MB, or `NA_real_`
//'   if it cannot be determined on the current platform.
//'
//' @examples
//' \dontrun{
//' bench_peak_rss_mb()
//'
//' # Allocate roughly 400 MB and observe the peak grow
//' x <- numeric(5e7)
//' bench_peak_rss_mb()
//' }
//' @keywords internal
// [[Rcpp::export]]
double bench_peak_rss_mb() {
#if defined(_WIN32)
  PROCESS_MEMORY_COUNTERS pmc;
  pmc.cb = sizeof(pmc);
  if (!GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
    return NA_REAL;
  }
  return static_cast<double>(pmc.PeakWorkingSetSize) / (1024.0 * 1024.0);
#else
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) != 0) {
    return NA_REAL;
  }
#if defined(__APPLE__)
  // macOS: bytes
  return static_cast<double>(ru.ru_maxrss) / (1024.0 * 1024.0);
#else
  // Linux / BSD: kilobytes
  return static_cast<double>(ru.ru_maxrss) / 1024.0;
#endif
#endif
}

//' System memory overview
//'
//' Returns the total and the available physical memory of the machine,
//' in gigabytes (GB, i.e. 1024^3 bytes).
//'
//' @details
//' The values are obtained directly from the operating system:
//' \itemize{
//'   \item **Linux**: `MemTotal` and `MemAvailable` from `/proc/meminfo`
//'     (falls back to `sysconf()` for the total if `/proc` is unreadable;
//'     `MemAvailable` requires Linux >= 3.14).
//'   \item **macOS**: total from `sysctl(hw.memsize)`; available is
//'     approximated as free plus inactive pages reported by
//'     `host_statistics64()`.
//'   \item **Windows**: `ullTotalPhys` and `ullAvailPhys` from
//'     `GlobalMemoryStatusEx()`.
//'   \item **Other Unix-alikes**: total from `sysconf()`; available only if
//'     `_SC_AVPHYS_PAGES` is supported.
//' }
//' Note that "available" is defined slightly differently by each
//' operating system, so cross-platform comparisons should be made with care.
//'
//' @return A named `numeric` vector of length 2 with elements `MemTotal`
//'   and `MemAvailable` (GB). Elements that cannot be determined are
//'   `NA_real_`.
//'
//' @examples
//' \dontrun{
//' bench_proc_meminfo()
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector bench_proc_meminfo() {
  double total = NA_REAL;
  double avail = NA_REAL;

#if defined(_WIN32)
  MEMORYSTATUSEX ms;
  ms.dwLength = sizeof(ms);
  if (GlobalMemoryStatusEx(&ms)) {
    total = static_cast<double>(ms.ullTotalPhys) / kGiB;
    avail = static_cast<double>(ms.ullAvailPhys) / kGiB;
  }
#elif defined(__APPLE__)
  uint64_t memsize = 0;
  size_t len = sizeof(memsize);
  if (sysctlbyname("hw.memsize", &memsize, &len, nullptr, 0) == 0) {
    total = static_cast<double>(memsize) / kGiB;
  }

  mach_port_t host = mach_host_self();
  vm_size_t page_size = 0;
  vm_statistics64_data_t vm;
  mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
  if (host_page_size(host, &page_size) == KERN_SUCCESS &&
      host_statistics64(host, HOST_VM_INFO64,
                        reinterpret_cast<host_info64_t>(&vm),
                        &count) == KERN_SUCCESS) {
    const double pages = static_cast<double>(vm.free_count) +
                         static_cast<double>(vm.inactive_count);
    avail = pages * static_cast<double>(page_size) / kGiB;
  }
#elif defined(__linux__)
  std::ifstream f("/proc/meminfo");
  if (f.is_open()) {
    std::string line;
    bool got_total = false, got_avail = false;
    while (std::getline(f, line) && !(got_total && got_avail)) {
      unsigned long long kb = 0;
      if (!got_total && parse_kb_line(line, "MemTotal:", kb)) {
        total = static_cast<double>(kb) / (1024.0 * 1024.0);
        got_total = true;
      } else if (!got_avail && parse_kb_line(line, "MemAvailable:", kb)) {
        avail = static_cast<double>(kb) / (1024.0 * 1024.0);
        got_avail = true;
      }
    }
  } else {
    const long pages = sysconf(_SC_PHYS_PAGES);
    const long psize = sysconf(_SC_PAGESIZE);
    if (pages > 0 && psize > 0) {
      total = static_cast<double>(pages) * static_cast<double>(psize) / kGiB;
    }
  }
#else
  const long psize = sysconf(_SC_PAGESIZE);
#if defined(_SC_PHYS_PAGES)
  const long pages = sysconf(_SC_PHYS_PAGES);
  if (pages > 0 && psize > 0) {
    total = static_cast<double>(pages) * static_cast<double>(psize) / kGiB;
  }
#endif
#if defined(_SC_AVPHYS_PAGES)
  const long apages = sysconf(_SC_AVPHYS_PAGES);
  if (apages > 0 && psize > 0) {
    avail = static_cast<double>(apages) * static_cast<double>(psize) / kGiB;
  }
#endif
#endif

  return Rcpp::NumericVector::create(Rcpp::Named("MemTotal") = total,
                                     Rcpp::Named("MemAvailable") = avail);
}

//' System load average
//'
//' Returns the 1-minute system load average.
//'
//' @details
//' On Linux, macOS and other Unix-alikes the value comes from
//' `getloadavg()` (on Linux this is the first field of `/proc/loadavg`).
//' Windows has no load-average concept, so `NA_real_` is returned there.
//' Inside containers the value usually reflects the host, not the container.
//'
//' @return A single `numeric` value: the 1-minute load average, or
//'   `NA_real_` if it is unavailable on the current platform.
//'
//' @examples
//' \dontrun{
//' bench_proc_load()
//' }
//' @keywords internal
// [[Rcpp::export]]
double bench_proc_load() {
#if defined(_WIN32)
  return NA_REAL;
#else
  double load[1];
  if (getloadavg(load, 1) != 1) {
    return NA_REAL;
  }
  return load[0];
#endif
}