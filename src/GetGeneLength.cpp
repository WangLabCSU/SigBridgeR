#include <Rcpp.h>
#include <zlib.h>

#include <algorithm>
#include <climits>
#include <cstring>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "SigBridgeR.h"

using namespace Rcpp;
using namespace std;

struct Interval {
  int start;
  int end;
};

struct Field {
  const char *begin;
  const char *end;

  size_t size() const { return static_cast<size_t>(end - begin); }
};

using GroupMap = unordered_map<string, vector<Interval>>;

// ============================================================
// GTF/GTF.GZ line reader
// ============================================================

class LineReader {
private:
  bool is_gz_;
  gzFile gz_fp_;
  ifstream file_;
  vector<char> buffer_;

public:
  explicit LineReader(const string &path)
      : is_gz_(false), gz_fp_(nullptr), buffer_(1 << 16) {

    const bool gz =
        path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;

    if (gz) {
      is_gz_ = true;
      gz_fp_ = gzopen(path.c_str(), "rb");

      if (gz_fp_ == nullptr) {
        stop("Cannot open gzip file: " + path);
      }

    } else {
      file_.open(path.c_str());

      if (!file_) {
        stop("Cannot open file: " + path);
      }
    }
  }

  ~LineReader() {
    if (gz_fp_ != nullptr) {
      gzclose(gz_fp_);
      gz_fp_ = nullptr;
    }
  }

  bool next(string &line) {

    line.clear();

    // plain text file
    if (!is_gz_) {
      if (!getline(file_, line)) {
        return false;
      }

      if (!line.empty() && line.back() == '\r') {
        line.pop_back();
      }

      return true;
    }

    // gzip file
    while (true) {

      char *p =
          gzgets(gz_fp_, buffer_.data(), static_cast<int>(buffer_.size()));

      if (p == nullptr) {
        if (!line.empty()) {
          break;
        }

        return false;
      }

      const size_t n = strlen(p);

      line.append(p, n);

      // newline indicates the current line has ended
      if (n > 0 && p[n - 1] == '\n') {
        break;
      }

      // the last line may not have a newline
      if (gzeof(gz_fp_)) {
        break;
      }
    }

    while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
      line.pop_back();
    }

    return true;
  }
};

// ============================================================
// Fast chromosome validation
// Keep 1-22, X, Y, MT
// ============================================================

inline bool valid_chr_fast(const Field &chr) {

  const size_t n = chr.size();

  if (n == 1 && (chr.begin[0] == 'X' || chr.begin[0] == 'Y')) {
    return true;
  }

  if (n == 2 && chr.begin[0] == 'M' && chr.begin[1] == 'T') {
    return true;
  }

  if (n == 0) {
    return false;
  }

  int value = 0;

  for (const char *p = chr.begin; p != chr.end; ++p) {

    if (*p < '0' || *p > '9') {
      return false;
    }

    const int digit = *p - '0';

    // prevent overflow and exclude numbers greater than 22
    if (value > 2 || (value == 2 && digit > 2)) {
      return false;
    }

    value = value * 10 + digit;
  }

  return value >= 1 && value <= 22;
}

// ============================================================
// Fast integer parsing
// ============================================================

inline bool parse_int_fast(const Field &field, int &value) {

  if (field.begin == field.end) {
    return false;
  }

  long long x = 0;

  for (const char *p = field.begin; p != field.end; ++p) {

    if (*p < '0' || *p > '9') {
      return false;
    }

    x = x * 10 + (*p - '0');

    if (x > INT_MAX) {
      return false;
    }
  }

  value = static_cast<int>(x);

  return true;
}

// ============================================================
// gene_id attribute parsing
// ============================================================

inline bool is_attr_separator(char c) {
  return c == ' ' || c == '\t' || c == ';';
}

inline bool get_gene_id_fast(const Field &attr, string &gene_id) {

  static const char target[] = "gene_id";
  static const size_t target_len = 7;

  if (attr.size() < target_len) {
    return false;
  }

  const char *hit = nullptr;

  for (const char *p = attr.begin; p + target_len <= attr.end; ++p) {

    if (memcmp(p, target, target_len) != 0) {
      continue;
    }

    // prevent false matches on fields such as transcript_gene_id
    const bool valid_left = (p == attr.begin) || is_attr_separator(*(p - 1));

    const bool valid_right =
        (p + target_len == attr.end) || is_attr_separator(*(p + target_len));

    if (valid_left && valid_right) {
      hit = p;
      break;
    }
  }

  if (hit == nullptr) {
    return false;
  }

  const char *quote1 = static_cast<const char *>(
      memchr(hit, '"', static_cast<size_t>(attr.end - hit)));

  if (quote1 == nullptr) {
    return false;
  }

  const char *quote2 = static_cast<const char *>(
      memchr(quote1 + 1, '"', static_cast<size_t>(attr.end - quote1 - 1)));

  if (quote2 == nullptr) {
    return false;
  }

  gene_id.assign(quote1 + 1, static_cast<size_t>(quote2 - quote1 - 1));

  return !gene_id.empty();
}

// ============================================================
// Process a single GTF record
// Does not build a vector<string>; only locates field positions
// ============================================================

inline void add_gtf_line_fast(const char *line, size_t line_size,
                              GroupMap &groups, bool canonical_chr_only) {

  if (line == nullptr || line_size == 0) {
    return;
  }

  if (line[0] == '#') {
    return;
  }

  // handle Windows line endings
  while (line_size > 0 &&
         (line[line_size - 1] == '\r' || line[line_size - 1] == '\n')) {
    --line_size;
  }

  if (line_size == 0) {
    return;
  }

  const char *line_begin = line;
  const char *line_end = line + line_size;
  const char *p = line_begin;

  Field fields[9];

  // the first 8 GTF columns are tab-separated
  for (int i = 0; i < 9; ++i) {

    if (i < 8) {

      const char *tab = static_cast<const char *>(
          memchr(p, '\t', static_cast<size_t>(line_end - p)));

      if (tab == nullptr) {
        return;
      }

      fields[i] = {p, tab};
      p = tab + 1;

    } else {
      fields[i] = {p, line_end};
    }
  }

  // column 3: feature
  const Field &feature = fields[2];

  if (feature.size() != 4 || memcmp(feature.begin, "exon", 4) != 0) {
    return;
  }

  // column 1: chromosome
  const Field &chr = fields[0];

  if (canonical_chr_only && !valid_chr_fast(chr)) {
    return;
  }

  // columns 4 and 5: start, end
  int start;
  int end;

  if (!parse_int_fast(fields[3], start)) {
    return;
  }

  if (!parse_int_fast(fields[4], end)) {
    return;
  }

  if (start <= 0 || end < start) {
    return;
  }

  // column 9: attributes
  string gene_id;

  if (!get_gene_id_fast(fields[8], gene_id)) {
    return;
  }

  // column 7: strand
  const Field &strand = fields[6];

  // key = gene_id + chromosome + strand
  string key;

  key.reserve(gene_id.size() + chr.size() + strand.size() + 2);

  key.append(gene_id);
  key.push_back('\t');
  key.append(chr.begin, chr.size());
  key.push_back('\t');
  key.append(strand.begin, strand.size());

  groups[key].push_back({start, end});
}

// ============================================================
// Internal read function
//
// Not exported to R.
// Reads and processes line by line; does not store all GTF lines.
// ============================================================

void read_gtf_lines(const string &path, GroupMap &groups, bool verbose,
                    bool canonical_chr_only) {

  LineReader reader(path);
  string line;

  size_t n = 0;

  while (reader.next(line)) {

    add_gtf_line_fast(line.data(), line.size(), groups, canonical_chr_only);

    ++n;

    if ((n & 131071) == 0) {
      checkUserInterrupt();
    }
  }

  if (verbose) {
    cli_emit("success", "Read " + std::to_string(n) + " GTF lines from file");
  }
}

// ============================================================
// Merge intervals and return a named numeric vector
// ============================================================

NumericVector finalize_gene_lengths(GroupMap &groups, bool verbose) {

  unordered_map<string, long long> gene_lengths;

  gene_lengths.reserve(groups.size());

  for (auto &item : groups) {

    vector<Interval> &intervals = item.second;

    if (intervals.empty()) {
      continue;
    }

    sort(intervals.begin(), intervals.end(),
         [](const Interval &a, const Interval &b) {
           if (a.start != b.start) {
             return a.start < b.start;
           }

           return a.end < b.end;
         });

    int current_start = intervals[0].start;
    int current_end = intervals[0].end;

    long long total = 0;

    for (size_t i = 1; i < intervals.size(); ++i) {

      const Interval &current = intervals[i];

      // merge overlapping intervals
      // whether adjacent intervals merge does not affect the final length
      if (current.start <= current_end) {

        if (current.end > current_end) {
          current_end = current.end;
        }

      } else {

        total += static_cast<long long>(current_end - current_start + 1);

        current_start = current.start;
        current_end = current.end;
      }
    }

    // add the final interval
    total += static_cast<long long>(current_end - current_start + 1);

    // key format is gene_id\tchr\tstrand
    const size_t tab = item.first.find('\t');

    if (tab == string::npos) {
      continue;
    }

    const string gene_id = item.first.substr(0, tab);

    gene_lengths[gene_id] += total;
  }

  NumericVector result(gene_lengths.size());
  CharacterVector names(gene_lengths.size());

  R_xlen_t i = 0;

  for (const auto &item : gene_lengths) {

    names[i] = item.first;
    result[i] = static_cast<double>(item.second);

    ++i;
  }

  result.attr("names") = names;

  if (verbose) {
    cli_emit("success", "Computed merged gene lengths for " +
                            std::to_string(result.size()) + " genes");
  }

  return result;
}

// ============================================================
// Compute gene length from a character vector passed from R
// ============================================================

NumericVector gene_length_from_r_lines(const CharacterVector &lines,
                                       bool verbose, bool canonical_chr_only) {

  const R_xlen_t n = lines.size();

  if (verbose) {
    cli_emit("info",
             "Processing " + std::to_string(n) + " GTF lines passed from R");
  }

  GroupMap groups;

  groups.reserve(65536);

  for (R_xlen_t i = 0; i < n; ++i) {

    SEXP element = STRING_ELT(lines, i);

    if (element == NA_STRING) {
      continue;
    }

    const char *line = CHAR(element);
    const size_t line_size = strlen(line);

    add_gtf_line_fast(line, line_size, groups, canonical_chr_only);

    if ((i & 131071) == 0) {
      checkUserInterrupt();
    }
  }

  return finalize_gene_lengths(groups, verbose);
}

// ============================================================
// R interface 1: pass a GTF file path
//
// This interface reads line by line directly and does not load
// the entire GTF file into memory.
// Suitable for:
//   gtf_file_to_gene_length("Homo_sapiens.GRCh37.75.gtf.gz")
// ============================================================

// [[Rcpp::export]]
NumericVector gtf_file_to_gene_length(std::string path, bool verbose = true,
                                      bool canonical_chr_only = false) {

  GroupMap groups;

  groups.reserve(65536);

  // read and process line by line internally
  read_gtf_lines(path, groups, verbose, canonical_chr_only);

  return finalize_gene_lengths(groups, verbose);
}

// ============================================================
// R interface 2: pass a character vector from R
//
// Suitable for:
//   x <- readLines("file.gtf")
//   r_gtf_lines_to_gene_length(x)
//
// Alternatively, pass an already-filtered GTF character vector.
// ============================================================

// [[Rcpp::export]]
NumericVector r_gtf_lines_to_gene_length(CharacterVector lines,
                                         bool verbose = true,
                                         bool canonical_chr_only = false) {

  return gene_length_from_r_lines(lines, verbose, canonical_chr_only);
}