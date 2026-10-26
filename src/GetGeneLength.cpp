#include <Rcpp.h>
#include <zlib.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "SigBridgeR.h"

using namespace Rcpp;

namespace {

struct Interval {
  std::int64_t start;
  std::int64_t end;
};

struct Field {
  const char *begin;
  const char *end;
  std::size_t size() const { return static_cast<std::size_t>(end - begin); }
};

using GroupMap = std::unordered_map<std::string, std::vector<Interval>>;

inline void r_warning(const std::string &msg) {
  Rcpp::Function w("warning");
  w(msg, Rcpp::Named("call.") = false);
}

// ---------------------------------------------------------------
// zlib 透明读取 .gz / 纯文本；gzread 分块 + memchr 切行（零拷贝）
// next() 返回的视图在下一次 next() 之前有效，不含 '\n'
// ---------------------------------------------------------------
class LineReader {
public:
  explicit LineReader(const std::string &path)
      : fp_(gzopen(path.c_str(), "rb")), buf_(1 << 20) {
    if (fp_ == nullptr) {
      stop("Cannot open file: " + path);
    }
    gzbuffer(fp_, 1 << 18);
  }
  ~LineReader() {
    if (fp_ != nullptr)
      gzclose(fp_);
  }
  LineReader(const LineReader &) = delete;
  LineReader &operator=(const LineReader &) = delete;

  bool next(const char *&line, std::size_t &len) {
    while (true) {
      const char *base = buf_.data();
      const char *nl = static_cast<const char *>(
          std::memchr(base + pos_, '\n', end_ - pos_));
      if (nl != nullptr) {
        line = base + pos_;
        len = static_cast<std::size_t>(nl - line);
        pos_ = static_cast<std::size_t>(nl - base) + 1;
        return true;
      }
      if (eof_) {
        if (pos_ == end_)
          return false;
        line = base + pos_;
        len = end_ - pos_; // 最后一行没有换行符
        pos_ = end_;
        return true;
      }
      refill();
    }
  }

private:
  void refill() {
    const std::size_t remain = end_ - pos_;
    if (pos_ > 0) {
      if (remain > 0)
        std::memmove(buf_.data(), buf_.data() + pos_, remain);
      pos_ = 0;
      end_ = remain;
    }
    if (end_ == buf_.size())
      buf_.resize(buf_.size() * 2); // 超长行
    const unsigned room = static_cast<unsigned>(
        std::min<std::size_t>(buf_.size() - end_, 1u << 30));
    const int got = gzread(fp_, buf_.data() + end_, room);
    if (got <= 0) {
      int errnum = Z_OK;
      const char *msg = gzerror(fp_, &errnum);
      if (got < 0 || (errnum != Z_OK && errnum != Z_STREAM_END)) {
        stop(std::string(
                 "Error while reading GTF (truncated or corrupted file?): ") +
             (msg != nullptr ? msg : "unknown"));
      }
      eof_ = true;
      return;
    }
    end_ += static_cast<std::size_t>(got);
  }

  gzFile fp_;
  std::vector<char> buf_;
  std::size_t pos_ = 0, end_ = 0;
  bool eof_ = false;
};

// ---------------------------------------------------------------
inline bool is_attr_separator(char c) {
  return c == ' ' || c == '\t' || c == ';';
}

// 按 `key value;` 逐项解析；值可带引号也可不带；只比较 key
inline bool get_gene_id_fast(const Field &attr, const char *&val,
                             std::size_t &val_len) {
  const char *p = attr.begin;
  const char *const e = attr.end;
  while (p < e) {
    while (p < e && is_attr_separator(*p))
      ++p;
    const char *key = p;
    while (p < e && !is_attr_separator(*p))
      ++p;
    const std::size_t key_len = static_cast<std::size_t>(p - key);
    while (p < e && (*p == ' ' || *p == '\t'))
      ++p;

    const char *v = p;
    std::size_t vl = 0;
    if (p < e && *p == '"') {
      v = ++p;
      while (p < e && *p != '"')
        ++p;
      if (p == e)
        return false; // 引号未闭合
      vl = static_cast<std::size_t>(p - v);
      ++p;
    } else {
      while (p < e && *p != ';')
        ++p;
      vl = static_cast<std::size_t>(p - v);
      while (vl > 0 && (v[vl - 1] == ' ' || v[vl - 1] == '\t'))
        --vl;
    }
    if (key_len == 7 && std::memcmp(key, "gene_id", 7) == 0) {
      val = v;
      val_len = vl;
      return vl > 0;
    }
  }
  return false;
}

// 保留 1-22, X, Y, MT/M；兼容 "chr" 前缀；拒绝前导零
inline bool valid_chr_fast(Field chr) {
  if (chr.size() > 3 && std::memcmp(chr.begin, "chr", 3) == 0)
    chr.begin += 3;
  const std::size_t n = chr.size();
  if (n == 0 || n > 2)
    return false;
  if (n == 1 &&
      (chr.begin[0] == 'X' || chr.begin[0] == 'Y' || chr.begin[0] == 'M'))
    return true;
  if (n == 2 && chr.begin[0] == 'M' && chr.begin[1] == 'T')
    return true;
  if (chr.begin[0] == '0')
    return false;
  int value = 0;
  for (const char *p = chr.begin; p != chr.end; ++p) {
    if (*p < '0' || *p > '9')
      return false;
    value = value * 10 + (*p - '0');
  }
  return value >= 1 && value <= 22;
}

inline bool parse_coord(const Field &f, std::int64_t &out) {
  const std::size_t n = f.size();
  if (n == 0 || n > 15)
    return false; // 15 位以内，int64 不会溢出
  std::int64_t x = 0;
  for (const char *p = f.begin; p != f.end; ++p) {
    if (*p < '0' || *p > '9')
      return false;
    x = x * 10 + (*p - '0');
  }
  out = x;
  return true;
}

// ---------------------------------------------------------------
struct GtfAccumulator {
  explicit GtfAccumulator(bool canonical_chr_only)
      : canonical(canonical_chr_only) {
    groups.reserve(65536);
  }
  GtfAccumulator(const GtfAccumulator &) = delete;
  GtfAccumulator &operator=(const GtfAccumulator &) = delete;

  void add_line(const char *line, std::size_t n) {
    while (n > 0 && (line[n - 1] == '\r' || line[n - 1] == '\n'))
      --n;
    if (n == 0 || line[0] == '#')
      return;

    Field f[9];
    const char *p = line;
    const char *const e = line + n;
    for (int i = 0; i < 8; ++i) {
      const char *tab = static_cast<const char *>(
          std::memchr(p, '\t', static_cast<std::size_t>(e - p)));
      if (tab == nullptr)
        return;
      f[i] = {p, tab};
      p = tab + 1;
    }
    f[8] = {p, e};

    if (f[2].size() != 4 || std::memcmp(f[2].begin, "exon", 4) != 0)
      return;
    if (canonical && !valid_chr_fast(f[0]))
      return;

    std::int64_t start = 0, end = 0;
    if (!parse_coord(f[3], start) || !parse_coord(f[4], end) || start <= 0 ||
        end < start) {
      ++n_malformed;
      return;
    }
    const char *gid = nullptr;
    std::size_t gid_len = 0;
    if (!get_gene_id_fast(f[8], gid, gid_len)) {
      ++n_malformed;
      return;
    }

    // key = gene_id \t chr \t strand（复用缓冲，无额外分配）
    key.clear();
    key.append(gid, gid_len);
    key.push_back('\t');
    key.append(f[0].begin, f[0].size());
    key.push_back('\t');
    key.append(f[6].begin, f[6].size());

    // 同一基因的 exon 通常连续：命中缓存则跳过哈希
    if (cur == nullptr || key != last_key) {
      cur = &groups[key]; // unordered_map 的元素引用在 rehash 后仍有效
      last_key.swap(key);
    }
    cur->push_back({start, end});
    ++n_exon;
  }

  bool canonical;
  GroupMap groups;
  std::size_t n_exon = 0;
  std::size_t n_malformed = 0;

private:
  std::string key, last_key;
  std::vector<Interval> *cur = nullptr;
};

// 合并区间 -> 按 gene_id 汇总 -> 按 gene_id 排序（输出顺序确定）
NumericVector finalize_gene_lengths(GroupMap &groups) {
  std::unordered_map<std::string, long long> sums;
  sums.reserve(groups.size());

  for (auto &item : groups) {
    std::vector<Interval> &iv = item.second;
    if (iv.empty())
      continue;
    std::sort(iv.begin(), iv.end(), [](const Interval &a, const Interval &b) {
      return a.start != b.start ? a.start < b.start : a.end < b.end;
    });

    long long total = 0;
    std::int64_t cs = iv[0].start, ce = iv[0].end;
    for (std::size_t i = 1; i < iv.size(); ++i) {
      if (iv[i].start <= ce) {
        if (iv[i].end > ce)
          ce = iv[i].end;
      } else {
        total += ce - cs + 1;
        cs = iv[i].start;
        ce = iv[i].end;
      }
    }
    total += ce - cs + 1;

    // 从尾部切 key，即使 gene_id 含 '\t' 也不会错位
    const std::string &key = item.first;
    const std::size_t p2 = key.rfind('\t');
    if (p2 == std::string::npos || p2 == 0)
      continue;
    const std::size_t p1 = key.rfind('\t', p2 - 1);
    if (p1 == std::string::npos)
      continue;
    sums[key.substr(0, p1)] += total;
  }

  std::vector<std::pair<std::string, long long>> items(sums.begin(),
                                                       sums.end());
  std::sort(items.begin(), items.end(),
            [](const std::pair<std::string, long long> &a,
               const std::pair<std::string, long long> &b) {
              return a.first < b.first;
            });

  const R_xlen_t n = static_cast<R_xlen_t>(items.size());
  NumericVector result(n);
  CharacterVector names(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    names[i] = items[i].first;
    result[i] = static_cast<double>(items[i].second);
  }
  result.attr("names") = names;
  return result;
}

NumericVector finish(GtfAccumulator &acc, bool verbose) {
  NumericVector res = finalize_gene_lengths(acc.groups);
  if (acc.n_malformed > 0) {
    r_warning(std::to_string(acc.n_malformed) +
              " exon lines with invalid coordinates or missing gene_id were "
              "skipped.");
  }
  if (res.size() == 0) {
    r_warning(
        "No exon with a valid gene_id was found. Check the file format and "
        "`canonical_chr_only` (chromosome names).");
  }
  if (verbose) {
    cli_emit("success", "Computed merged gene lengths for " +
                            std::to_string(res.size()) + " genes from " +
                            std::to_string(acc.n_exon) + " exons");
  }
  return res;
}

} // namespace

// [[Rcpp::export]]
NumericVector gtf_file_to_gene_length(const std::string &path,
                                      bool verbose = true,
                                      bool canonical_chr_only = false) {
  GtfAccumulator acc(canonical_chr_only);
  std::size_t n_lines = 0;
  {
    LineReader reader(path);
    const char *line = nullptr;
    std::size_t len = 0;
    while (reader.next(line, len)) {
      acc.add_line(line, len);
      if ((++n_lines & 0x1FFFF) == 0)
        checkUserInterrupt();
    }
  } // 先关闭文件，再做汇总
  if (verbose) {
    cli_emit("info", "Read " + std::to_string(n_lines) + " GTF lines");
  }
  return finish(acc, verbose);
}

// [[Rcpp::export]]
NumericVector r_gtf_lines_to_gene_length(const CharacterVector &lines,
                                         bool verbose = true,
                                         bool canonical_chr_only = false) {
  const R_xlen_t n = lines.size();
  GtfAccumulator acc(canonical_chr_only);
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP el = STRING_ELT(lines, i);
    if (el == NA_STRING)
      continue;
    acc.add_line(CHAR(el), static_cast<std::size_t>(Rf_length(el)));
    if (((i + 1) & 0x1FFFF) == 0)
      checkUserInterrupt();
  }
  return finish(acc, verbose);
}