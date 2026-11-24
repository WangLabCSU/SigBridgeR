---
name: reference-coding
description: Retrieve how peer projects solved a problem before writing code in SigBridgeR. Query the local XERJ corpora (this package, the 12 upstream screening algorithms, Seurat, Rcpp/sparse-matrix code) with xerj code / xerj search, read the cited file:line, respect the licence banner, then implement. Use when adding or changing a screening method, adapting an upstream algorithm, touching Seurat object manipulation, or writing C++ in src/.
---

# Reference coding with XERJ

## Purpose

SigBridgeR re-implements and adapts a dozen published screening algorithms. Before
inventing logic, retrieve how the upstream project (or Seurat) already solved it,
cite `file:line`, and write our own code against the project's conventions.

XERJ is optional developer tooling: nothing here is needed to build or check the
package. Commands assume `xerj` is on `PATH` and a local node answers on
`http://localhost:9200` (see Bootstrap).

## Corpora

| Need | Command | Covers |
| --- | --- | --- |
| This package's own code | `xerj search --prefix sigbridger "<words>"` / `xerj def --prefix sigbridger <symbol>` | `R/`, `src/`, `tests/`, `vignettes/`, `man/` |
| Upstream screening algorithms | `xerj code scpheno "<what you need>"` | Scissor, scPAS, scAB, ScPP, DEGAS, DEGASv2, LP_SGL, PIPET, SIDISH, SCIPAC, TiRank, PhenoMapR |
| Seurat object and assay APIs | `xerj code seurat "<what you need>"` | satijalab/seurat, satijalab/seurat-object |
| C++ / sparse-matrix plumbing | `xerj code rcpp-sc "<what you need>"` | RcppArmadillo, assorthead, beachmat |
| Own-lab utility package | `xerj code sigbridger-utils "<what you need>"` | WangLabCSU/SigBridgeRUtils |

Add `-k 8` for more passages, `--lang r` or `--lang cpp` to filter, `--full 1500`
for longer context. `xerj code` prints a licence header per hit; read it. Corpora
are pinned shallow clones under `~/.xerj-code/corpora/<name>/`; open the cited
file there when the passage is not enough.

## Bootstrap (once per machine)

```bash
xerj --insecure --data-dir ~/xerj-data &        # data dir must be outside the repo
xerj autoindex . --prefix sigbridger --state-dir ~/.xerj/autoindex/sigbridger

xerj corpus add scpheno \
  https://github.com/sunduanchen/Scissor https://github.com/aiminXie/scPAS \
  https://github.com/jinworks/scAB https://github.com/WangX-Lab/ScPP \
  https://github.com/tsteelejohnson91/DEGAS https://github.com/ElliotLiu1997/DEGASv2 \
  https://github.com/hongmeizhanghm/LP_SGL https://github.com/ruan2ruan/PIPET \
  https://github.com/mcgilldinglab/SIDISH https://github.com/mianaz/SCIPAC \
  https://github.com/LenisLin/TiRank https://github.com/brooksbenard/PhenoMapR
xerj corpus add seurat https://github.com/satijalab/seurat https://github.com/satijalab/seurat-object
xerj corpus add rcpp-sc https://github.com/RcppCore/RcppArmadillo \
  https://github.com/LTLA/assorthead https://github.com/LTLA/beachmat
xerj corpus add sigbridger-utils https://github.com/WangLabCSU/SigBridgeRUtils --own

xerj corpus index scpheno && xerj corpus index seurat
xerj corpus index rcpp-sc && xerj corpus index sigbridger-utils
```

`xerj corpus add` shells out to `git clone`, so a firewall or proxy is configured
the git way: `export https_proxy=http://127.0.0.1:7890` before the command (a
local proxy can be 90x faster than a throttled direct GitHub connection).

If `xerj code <corpus>` says the corpus is unknown, the Bootstrap step for it was
not run; fall back to `xerj search --prefix sigbridger` and say so.

## Workflow

1. State the mechanism you need, in the words the source would use — e.g.
   "network penalty elastic net survival", "rank matrix correlation bulk single cell",
   "sparse column mean by group". Retrieval is lexical: it matches the identifiers
   and terms you type, not meaning.
2. Run the matching `xerj code` (upstream) and `xerj search --prefix sigbridger`
   (our side, to see whether it already exists).
3. Open the cited `file:line` in `~/.xerj-code/corpora/<corpus>/` when needed.
4. Decide and implement in this repo's style: S7 classes, `rlang` errors,
   `ScreenMethod` contract, `check_*` validators.
5. If retrieval returns nothing relevant, say so and proceed normally — do not
   pretend a citation exists.

## Licence rule

- `GPL`/`AGPL` (most of `scpheno`): adapt the **approach**, never copy code.
  SigBridgeR is GPL-3, but vendored code needs its author and `Copyright` entry in
  `DESCRIPTION` — see the existing `R/vendor-*.R` files for the pattern.
- `MIT`/`BSD`/`Apache`: reuse allowed with attribution.
- `NONE-FOUND`/`UNKNOWN`: approach only.
- `approach-only licence ... read, do not copy` in a hit header is binding.

## Refresh

```bash
# upstream corpora: pull, re-pin with the Bootstrap command, re-index
for d in ~/.xerj-code/corpora/scpheno/*/; do git -C "$d" pull --ff-only; done
xerj corpus add scpheno <same git URLs as Bootstrap>
xerj corpus index scpheno

# this package
xerj autoindex . --prefix sigbridger --state-dir ~/.xerj/autoindex/sigbridger
```

Corpus registry: `xerj corpus list`.
