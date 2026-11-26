# SigBridgeR: Integrative Framework and Toolkit for Single-Cell Screening of Phenotype-Associated Cells

SigBridgeR is an integrative toolkit designed to identify
phenotype-associated cell subpopulations by combining phenotype(e.g.
survival, drug sensitivity), bulk expression and single-cell RNA-seq
data. It leverages multiple algorithms to robustly link cell features
with clinical or functional phenotypes. The package provides a unified
pipeline for cross-modal data analysis, enabling the discovery of
biologically and clinically relevant cell states in heterogeneous
samples.

## Main functions

The package includes multiple algorithms for integrative analysis of
single-cell and bulk data to identify phenotype-associated cell
populations. (see function
[`Screen`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md))

## Data requirements

- Single-cell RNA-seq data (Seurat object format)

- Bulk expression data

- Phenotype data (survival, drug sensitivity, etc.)

## Tutorial

For a detailed tutorial, please visit:
<https://wanglabcsu.github.io/SigBridgeR/>

## Citation

Use `citation("SigBridgeR")` to cite SigBridgeR in publications.

## License

SigBridgeR is licensed under the GPL version 3.0. Please see the LICENSE
file for details.

## See also

Useful links:

- <https://wanglabcsu.github.io/sigbridger>

- Report bugs at <https://github.com/WangLabCSU/SigBridgeR/issues>

## Author

Yuxi Yang <15364051195@163.com> (ORCID: 0009-0006-1329-1224) Shixiang
Wang <shixiang1994wang@gmail.com> (ORCID: 0000-0001-9855-7357) Zeyu Yan
<3191647388@qq.com> Duanchen Sun <sund@ohsu.edu> (ORCID: NA) Travis
Johnson <johnstrs@iu.edu> (ORCID: 0000-0002-4628-2256) Zhi Huang
<huang898@purdue.edu> (ORCID: NA) Qinran Zhang <qr_zhang@whu.edu.cn>
(ORCID: 0000-0001-8920-586X) Hongmei Zhang <zhanghmmail@163.com> (ORCID:
0009-0003-0914-1539) Dailin Gan <dgan@nd.edu> (ORCID: NA) Jun Li (ORCID:
0000-0003-4353-5761) Xinjia Ruan <ruan.cpu@foxmail.com> (ORCID:
0000-0003-2309-2641) Aimin Xie <aiminyy1993@gmail.com> (ORCID: NA)
Xiaosheng Wang <xiaosheng.wang@cpu.edu.cn> (ORCID: 0000-0002-7199-7093)
Yasmin Jolasun <yasmin.jolasun@mail.mcgill.ca> (ORCID:
0009-0009-2559-1015) Jun Ding <jun.ding@mcgill.ca> (ORCID:
0000-0001-5183-6885) Yuxiang Lin (ORCID: 0000-0001-8854-2995)
