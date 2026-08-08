"""
Annotate cells using CellTypist.

This module is intended to be called from R via reticulate.
"""

from __future__ import annotations

import datetime
import os
from pathlib import Path
from typing import Literal

import celltypist
from anndata import AnnData
from celltypist import models

LogLevel = Literal["info", "success", "warning", "error", "debug"]


def get_dirname(path_str: str) -> str:
    """
    获取路径的 dirname（目录名）

    Returns:
        str: 路径中的目录部分，纯文件名则返回空字符串
    """
    path = Path(path_str)

    # 如果是纯文件名（没有父目录），返回空字符串
    if path.parent == Path("."):
        return ""

    return str(path.parent)


def is_path(input_str: str) -> bool:
    """
    检查输入字符串是否是文件路径

    Returns:
        bool: True 如果是路径，False 如果是纯名称
    """
    path = Path(input_str)
    return (
        path.is_absolute()
        or "/" in input_str
        or "\\" in input_str
        or (len(input_str) >= 2 and input_str[1] == ":")
    )


def ts_print(
    message: str,
    symbol: LogLevel | None = None,
    color: bool = True,
) -> None:
    """
    带时间戳的信息输出函数

    Args:
        message: 要输出的消息内容
        symbol: CLI 符号类型，可选值: 'info', 'success', 'warning', 'error', 'debug'，默认无符号
        color: 是否启用 ANSI 颜色输出（默认 True）
    """
    timestamp = datetime.datetime.now(tz=datetime.timezone.utc).strftime(
        "[%Y/%m/%d %H:%M:%S]"
    )

    symbols = {
        "info": ("ℹ", "\033[36m" if color else ""),  # Cyan
        "success": ("✔", "\033[32m" if color else ""),  # Green
        "warning": ("⚠", "\033[33m" if color else ""),  # Yellow
        "error": ("✖", "\033[31m" if color else ""),  # Red
        "debug": ("◼", "\033[35m" if color else ""),  # Magenta
    }

    t_prefix: str = f"{timestamp} "
    if symbol in symbols:
        sym, col = symbols[symbol]
        symbol_prefix: str = f"{col}{sym} \033[0m" if color else f"{sym} "
    else:
        symbol_prefix: str = ""

    print(symbol_prefix + f"{t_prefix}{message}")


def configure_model_path(model: str) -> str:
    """
    Configure CellTypist's model directory when a model path is supplied.

    Returns
    -------
    str
        The model filename used by CellTypist.
    """
    model_path = Path(model)

    if not is_path(model):
        return model

    if not model_path.exists():
        raise FileNotFoundError(f"CellTypist model was not found: {model}")

    os.environ["CELLTYPIST_FOLDER"] = str(model_path.parent)
    return model_path.name


def annotate_celltypist(
    adata: AnnData,
    model: str,
    *,
    transpose_input: bool = False,
    gene_file: str | None = None,
    cell_file: str | None = None,
    majority_voting: bool = False,
    over_clustering: (
        str
        | list
        | tuple
        | None
        # numpy.ndarray,
        # pandas.core.series.Series,
        # pandas.core.indexes.base.Index,
        # NoneType,
    ) = None,
    use_GPU: bool = False,
    mode: Literal["best match", "prob match"] = "best match",
    p_thres: float = 0.5,
    min_prop: float = 0,
    download: bool = False,
    force_update: bool = False,
    verbose: bool = True,
):
    """
    Annotate cells in an AnnData object using CellTypist.

    Parameters
    ----------
    adata
        Input data as an :class:`~anndata.AnnData` object already loaded in
        memory. Genes should be gene symbols; non-expressed genes are preferred
        to be provided as well.
    model
        Model used to predict the input cells. Can be a
        :class:`~celltypist.models.Model` object that wraps the logistic
        Classifier and the StandardScaler, the path to the desired model file,
        or a built-in model name (e.g. ``'Immune_All_Low'``). To see all
        available models and their descriptions, use
        :func:`~celltypist.models.models_description`.
    transpose_input
        Whether to transpose the input matrix. Set to ``True`` if the input is
        provided in a gene-by-cell format. (Default: ``False``)
    gene_file
        Path to the file which stores each gene per line, corresponding to the
        genes used in the provided mtx file. Ignored if the input is not
        provided in the mtx format. (Default: ``None``)
    cell_file
        Path to the file which stores each cell per line, corresponding to the
        cells used in the provided mtx file. Ignored if the input is not
        provided in the mtx format. (Default: ``None``)
    majority_voting
        Whether to refine the predicted labels by running the majority voting
        classifier after over-clustering. (Default: ``False``)
    over_clustering
        Over-clustering result, which can be provided in several ways:
        1) an input plain file with the over-clustering result of one cell per
        line; 2) a string key specifying an existing metadata column in the
        AnnData (pre-created by the user); 3) a python list, tuple, numpy
        array, pandas series or index representing the over-clustering result
        of the input cells; 4) if none of the above is provided, a heuristic
        over-clustering approach will be used according to the size of the
        input data. Ignored if ``majority_voting`` is ``False``.
        (Default: ``None``)
    use_GPU
        Whether to use GPU for over-clustering on the basis of
        rapids-singlecell. This argument is only relevant when
        ``majority_voting = True``. (Default: ``False``)
    mode
        The way cell prediction is performed. For each query cell, the default
        (``'best match'``) is to choose the cell type with the largest
        score/probability as the final prediction. Setting to ``'prob match'``
        will enable a multi-label classification, which assigns 0 (i.e.,
        unassigned), 1, or >=2 cell type labels to each query cell.
        (Default: ``'best match'``)
    p_thres
        Probability threshold for the multi-label classification. Ignored if
        ``mode`` is ``'best match'``. (Default: ``0.5``)
    min_prop
        For the dominant cell type within a subcluster, the minimum proportion
        of cells required to support naming of the subcluster by this cell
        type. Ignored if ``majority_voting`` is ``False``. Subclusters that
        fail to pass this proportion threshold will be assigned
        ``'Heterogeneous'``. (Default: ``0``)
    download
        Whether to download the CellTypist model files first (e.g. when the
        model is not available locally). (Default: ``False``)
    force_update
        Whether to force re-downloading the CellTypist model files.
        (Default: ``False``)
    verbose
        Whether to print progress messages. (Default: ``True``)

    Returns
    -------
    tuple
        ``predicted_labels``, ``decision_matrix`` and ``probability_matrix``.
    """
    if verbose:
        ts_print(f"celltypist version: {celltypist.__version__}", symbol="info")

    if not isinstance(adata, AnnData):
        raise TypeError(f"adata must be an AnnData object, got {type(adata)!r}")

    if not isinstance(model, str) or not model.strip():
        raise TypeError("model must be a non-empty string")

    model_name: str = configure_model_path(model)

    if download:
        if verbose:
            ts_print(
                f"Downloading CellTypist models: {model_name}",
                symbol="info",
            )

        models.download_models(force_update=force_update, model=model_name)

        if verbose:
            ts_print(
                f"Models saved to `{models.models_path}`",
                symbol="success",
            )

    loaded_model = models.Model.load(model=model_name)  # celltypist.models.Model

    if verbose:
        ts_print(
            f"Annotating cells with model `{model_name}`",
            symbol="info",
        )

    predictions = celltypist.annotate(
        adata,
        model=loaded_model,
        transpose_input=transpose_input,
        gene_file=gene_file,
        cell_file=cell_file,
        mode=mode,
        p_thres=p_thres,
        majority_voting=majority_voting,
        over_clustering=over_clustering,
        use_GPU=use_GPU,
        min_prop=min_prop,
    )

    if verbose:
        ts_print(
            "CellTypist annotation completed.",
            symbol="success",
        )

    return (
        predictions.predicted_labels,
        predictions.decision_matrix,
        predictions.probability_matrix,
    )


if __name__ == "__main__":
    raise SystemExit("This module is intended to be called from R via reticulate.")
# R 端模块
