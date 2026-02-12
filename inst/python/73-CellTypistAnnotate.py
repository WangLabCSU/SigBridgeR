# ----------------------------------------------------------------------------------------------------
# ! Don't run this script directly, use `R/73-CellTypistAnnotate.R` instead
# ----------------------------------------------------------------------------------------------------

from anndata import AnnData
import celltypist
from celltypist import models
from datetime import datetime
from typing import Literal, Optional, Callable, Any
import inspect
import os
import pandas as pd
from pathlib import Path
import re


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
    # 检查是否是绝对路径
    is_absolute = path.is_absolute()
    # 检查字符串中是否包含路径分隔符
    contains_separator = "/" in input_str or "\\" in input_str
    # 检查是否包含 Windows 盘符 (如 C:, D:)
    has_windows_drive = bool(re.match(r"^[A-Za-z]:", input_str))
    # 满足任一条件即为路径
    return is_absolute or contains_separator or has_windows_drive


def ts_print(
    message: str,
    symbol: Optional[Literal["info", "success", "warning", "error", "debug"]] = None,
    color: bool = True,
) -> None:
    """
    带时间戳的信息输出函数

    Args:
        message: 要输出的消息内容
        symbol: CLI 符号类型，可选值: 'info', 'success', 'warning', 'error', 'debug'，默认无符号
        color: 是否启用 ANSI 颜色输出（默认 True）
    """
    timestamp = datetime.now().strftime("[%Y/%m/%d %H:%M:%S]")

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


def filter_args_4_func(input_dict: dict[str, Any], func: Callable) -> dict[str, Any]:
    """
    从 input_dict 中提取出所有键名与 func 的参数名匹配的键值对。

    Parameters
    ----------
    input_dict : dict
        输入字典，包含可能的参数键值对。
    func : callable
        目标函数。
    Returns
    -------
    dict
        仅包含 func 所需参数的子集字典，可用于 **kwargs 调用。
    """
    sig = inspect.signature(func)
    func_params = set(sig.parameters.keys())
    return {k: v for k, v in input_dict.items() if k in func_params}


def main() -> pd.DataFrame:
    # * recept them from R script
    adata: AnnData = globals().get("adata")
    model: str = globals().get("model")
    verbose: bool = globals().get("verbose")
    download: bool = globals().get("download")
    force_update: bool = globals().get("force_update")

    if type(model) is not str:
        raise ValueError(
            f"Unknown type of model {type(model)}, expected <class 'str'>)"
        )

    if is_path(model):
        model: str = Path(model).name
        os.environ["CELLTYPIST_FOLDER"] = get_dirname(model)

    if download:
        if verbose:
            ts_print(message=f"Download CellTypist models: {model}", symbol="info")

        # ! tried but failed with native `celltypist.models.download_models` when
        # ! asked by R script
        models.download_models(force_update=force_update)

        if verbose:
            ts_print(
                message=f"Models saved to `{models.models_path}`", symbol="success"
            )

    model: celltypist.models.Model = models.Model.load(model=model)

    kwarg: dict[str, Any] = filter_args_4_func(
        input_dict=globals(), func=celltypist.annotate
    )

    predictions = celltypist.annotate(adata, **kwarg)

    return (
        # the predicted cell type labels
        predictions.predicted_labels,
        # the matrix representing the decision score of each cell belonging to a given cell type
        predictions.decision_matrix,
        # the matrix representing the probability each cell belongs to a given cell type
        # (transformed from decision matrix by the sigmoid function)
        predictions.probability_matrix,
    )


predicted_labels, decision_matrix, probability_matrix = main()
