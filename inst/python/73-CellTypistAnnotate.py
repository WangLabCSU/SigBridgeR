# -------------------------------------------------------------------------------------------------------
# ! Don't run this script directly, use `R/73-CellTypistAnnotate.R` instead
# -------------------------------------------------------------------------------------------------------

from anndata import AnnData
import celltypist
from celltypist import models
from datetime import datetime
from typing import Literal, Optional, Callable, Any
import inspect
import sys
import pandas as pd

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
        "info": ("💡", "\033[36m" if color else ""),  # Cyan
        "success": ("✅", "\033[32m" if color else ""),  # Green
        "warning": ("⚠️", "\033[33m" if color else ""),  # Yellow
        "error": ("❌", "\033[31m" if color else ""),  # Red
        "debug": ("🐞", "\033[35m" if color else ""),  # Magenta
    }

    prefix = f"{timestamp} "
    if symbol in symbols:
        sym, col = symbols[symbol]
        prefix += f"{col}{sym} \033[0m" if color else f"{sym} "

    print(f"{prefix}{message}" if not color else f"{prefix}\033[0m{message}")


def filter_args_4_func(dict: dict[str, Any], func: Callable) -> dict[str, Any]:
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
    return {k: v for k, v in dict.items() if k in func_params}


def main() -> pd.DataFrame:
    # * recept from R script
    model: str = globals().get("model")
    adata: AnnData = globals().get("adata")
    verbose: bool = globals().get("verbose")
    download: bool = globals().get("download")
    kwarg = filter_args_4_func(dict=globals(), func=celltypist.annotate)

    if download:
        if verbose:
            ts_print(message="Download CellTypist models", symbol="info")
        models.download_models(model=model, force_update=True)
        if verbose:
            ts_print(message=f"Models saved to `{models.models_path}`")

    if type(model) is str:
        predictions = celltypist.annotate(adata, model=model, **kwarg)
        return predictions.predicted_labels
    else:
        sys.exit(f"Unknown type of model {type(model)}\n> Expected <class 'str'>")

res: pd.DataFrame = main()
