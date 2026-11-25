    library(SigBridgeR)

# Extend Cell Type Annotation Methods

By default, SigBridgeR includes three built-in cell type annotation
algorithms: **SingleR**, **CellTypist**, and **mLLMCelltype**-all of
which require additional installation and dependencies (See
[README](https://wanglabcsu.github.io/SigBridgeR/index.html)). They can
all be invoked through the unified interface function `SCAnnotate`.

Extending cell annotation method is easy:

    my_method <- function(seurat_obj, ...) {
      # TODO: Annotate cell type
      return(seurat_obj)
    }
    my_sc_anno <- AnnotationMethod(
      method_name = "my_method",
      executor = my_method
    )
    Register(
      my_method = my_sc_anno
    )

To check:

    names(SCAnnotateStrategy)

    InterceptStrategy("SCAnnotateStrategy")

Here are some recommendations for user-defined cell annotation
functions.

1.  The first input and the first output should both be Seurat objects.
2.  If multiple objects need to be returned, it is recommended to use a
    `list` format. And the first element of the list should be the
    Seurat object.

## Session Info
