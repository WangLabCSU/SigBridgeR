import scanpy as sc
import os
import celltypist
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import logging
import traceback


def SetupLogging():
    """Set the log format"""
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "[%(asctime)s] %(message)s", datefmt="%Y/%m/%d %H:%M:%S"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


def CreateDir(*dirs):
    """Create the required directories"""
    for d in dirs:
        os.makedirs(d, exist_ok=True)


def DrawDotplot(predictions, plot_dir):
    """Generate and save the dot plot"""
    try:
        dotplot = celltypist.dotplot(
            predictions,
            use_as_reference="seurat_clusters",
            use_as_prediction="predicted_labels",
            figsize=(10, 5),
            return_fig=True,
            show=False,
        )
        dotplot_path = os.path.join(plot_dir, "celltypist_dotplot.pdf")
        dotplot.savefig(dotplot_path, bbox_inches="tight")
        return dotplot_path
    except Exception as e:
        logger.error(f"Error generating dotplot: {str(e)}")
        raise


def DrawUMAP(adata, plot_dir):
    """Generate and save the UMAP"""
    try:
        fig = plt.subplots(figsize=(8, 6))
        sc.pl.umap(
            adata,
            color=["clusterName", "majority_voting"],
            legend_loc="on data",
            legend_fontsize=8,
            show=False,
        )
        umap_path = os.path.join(plot_dir, "celltypist_umap.pdf")
        fig.savefig(umap_path, bbox_inches="tight", dpi=300)
        return umap_path
    except Exception as e:
        logger.error(f"Error generating UMAP: {str(e)}")
        raise


def Main():
    global path, model_path, data_output_dir, plot_output_dir, n_comps_required

    # Verify the necessary parameters
    required_vars = ["path", "model_path"]
    for var in required_vars:
        if var not in globals():
            raise ValueError(f"Required argument '{var}' not found")

    # Set arguments
    data_output_dir = data_output_dir if "data_output_dir" in globals() else "."
    plot_output_dir = plot_output_dir if "plot_output_dir" in globals() else "."
    n_comps_required = n_comps_required if "n_comps_required" in globals() else 50

    CreateDir(data_output_dir, plot_output_dir)

    try:
        # 1. Load data
        logger.info("Loading data...")
        adata = sc.read_h5ad(path)

        # 2. Preprocess
        logger.info("Normalizing data...")
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # 3. Check and run PCA
        if "X_pca" not in adata.obsm or adata.obsm["X_pca"].shape[1] < n_comps_required:
            logger.info(f"Computing PCA with {n_comps_required} components...")
            sc.pp.pca(adata, n_comps=n_comps_required)
        logger.info(f"PCA shape: {adata.obsm['X_pca'].shape}")

        # 4. CellTypist annotation
        logger.info("Running CellTypist annotation...")
        predictions = celltypist.annotate(adata, model=model_path, majority_voting=True)

        # 5. optional: dotplot
        dotplot_path = DrawDotplot(predictions, plot_output_dir)
        logger.info(f"Dotplot saved to {dotplot_path}")

        # 6. Transform the prediction results
        adata = predictions.to_adata()
        adata.obs["clusterName"] = adata.obs["majority_voting"]

        # 7. optional: UMAP
        umap_path = DrawUMAP(adata, plot_output_dir)
        logger.info(f"UMAP plot saved to {umap_path}")

        # 8. output
        logger.info("Preparing output data...")
        if "cellID" not in adata.obs.columns:
            adata.obs["cellID"] = adata.obs.index.astype(str)

        cluster_colors = adata.uns.get("clusterName_colors", [])
        cluster_hex = (
            [mcolors.rgb2hex(c) for c in cluster_colors]
            if cluster_colors
            else ["#000000"] * len(adata.obs["clusterName"].cat.categories)
        )

        color_df = pd.DataFrame(
            {
                "clusterName": adata.obs["clusterName"].cat.categories,
                "clusterColor": cluster_hex,
            }
        )

        result_df = pd.merge(
            adata.obs[["cellID", "clusterName", "predicted_labels", "conf_score"]],
            color_df,
            on="clusterName",
            how="left",
        )

        # 9. Save
        output_path = os.path.join(data_output_dir, "celltypist.annotation.csv")
        result_df.to_csv(output_path, index=False)
        logger.info(f"Annotation results (for Seurat meta.data) saved to {output_path}")

        logger.info("\033[92mCellTypist analysis completed successfully\033[0m")

    except Exception as e:
        logger.error(f"\033[91mError occurred: {str(e)}\033[0m")
        logger.error(traceback.format_exc())
        raise


# Initialize log
logger = SetupLogging()

if __name__ == "__main__":
    Main()
