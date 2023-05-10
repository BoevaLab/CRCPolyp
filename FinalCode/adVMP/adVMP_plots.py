import pandas as pd
import numpy as np
import pathlib as pl
import matplotlib.pyplot as plt
import os

from typing import Optional, List, Dict
from sklearn.metrics import RocCurveDisplay, PrecisionRecallDisplay
from statannotations.Annotator import Annotator

import seaborn as sns

from utils.plotting import transform_plot_ax
from .adVMP_discovery import get_polyp_size_nr_link

from typing import Optional, List, Dict
import os


def get_performance_plots(
    heatmap_df: pd.DataFrame,
    fig_dir: pl.Path,
    hue_palette_worm: Dict,
    ftsize: int = 10,
    leg_ftsize: int = 10,
    hit_lim: int = 4,
    hue_worm: str = "Ad",
    figsize: int = 12,
    hue_order: Optional[List] = None,
    rocauc: bool = True,
    leg_title: str = "Adenoma",
) -> None:
    os.makedirs(fig_dir, exist_ok=True)
    cols = heatmap_df.columns.str.startswith("cg") + heatmap_df.columns.str.startswith(
        "chr"
    )
    ax = sns.clustermap(
        heatmap_df.loc[:, cols].clip(hit_lim, -hit_lim).T,
        cmap="vlag",
        center=0,
        mask=heatmap_df.loc[:, cols].abs().T < hit_lim,
        col_cluster=False,
    )
    ax.ax_heatmap.set_xticks([])
    ax.ax_heatmap.set_yticks([])
    ax.ax_heatmap.set_xlabel("Samples", fontsize=ftsize)
    ax.ax_heatmap.set_ylabel("adVMPs", fontsize=ftsize)
    ax.cax.set_visible(False)
    ax.figure.savefig(fig_dir / "heatmap_probes.png", dpi=250, bbox_inches="tight")

    fig, ax = plt.subplots(1, 1, figsize=(3 * figsize, figsize))
    sns.scatterplot(
        data=heatmap_df,
        x="Order",
        y="Hit fraction",
        hue=hue_worm,
        hue_order=hue_order,
        palette=hue_palette_worm,
    )
    transform_plot_ax(ax, legend_title=leg_title, ftsize=ftsize, leg_ftsize=leg_ftsize)
    fig.savefig(fig_dir / "worm_plot_ad.svg", bbox_inches="tight")

    if rocauc:
        fig, ax = plt.subplots(1, 1, figsize=(figsize, figsize))
        RocCurveDisplay.from_predictions(
            heatmap_df["Ad"].astype(int).ravel(),
            heatmap_df["Hit fraction"].ravel(),
            ax=ax,
            c="black",
        )
        plt.plot(
            np.linspace(0, 1, 100),
            np.linspace(0, 1, 100),
            c="r",
        )
        transform_plot_ax(ax, legend_title="", ftsize=ftsize, leg_ftsize=leg_ftsize)
        fig.savefig(fig_dir / "roc_auc_ad.svg", bbox_inches="tight")

        fig, ax = plt.subplots(1, 1, figsize=(figsize, figsize))
        PrecisionRecallDisplay.from_predictions(
            heatmap_df["Ad"].astype(int).ravel(),
            heatmap_df["Hit fraction"].ravel(),
            ax=ax,
            c="black",
        )
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        transform_plot_ax(ax, legend_title="", ftsize=ftsize, leg_ftsize=leg_ftsize)
        fig.savefig(fig_dir / "pr_auc_ad.svg", bbox_inches="tight")

        fig, ax = plt.subplots(1, 1, figsize=(figsize, figsize))
        ax = sns.boxplot(
            data=heatmap_df,
            x=hue_worm,
            y="Hit fraction",
            ax=ax,
            palette=hue_palette_worm,
            hue_order=[],
        )
        annot = Annotator(
            ax,
            pairs=[(0, 1)],
            data=heatmap_df,
            x="Ad",
            y="Hit fraction",
            order=[0, 1],
        )
        annot.configure(
            test="Mann-Whitney",
            loc="inside",
            show_test_name=False,
            verbose=2,
            comparisons_correction="BH",
            correction_format="replace",
        )
        annot.apply_test()
        ax, _ = annot.annotate()
        transform_plot_ax(ax, legend_title="", ftsize=ftsize, leg_ftsize=leg_ftsize)
        ax.set_xlabel("Adenoma (right)", fontsize=ftsize)
        ax.figure.savefig(fig_dir / "ad_dist_boxplot.svg", bbox_inches="tight")


def plot_polyp_size_nr_link(
    heatmap_df: pd.DataFrame,
    fig_dir: pl.Path,
    palette_size: Dict,
    palette_nr: Dict,
    ftsize: int = 10,
    leg_ftsize: int = 10,
) -> pd.DataFrame:
    fig, ax = plt.subplots(1, 1)
    sns.boxplot(
        data=heatmap_df,
        x="Polyp Size cat",
        y="Hit fraction",
        order=["None", "5mm", ">=6mm"],
        palette=palette_size,
        ax=ax,
    )
    annot = Annotator(
        ax,
        pairs=[("None", "5mm"), (">=6mm", "5mm"), (">=6mm", "None")],
        data=heatmap_df,
        x="Polyp Size cat",
        y="Hit fraction",
        order=["None", "5mm", ">=6mm"],
    )
    annot.configure(
        test="Mann-Whitney",
        loc="inside",
        show_test_name=False,
        verbose=2,
        comparisons_correction=None,
        fontsize=ftsize,
    )
    annot.apply_test()
    ax, test_results = annot.annotate()
    vc = pd.cut(
        heatmap_df["size_py_rght"],
        bins=[-1, 4, 5.5, 100],
        labels=["None", "5mm", ">=6mm"],
    ).value_counts()
    ax.set_xticklabels(
        [
            f"None\nN={vc.loc['None']}",
            f"5mm\nN={vc.loc['5mm']}",
            f">=6mm\nN={vc.loc['>=6mm']}",
        ]
    )
    transform_plot_ax(ax, legend_title="", ftsize=ftsize, leg_ftsize=leg_ftsize)
    fig.savefig(fig_dir / "polyp_size_hit_fraction.svg", bbox_inches="tight")

    fig, ax = plt.subplots(1, 1)
    sns.boxplot(
        data=heatmap_df,
        x="Polyp Nr Right",
        y="Hit fraction",
        order=["0", "1", ">=2"],
        palette=palette_nr,
        ax=ax,
    )
    annot = Annotator(
        ax,
        pairs=[
            ("0", "1"),
            ("0", ">=2"),
            ("1", ">=2"),
        ],
        data=heatmap_df,
        x="Polyp Nr Right",
        y="Hit fraction",
        order=["0", "1", ">=2"],
    )
    annot.configure(
        test="Mann-Whitney",
        loc="inside",
        show_test_name=False,
        verbose=2,
        comparisons_correction=None,
        fontsize=ftsize,
    )
    annot.apply_test()
    ax, test_results = annot.annotate()
    vc = heatmap_df["Polyp Nr Right"].value_counts()
    ax.set_xticklabels(
        [f"0\nN={vc.loc['0']}", f"1\nN={vc.loc['1']}", f">=2\nN={vc.loc['>=2']}"]
    )
    transform_plot_ax(ax, legend_title="", ftsize=ftsize, leg_ftsize=leg_ftsize)
    fig.savefig(fig_dir / "polyp_nr_right_hit_fraction.svg", bbox_inches="tight")

    fig, ax = plt.subplots(1, 1)
    sns.boxplot(
        data=heatmap_df,
        x="Polyp Nr Total",
        y="Hit fraction",
        order=["0", "1", ">=2"],
        palette=palette_nr,
        ax=ax,
    )
    annot = Annotator(
        ax,
        pairs=[
            ("0", "1"),
            ("0", ">=2"),
            ("1", ">=2"),
        ],
        data=heatmap_df,
        x="Polyp Nr Total",
        y="Hit fraction",
        order=["0", "1", ">=2"],
    )
    annot.configure(
        test="Mann-Whitney",
        loc="inside",
        show_test_name=False,
        verbose=2,
        comparisons_correction=None,
        fontsize=ftsize,
    )
    annot.apply_test()
    ax, test_results = annot.annotate()
    vc = heatmap_df["Polyp Nr Total"].value_counts()
    ax.set_xticklabels(
        [f"0\nN={vc.loc['0']}", f"1\nN={vc.loc['1']}", f">=2\nN={vc.loc['>=2']}"]
    )
    transform_plot_ax(ax, legend_title="", ftsize=ftsize, leg_ftsize=leg_ftsize)
    fig.savefig(fig_dir / "polyp_nr_total_hit_fraction.svg", bbox_inches="tight")
