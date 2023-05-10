import pandas as pd
import numpy as np
import pathlib as pl
import seaborn as sns
import matplotlib.pyplot as plt

from tqdm import tqdm
from scipy.stats import median_abs_deviation
from sklearn.metrics import roc_auc_score
from typing import Tuple, Optional, Dict

from utils.plotting import transform_plot_ax
from .adVMP_crossval import format_p


def get_significance(
    data: pd.DataFrame,
    clin: pd.DataFrame,
    union_cpgs: np.ndarray,
    background_cpgs: Optional[np.ndarray] = None,
    nadj: bool = True,
    hit_limit: float = 2,
    risk_col: str = "Ad_risk",
    exclude_one: bool = True,
    n_iter: int = 500,
) -> Tuple[np.ndarray, np.ndarray]:
    rdn_scores = []
    higher = []
    for i in tqdm(range(n_iter)):
        if background_cpgs is None:
            rdn_cpgs = np.random.choice(
                data.columns,
                size=(len(data.columns.intersection(union_cpgs))),
                replace=False,
            )
        else:
            rdn_cpgs = np.random.choice(
                data.columns.intersection(background_cpgs),
                size=(len(data.columns.intersection(union_cpgs))),
                replace=False,
            )
        selcpgs = rdn_cpgs

        median_test = data.loc[:, data.columns.intersection(selcpgs)].median()
        mad_test = median_abs_deviation(data.loc[:, data.columns.intersection(selcpgs)])
        mad_test = pd.Series(mad_test, index=data.columns.intersection(selcpgs))

        test_z_score = data[data.columns.intersection(selcpgs)]
        test_z_score = (test_z_score - median_test) / mad_test

        hit_test = (test_z_score.abs() > hit_limit).astype(int)
        hit_fraction = hit_test.sum(axis=1) / hit_test.shape[1]
        hit_fraction.name = "Hit fraction"

        heatmap_df = pd.concat([test_z_score, clin[risk_col], hit_fraction], axis=1)

        if nadj:
            red_df = heatmap_df[heatmap_df[risk_col].isin([0, 1])]
        else:
            red_df = heatmap_df.copy()
            if exclude_one:
                red_df = red_df[red_df[risk_col] != 1]
            else:
                red_df[risk_col] = red_df[risk_col].replace({i: 0 for i in range(2)})
            red_df[risk_col] = red_df[risk_col].replace(
                {i: 1 for i in range(2, int(red_df[risk_col].max()) + 1)}
            )

        rdn_scores.append(
            roc_auc_score(
                y_true=red_df[risk_col].ravel(), y_score=red_df["Hit fraction"].ravel()
            )
        )

    return rdn_scores


def get_comparison_rdn_val(
    background_cpgs: np.ndarray,
    figdir: pl.Path,
    ref: float,
    phenotypes: np.ndarray,
    union_cpgs: Dict,
    data: pd.DataFrame,
    clin: pd.DataFrame,
    n_iter: int = 100,
    nadj: bool = True,
    hit_limit: int = 4,
    risk_col: str = "Ad_risk",
    exclude_one: bool = True,
    age_col: Optional[str] = "Age at visit",
) -> None:
    avg_perf = get_significance(
        data=data,
        clin=clin,
        union_cpgs=union_cpgs,
        background_cpgs=None,
        nadj=nadj,
        hit_limit=hit_limit,
        risk_col=risk_col,
        exclude_one=exclude_one,
        n_iter=n_iter,
    )
    avg_perf_back = get_significance(
        data=data,
        clin=clin,
        union_cpgs=union_cpgs,
        background_cpgs=background_cpgs,
        nadj=nadj,
        hit_limit=hit_limit,
        risk_col=risk_col,
        exclude_one=exclude_one,
        n_iter=n_iter,
    )

    pval_rdn = np.sum(ref < np.array(avg_perf)) / len(avg_perf)
    pval_rdn_bck = np.sum(ref < np.array(avg_perf_back)) / len(avg_perf_back)
    df = pd.concat(
        [
            pd.DataFrame(avg_perf, columns=["Random CpGs"]),
            pd.DataFrame(avg_perf_back, columns=["Random Highly Variable CpGs"]),
        ],
        axis=1,
    )

    if age_col is not None:
        if nadj:
            patloc = np.where((phenotypes == 0) | (phenotypes == 1))[0]
            cdf = clin.iloc[patloc]
            age_ref = roc_auc_score(
                y_score=cdf[age_col].ravel(), y_true=phenotypes[patloc]
            )
        else:
            if exclude_one:
                tokeep = np.where(phenotypes != 1)[0]
                cdf = clin.iloc[tokeep].copy()
                redphen = phenotypes[tokeep]
                redphen = pd.Series(redphen).replace(
                    {i: 1 for i in range(2, int(np.max(phenotypes)) + 1)}
                )
                age_ref = roc_auc_score(y_score=cdf[age_col].ravel(), y_true=redphen)
            else:
                redphen = phenotypes
                redphen = pd.Series(redphen).replace({i: 0 for i in range(2)}).ravel()
                redphen = (
                    pd.Series(redphen)
                    .replace({i: 1 for i in range(2, int(np.max(phenotypes)) + 1)})
                    .ravel()
                )
                age_ref = roc_auc_score(y_score=clin[age_col].ravel(), y_true=redphen)

    fig, ax = plt.subplots(1, 1)
    sns.boxplot(data=df, ax=ax, palette="bright")
    ax.set_ylim([0.4, 1])
    xmin, xmax = ax.get_xlim()
    # vertical lines
    ax.hlines(y=0.5, xmin=xmin, xmax=xmax, colors="r", lw=2, label="Random")
    if age_col is not None:
        ax.hlines(
            y=age_ref, xmin=xmin, xmax=xmax, colors="orange", ls="--", lw=2, label="Age"
        )
    ax.hlines(y=ref, xmin=xmin, xmax=xmax, colors="b", ls="--", lw=2, label="adVMP")

    ax.text(-0.1, 1.1 * ref, format_p(pval_rdn, len(avg_perf)), fontdict={"size": 15})
    ax.text(
        0.9,
        min(1.1 * ref, 0.8),
        format_p(pval_rdn_bck, len(avg_perf_back)),
        fontdict={"size": 15},
    )

    transform_plot_ax(ax, legend_title="")
    fig.savefig(figdir / "comparison_rdn_cpgs.svg", bbox_inches="tight")
