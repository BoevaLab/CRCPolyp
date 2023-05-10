import os
import numpy as np
import pathlib as pl
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from typing import Dict, Tuple, Optional
from tqdm import tqdm
from collections import defaultdict
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import average_precision_score, roc_auc_score
from scipy.stats import median_abs_deviation

from .adVMP_discovery import get_hyper_vDMC
from utils.plotting import transform_plot_ax


def get_stratified_hyper_DMC(
    y: np.ndarray,
    EPIC_m: pd.DataFrame,
    result_dir: pl.Path,
    n_splits: int = 5,
    rs: int = 0,
) -> None:
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=rs)
    os.makedirs(result_dir, exist_ok=True)

    for i, (train_index, test_index) in enumerate(skf.split(EPIC_m, y)):
        print(f"Starting fold {i}...")
        os.makedirs(result_dir / f"fold{i}", exist_ok=True)

        EPIC_train, EPIC_test = EPIC_m.iloc[train_index], EPIC_m.iloc[test_index]
        EPIC_y_train, _ = y[train_index], y[test_index]

        test_results = get_hyper_vDMC(methylation=EPIC_train, phenotypes=EPIC_y_train)

        test_results.to_csv(result_dir / f"fold{i}" / "adVMP_right.csv")

        pd.Series(EPIC_train.index).to_csv(
            result_dir / f"fold{i}" / "sample_train_set.csv"
        )
        pd.Series(EPIC_test.index).to_csv(
            result_dir / f"fold{i}" / "sample_test_set.csv"
        )


def get_fold_specific_ensembling_cpgs(
    test_results: Dict[str, pd.DataFrame],
    noens: bool = False,
) -> Dict[str, np.ndarray]:
    sel_probes = {}
    for ds in test_results:
        sel_probes[ds] = defaultdict(list)
        for i, fold in enumerate(test_results[ds]):
            df = test_results[ds][fold]
            sign = df[(df["q"] < 0.001) & (df["ttest_p"] < 0.05) & (df.diffV > 0)]
            sel_probes[ds][fold].append(sign.index)
    if noens:
        sel_probes = {ds: dict(sel_probes[ds]) for ds in sel_probes}
        return sel_probes

    datasets = sorted(list(test_results.keys()))

    fold_probes = {
        ds: {fold: defaultdict(list) for fold in sel_probes[ds]} for ds in sel_probes
    }
    for ds1 in datasets:
        for fold1 in fold_probes[ds1]:
            fold_spec_probes = sel_probes[ds1][fold1]
            intersections = []
            for ds2 in datasets:
                if ds1 == ds2:
                    continue
                else:
                    for fold2 in sel_probes[ds2]:
                        intersections.append(
                            np.intersect1d(fold_spec_probes, sel_probes[ds2][fold2])
                        )

            intersections = np.unique(np.concatenate(intersections))
            fold_probes[ds1][fold1] = intersections
    return fold_probes


def get_hit_fraction(
    EPIC_m: pd.DataFrame,
    selcpgs: np.ndarray,
    median: Optional[float] = None,
    mad: Optional[float] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if median is None:
        median = EPIC_m[EPIC_m.columns.intersection(selcpgs)].median()
    if mad is None:
        mad = median_abs_deviation(EPIC_m[EPIC_m.columns.intersection(selcpgs)])
        mad = pd.Series(mad, index=EPIC_m.columns.intersection(selcpgs))

    EPIC_z_score = EPIC_m[EPIC_m.columns.intersection(selcpgs)]
    EPIC_z_score = (EPIC_z_score - median) / mad

    hit_EPIC = (EPIC_z_score.abs() > 4).astype(int)
    hit_fraction = hit_EPIC.sum(axis=1) / hit_EPIC.shape[1]
    hit_fraction.name = "Hit fraction"
    return EPIC_z_score, hit_fraction, median, mad


def get_fold_results(
    EPIC_m: pd.DataFrame,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    selcpgs: np.ndarray,
    y: pd.DataFrame,
    estimate_copa: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    EPIC_train, EPIC_test = EPIC_m.iloc[train_idx], EPIC_m.iloc[test_idx]
    EPIC_y_train, EPIC_y_test = y[train_idx], y[test_idx]

    EPIC_z_score_train, hit_fraction_train, median, mad = get_hit_fraction(
        EPIC_m=EPIC_train, selcpgs=selcpgs
    )
    if estimate_copa:
        EPIC_z_score_test, hit_fraction_test, _, _ = get_hit_fraction(
            EPIC_m=EPIC_test, selcpgs=selcpgs, median=median, mad=mad
        )
    else:
        EPIC_z_score_test, hit_fraction_test, _, _ = get_hit_fraction(
            EPIC_m=EPIC_test, selcpgs=selcpgs
        )

    heatmap_df_train = pd.concat(
        [
            EPIC_z_score_train,
            pd.DataFrame(EPIC_y_train, index=EPIC_train.index, columns=["Ad"]),
            hit_fraction_train,
        ],
        axis=1,
    )

    heatmap_df_test = pd.concat(
        [
            EPIC_z_score_test,
            pd.DataFrame(EPIC_y_test, index=EPIC_test.index, columns=["Ad"]),
            hit_fraction_test,
        ],
        axis=1,
    )

    roctrain = roc_auc_score(
        y_true=heatmap_df_train["Ad"].astype(int).ravel(),
        y_score=heatmap_df_train["Hit fraction"].ravel(),
    )
    roctest = roc_auc_score(
        y_true=heatmap_df_test["Ad"].astype(int).ravel(),
        y_score=heatmap_df_test["Hit fraction"].ravel(),
    )
    auprctrain = average_precision_score(
        y_true=heatmap_df_train["Ad"].astype(int).ravel(),
        y_score=heatmap_df_train["Hit fraction"].ravel(),
    )
    auprctest = average_precision_score(
        y_true=heatmap_df_test["Ad"].astype(int).ravel(),
        y_score=heatmap_df_test["Hit fraction"].ravel(),
    )

    return (
        heatmap_df_train,
        heatmap_df_test,
        np.array([roctrain, roctest, auprctrain, auprctest]),
    )


def get_crossval_performance(
    ds_dir: pl.Path,
    EPIC_b: pd.DataFrame,
    union_cpgs_fold_spec: Dict,
    EPIC_phenotypes: np.ndarray,
    estimate_copa: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    all_stats = []
    crossval_hit_fraction = []
    for fold in ds_dir.iterdir():
        if fold.stem == ".DS_Store":
            continue
        test_samples = pd.read_csv(fold / "sample_test_set.csv", index_col=0).astype(
            str
        )
        train_samples = pd.read_csv(fold / "sample_train_set.csv", index_col=0).astype(
            str
        )

        _, idx_train, _ = np.intersect1d(
            EPIC_b.index, train_samples.values.ravel(), return_indices=True
        )
        _, idx_test, _ = np.intersect1d(
            EPIC_b.index, test_samples.values.ravel(), return_indices=True
        )

        _, heatmap_df_test, stats = get_fold_results(
            EPIC_m=EPIC_b,
            train_idx=idx_train,
            test_idx=idx_test,
            selcpgs=union_cpgs_fold_spec[fold.stem],
            y=EPIC_phenotypes.astype(int),
            estimate_copa=estimate_copa,
        )
        all_stats.append(stats)

        crossval_hit_fraction.append(heatmap_df_test[["Ad", "Hit fraction"]])

    all_stats = pd.DataFrame(
        all_stats,
        index=[f"fold{i+1}" for i in range(len(all_stats))],
        columns=["ROC AUC train", "ROC AUC test", "PR AUC train", "PR AUC test"],
    ).T

    all_stats["Mean"] = all_stats.mean(axis=1)
    all_stats["Std"] = all_stats.std(axis=1)

    crossval_hit_fraction = pd.concat(crossval_hit_fraction)
    crossval_hit_fraction["Ad_plot"] = crossval_hit_fraction["Ad"].replace(
        {0: "No", 1: "Yes"}
    )
    crossval_hit_fraction = crossval_hit_fraction.sort_values(by="Hit fraction")
    crossval_hit_fraction["Order"] = np.arange(crossval_hit_fraction.shape[0])

    return all_stats, crossval_hit_fraction


def get_rdn_cpgs_oneiter(
    union_cpgs_fold_spec: Dict,
    data: pd.DataFrame,
    background_cpgs: Optional[np.ndarray] = None,
) -> Dict:
    rdn_cpgs_pf = {}
    for fold in union_cpgs_fold_spec:
        if background_cpgs is None:
            rdn_cpgs = np.random.choice(
                data.columns,
                size=(len(data.columns.intersection(union_cpgs_fold_spec[fold]))),
                replace=False,
            )
        else:
            rdn_cpgs = np.random.choice(
                data.columns.intersection(background_cpgs),
                size=(len(data.columns.intersection(union_cpgs_fold_spec[fold]))),
                replace=False,
            )
        rdn_cpgs_pf[fold] = rdn_cpgs
    return rdn_cpgs_pf


def get_rdn_cpgs(
    union_cpgs_fold_spec: Dict,
    data: pd.DataFrame,
    background_cpgs: Optional[np.ndarray] = None,
    n_iter: int = 500,
) -> Dict:
    rdn_cpgs_iters = {}
    for i in tqdm(range(n_iter)):
        rdn_cpgs_iters[i] = get_rdn_cpgs_oneiter(
            union_cpgs_fold_spec=union_cpgs_fold_spec,
            data=data,
            background_cpgs=background_cpgs,
        )
    return rdn_cpgs_iters


def get_crossval_performance_rdn(
    ds_dir: pl.Path,
    EPIC_b: pd.DataFrame,
    union_cpgs_fold_spec: Dict,
    EPIC_phenotypes: np.ndarray,
    estimate_copa: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    all_stats = []
    for fold in ds_dir.iterdir():
        if fold.stem == ".DS_Store":
            continue
        test_samples = pd.read_csv(fold / "sample_test_set.csv", index_col=0).astype(
            str
        )
        train_samples = pd.read_csv(fold / "sample_train_set.csv", index_col=0).astype(
            str
        )

        _, idx_train, _ = np.intersect1d(
            EPIC_b.index, train_samples.values.ravel(), return_indices=True
        )
        _, idx_test, _ = np.intersect1d(
            EPIC_b.index, test_samples.values.ravel(), return_indices=True
        )

        _, _, stats = get_fold_results(
            EPIC_m=EPIC_b,
            train_idx=idx_train,
            test_idx=idx_test,
            selcpgs=union_cpgs_fold_spec[fold.stem],
            y=EPIC_phenotypes.astype(int),
            estimate_copa=estimate_copa,
        )
        all_stats.append(stats)

    all_stats = pd.DataFrame(
        all_stats,
        index=[f"fold{i+1}" for i in range(len(all_stats))],
        columns=["ROC AUC train", "ROC AUC test", "PR AUC train", "PR AUC test"],
    ).T

    all_stats["Mean"] = all_stats.mean(axis=1)

    return all_stats.loc["ROC AUC test", "Mean"]


def format_p(p, size):
    if p == 0:
        return f"p<{1/size}"
    else:
        return f"p={p}"


def get_comparison_rdn(
    background_cpgs: np.ndarray,
    figdir: pl.Path,
    ref: float,
    ds_dir: pl.Path,
    phenotypes: np.ndarray,
    union_cpgs_fold_spec: Dict,
    data: pd.DataFrame,
    clin: pd.DataFrame,
    n_iter: int = 100,
) -> None:
    rdn_cpgs = get_rdn_cpgs(
        union_cpgs_fold_spec=union_cpgs_fold_spec, data=data, n_iter=n_iter
    )
    rdn_cpgs_bck = get_rdn_cpgs(
        union_cpgs_fold_spec=union_cpgs_fold_spec,
        data=data,
        n_iter=n_iter,
        background_cpgs=background_cpgs,
    )
    avg_perf = []
    for i in rdn_cpgs:
        score = get_crossval_performance_rdn(
            ds_dir=ds_dir,
            EPIC_b=data,
            union_cpgs_fold_spec=rdn_cpgs[i],
            EPIC_phenotypes=phenotypes,
            estimate_copa=True,
        )
        avg_perf.append(score)

    avg_perf_back = []
    for i in rdn_cpgs_bck:
        score = get_crossval_performance_rdn(
            ds_dir=ds_dir,
            EPIC_b=data,
            union_cpgs_fold_spec=rdn_cpgs_bck[i],
            EPIC_phenotypes=phenotypes,
            estimate_copa=True,
        )
        avg_perf_back.append(score)

    pval_rdn = np.sum(ref < np.array(avg_perf)) / len(avg_perf)
    pval_rdn_bck = np.sum(ref < np.array(avg_perf_back)) / len(avg_perf_back)
    df = pd.concat(
        [
            pd.DataFrame(avg_perf, columns=["Random CpGs"]),
            pd.DataFrame(avg_perf_back, columns=["Random Highly Variable CpGs"]),
        ],
        axis=1,
    )

    age_ref = roc_auc_score(y_score=clin["Age at visit"].ravel(), y_true=phenotypes)
    fig, ax = plt.subplots(1, 1)
    sns.boxplot(data=df, ax=ax, palette="bright")
    ax.set_ylim([0.4, 1])
    xmin, xmax = ax.get_xlim()
    # vertical lines
    ax.hlines(y=0.5, xmin=xmin, xmax=xmax, colors="r", lw=2, label="Random")
    ax.hlines(
        y=age_ref, xmin=xmin, xmax=xmax, colors="orange", ls="--", lw=2, label="Age"
    )
    ax.hlines(y=ref, xmin=xmin, xmax=xmax, colors="b", ls="--", lw=2, label="adVMP")

    ax.text(-0.1, 1.1 * ref, format_p(pval_rdn, len(avg_perf)), fontdict={"size": 15})
    ax.text(
        0.9,
        1.1 * ref,
        format_p(pval_rdn_bck, len(avg_perf_back)),
        fontdict={"size": 15},
    )

    transform_plot_ax(ax, legend_title="")
    fig.savefig(figdir / "comparison_rdn_cpgs.svg", bbox_inches="tight")
