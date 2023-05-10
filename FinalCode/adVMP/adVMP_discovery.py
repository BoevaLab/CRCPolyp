import pandas as pd
import numpy as np

from typing import Tuple, Optional
from scipy.stats import bartlett, levene, ttest_ind, median_abs_deviation
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import roc_auc_score
from tqdm import tqdm


def get_hyper_vDMC(
    methylation: pd.DataFrame, phenotypes: np.ndarray, stat_test=bartlett
) -> pd.DataFrame:
    bartlett_list = []
    ttest_list = []
    hyper_list = []
    for cg in tqdm(methylation.columns):
        df = methylation.loc[:, cg]
        bartlett_list.append(stat_test(df[phenotypes == 0], df[phenotypes == 1])[1])
        ttest_list.append(ttest_ind(df[phenotypes == 0], df[phenotypes == 1])[1])
        hyper = df[phenotypes == 1].std() - df[phenotypes == 0].std()
        hyper_list.append(hyper)

    qvalues = multipletests(bartlett_list, method="fdr_bh")[1]

    test_results = pd.DataFrame(
        np.array([qvalues, np.array(hyper_list), np.array(ttest_list)]),
        index=["q", "diffV", "ttest_p"],
        columns=methylation.columns,
    ).T
    return test_results


def get_significance(
    ref: float,
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
        if rdn_scores[-1] > ref:
            higher.append(1)
        else:
            higher.append(0)
    return rdn_scores, higher


def estimate_median_mad(
    EPIC_m: pd.DataFrame,
    phenotypes: np.ndarray,
    selcpgs: np.ndarray,
    make_balanced: bool = False,
) -> Tuple[pd.Series, pd.Series]:
    ad_cases = np.where(phenotypes == 1)[0]
    nad_cases = np.where(phenotypes == 0)[0]

    if make_balanced:
        n_ad = len(ad_cases)
        pats = ad_cases
        rdn_sel = np.random.choice(nad_cases, size=n_ad, replace=False)
        pats = np.concatenate([pats, rdn_sel])
        df = EPIC_m.iloc[pats]
    else:
        df = EPIC_m

    median = df[df.columns.intersection(selcpgs)].median()
    mad = median_abs_deviation(df[df.columns.intersection(selcpgs)])
    mad = pd.Series(mad, index=df.columns.intersection(selcpgs))

    return median, mad


def get_heatmap_df(
    selcpgs: np.ndarray,
    EPIC_m: pd.DataFrame,
    phenotypes: np.ndarray,
    bal: bool = False,
    lim_hit: int = 4,
) -> Tuple[pd.DataFrame, pd.Series]:
    median, mad = estimate_median_mad(
        EPIC_m=EPIC_m, phenotypes=phenotypes, selcpgs=selcpgs, make_balanced=bal
    )

    EPIC_z_score = EPIC_m[EPIC_m.columns.intersection(selcpgs)]
    EPIC_z_score = (EPIC_z_score - median) / mad

    hit_EPIC = (EPIC_z_score.abs() > lim_hit).astype(int)
    hit_fraction = hit_EPIC.sum(axis=1) / hit_EPIC.shape[1]
    hit_fraction.name = "Hit fraction"

    heatmap_df = pd.concat(
        [
            EPIC_z_score,
            pd.DataFrame(phenotypes.astype(int), index=EPIC_m.index, columns=["Ad"]),
            hit_fraction,
        ],
        axis=1,
    )
    heatmap_df["Ad_plot"] = heatmap_df["Ad"].replace({0: "No", 1: "Yes"})
    heatmap_df = heatmap_df.sort_values(by="Hit fraction")
    heatmap_df["Order"] = np.arange(heatmap_df.shape[0])

    return heatmap_df, hit_fraction


def get_polyp_size_nr_link(
    EPIC_clin: pd.DataFrame, heatmap_df: pd.DataFrame
) -> pd.DataFrame:
    EPIC_clin["polyps_total_nr"] = EPIC_clin["polyps_total_nr"].replace({" ": 0})
    EPIC_clin["polyps_total_size"] = EPIC_clin["polyps_total_size"].replace({" ": 0})

    polyp_info = EPIC_clin[
        ["polyps_total_nr", "polyps_total_size", "polyps_right_nr", "size_py_rght"]
    ].astype(float)

    heatmap_df = heatmap_df.copy()
    heatmap_df = pd.concat([heatmap_df, polyp_info], axis=1)
    heatmap_df["Polyp Nr Right"] = heatmap_df["polyps_right_nr"].apply(
        lambda x: str(int(x)) if x < 2 else ">=2"
    )
    heatmap_df["Polyp Nr Total"] = heatmap_df["polyps_total_nr"].apply(
        lambda x: str(int(x)) if x < 2 else ">=2"
    )

    heatmap_df["Polyp Size cat"] = pd.cut(
        heatmap_df["size_py_rght"],
        bins=[-1, 4, 5.5, 100],
        labels=["None", "5mm", ">=6mm"],
    )

    heatmap_df["Polyp Size Total cat"] = pd.cut(
        heatmap_df["polyps_total_size"],
        bins=[-1, 4, 5.5, 100],
        labels=["None", "5mm", ">=6mm"],
    )

    return heatmap_df
