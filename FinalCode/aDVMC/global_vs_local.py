import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from .adVMP_discovery import get_heatmap_df
from .adVMP_plots import transform_plot_ax

MEANINGFUL_GROUPS = {
    "1_TssA": "Active promoter",
    "2_TssAFlnk": "Active promoter",
    "3_TxFlnk": "Active promoter",
    "4_Tx": "Transcription (body)",
    "5_TxWk": "Transcription (body)",
    "6_EnhG": "Enhancer",
    "7_Enh": "Enhancer",
    "8_ZNF/Rpts": "ZNF/Repeats",
    "10_TssBiv": "Bivalent promoter",
    "11_BivFlnk": "Bivalent promoter",
    "12_EnhBiv": "Bivalent enhancer",
    "13_ReprPC": "Repressed polycomb",
    "14_ReprPCWk": "Repressed polycomb",
    "15_Quies": "Quiescent",
}


def check_global_dysregulation_pergene(
    epic_manifest: pd.DataFrame,
    gene: str,
    EPIC2_b: pd.DataFrame,
    EPIC3_b: pd.DataFrame,
    EPIC4_b: pd.DataFrame,
    EPIC2_phenotypes: np.ndarray,
    EPIC3_phenotypes: np.ndarray,
    EPIC4_phenotypes: np.ndarray,
    binary1: pd.DataFrame,
    binary2: pd.DataFrame,
    binary3: pd.DataFrame,
    union_cpgs: np.ndarray,
) -> bool:
    background_cpgs = EPIC2_b.columns.intersection(EPIC4_b.columns)

    found_df = epic_manifest.loc[
        epic_manifest["UCSC_RefGene_Name"].str.contains(gene).fillna(False),
        ["UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR", "MAPINFO", "State"],
    ]

    found_df = found_df.sort_values(by="MAPINFO")

    found_df = clean_cg(found_df=found_df, gene=gene)

    found_df = found_df.loc[found_df.index.intersection(background_cpgs)]

    posmin, posmax = found_df["MAPINFO"].min(), found_df["MAPINFO"].max()
    found_df["relpos"] = (found_df["MAPINFO"] - posmin) / (posmax - posmin)

    signcpg = found_df.index.intersection(union_cpgs)[0]

    variable_pat1 = binary1[binary1.loc[:, signcpg] == 1].index
    variable_pat2 = binary2[binary2.loc[:, signcpg] == 1].index
    variable_pat3 = binary3[binary3.loc[:, signcpg] == 1].index

    dys1, dys2, dys3 = get_global_local_status(
        found_df=found_df,
        signcpg=signcpg,
        EPIC2_b=EPIC2_b,
        EPIC3_b=EPIC3_b,
        EPIC4_b=EPIC4_b,
        EPIC2_phenotypes=EPIC2_phenotypes,
        EPIC3_phenotypes=EPIC3_phenotypes,
        EPIC4_phenotypes=EPIC4_phenotypes,
        variable_pat1=variable_pat1,
        variable_pat2=variable_pat2,
        variable_pat3=variable_pat3,
    )

    isglobal = global_indicator(dys1, dys2, dys3)

    return isglobal


from typing import Tuple


def get_dysregulation_perc(
    flankingcpgs: pd.Index,
    EPIC_b: pd.DataFrame,
    EPIC_phenotypes: pd.DataFrame,
    variable_pat: pd.Index,
    hit_limit: int = 4,
) -> float:
    heatmap_flanking, _ = get_heatmap_df(
        selcpgs=flankingcpgs, EPIC_m=EPIC_b, phenotypes=EPIC_phenotypes, bal=True
    )

    n_dysregulated = (heatmap_flanking.loc[:, flankingcpgs].abs() > hit_limit).sum(
        axis=1
    )

    perc_mult_dysregulated = (n_dysregulated.loc[variable_pat] > 1).sum() / len(
        variable_pat
    )

    return perc_mult_dysregulated


def get_global_local_status(
    found_df: pd.DataFrame,
    signcpg: str,
    EPIC2_b: pd.DataFrame,
    EPIC3_b: pd.DataFrame,
    EPIC4_b: pd.DataFrame,
    EPIC2_phenotypes: np.ndarray,
    EPIC3_phenotypes: np.ndarray,
    EPIC4_phenotypes: np.ndarray,
    variable_pat1: pd.Index,
    variable_pat2: pd.Index,
    variable_pat3: pd.Index,
    flanking_region: int = 2500,
) -> Tuple[float, float, float]:
    pos = found_df.loc[signcpg].MAPINFO
    flankingpos = pos - flanking_region, pos + flanking_region

    flankingcpgs = found_df[
        (found_df["MAPINFO"] >= flankingpos[0])
        & (found_df["MAPINFO"] <= flankingpos[1])
    ].index

    if len(flankingcpgs) < 2:
        print("There is only 1 other CpG in this region...")
        return 0, 0, 0

    dys1 = get_dysregulation_perc(
        flankingcpgs=flankingcpgs,
        EPIC_b=EPIC2_b,
        EPIC_phenotypes=EPIC2_phenotypes,
        variable_pat=variable_pat1,
    )
    dys2 = get_dysregulation_perc(
        flankingcpgs=flankingcpgs,
        EPIC_b=EPIC3_b,
        EPIC_phenotypes=EPIC3_phenotypes,
        variable_pat=variable_pat2,
    )
    dys3 = get_dysregulation_perc(
        flankingcpgs=flankingcpgs,
        EPIC_b=EPIC4_b,
        EPIC_phenotypes=EPIC4_phenotypes,
        variable_pat=variable_pat3,
    )

    return dys1, dys2, dys3


def global_indicator(dys1: float, dys2: float, dys3: float) -> bool:
    if (
        ((dys1 >= 0.5) and (dys2 >= 0.5))
        or ((dys1 >= 0.5) and (dys3 >= 0.5))
        or ((dys2 >= 0.5) and (dys3 >= 0.5))
    ):
        return True
    else:
        return False


def get_lineplot(
    EPIC_b: pd.DataFrame,
    found_df: pd.DataFrame,
    variable_pat: np.ndarray,
    nvariable_pat: np.ndarray,
    promoter_reg: pd.DataFrame,
    background_cpgs: pd.Index,
    title: str,
    promoter_only: bool,
) -> plt.Axes:
    sel_cgs_desc_var = (
        EPIC_b.loc[variable_pat, found_df.index.intersection(background_cpgs)]
        .describe()
        .T
    )
    sel_cgs_desc_nvar = (
        EPIC_b.loc[nvariable_pat, found_df.index.intersection(background_cpgs)]
        .describe()
        .T
    )

    promoter_reg = promoter_reg.loc[found_df.index.intersection(background_cpgs)]

    sel_cgs_desc_var = pd.concat(
        [sel_cgs_desc_var, found_df["relpos"]], axis=1, join="inner"
    )
    sel_cgs_desc_nvar = pd.concat(
        [sel_cgs_desc_nvar, found_df["relpos"]], axis=1, join="inner"
    )

    sel_cgs_desc_var = sel_cgs_desc_var.reset_index().reset_index()
    sel_cgs_desc_nvar = sel_cgs_desc_nvar.reset_index().reset_index()

    fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    sns.lineplot(
        data=sel_cgs_desc_var,
        x="relpos",
        y="50%",
        ci=None,
        color="r",
        label=f"Var. N={len(variable_pat)}",
    )
    sns.scatterplot(data=sel_cgs_desc_var, x="relpos", y="50%", color="r")
    ax.fill_between(
        sel_cgs_desc_var.relpos,
        sel_cgs_desc_var["25%"],
        sel_cgs_desc_var["75%"],
        alpha=0.2,
        color="r",
    )
    sns.lineplot(
        data=sel_cgs_desc_nvar,
        x="relpos",
        y="50%",
        ci=None,
        color="b",
        label=f"N. Var. N={len(nvariable_pat)}",
    )
    sns.scatterplot(data=sel_cgs_desc_nvar, x="relpos", y="50%", color="b")
    ax.fill_between(
        sel_cgs_desc_nvar.relpos,
        sel_cgs_desc_nvar["25%"],
        sel_cgs_desc_nvar["75%"],
        alpha=0.2,
        color="b",
    )
    transform_plot_ax(ax, legend_title="Variable patients", leg_ftsize=10)
    ax.set_xlabel("")
    ax.set_title(title)

    if not promoter_only:
        for x, cg in enumerate(promoter_reg.index):
            if promoter_reg.loc[cg] == "Active promoter":
                ax.text(
                    sel_cgs_desc_var.iloc[x]["relpos"], 1, "__", color="orange", size=15
                )
            elif promoter_reg.loc[cg] == "Bivalent promoter":
                ax.text(
                    sel_cgs_desc_var.iloc[x]["relpos"], 1, "__", color="green", size=15
                )

    return ax


def clean_cg(found_df: pd.DataFrame, gene: str) -> pd.DataFrame:
    cg_to_keep = []
    listvals = found_df["UCSC_RefGene_Name"].str.split(";")
    for cg in listvals.index:
        if gene in listvals.loc[cg]:
            cg_to_keep.append(cg)
    return found_df.loc[cg_to_keep]


def get_full_cg_info_gene(
    epic_manifest: pd.DataFrame,
    gene: str,
    union_cpgs: np.ndarray,
    binary1: pd.DataFrame,
    binary2: pd.DataFrame,
    binary3: pd.DataFrame,
    EPIC2_b: pd.DataFrame,
    EPIC3_b: pd.DataFrame,
    EPIC4_b: pd.DataFrame,
    promoter_only: bool = True,
) -> Tuple[plt.Axes, plt.Axes, plt.Axes]:
    background_cpgs = EPIC2_b.columns.intersection(EPIC4_b.columns)

    found_df = epic_manifest.loc[
        epic_manifest["UCSC_RefGene_Name"].str.contains(gene).fillna(False),
        ["UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR", "MAPINFO", "State"],
    ]

    found_df = found_df.sort_values(by="MAPINFO")

    found_df = clean_cg(found_df=found_df, gene=gene)

    promoter_reg = found_df["State"].replace(MEANINGFUL_GROUPS)

    posmin, posmax = found_df["MAPINFO"].min(), found_df["MAPINFO"].max()
    found_df["relpos"] = (found_df["MAPINFO"] - posmin) / (posmax - posmin)

    if promoter_only:
        found_df = found_df[promoter_reg.isin(["Active promoter", "Bivalent promoter"])]
        promoter_reg = promoter_reg[
            promoter_reg.isin(["Active promoter", "Bivalent promoter"])
        ]

    print(found_df.index.intersection(union_cpgs))
    signcpg = found_df.index.intersection(union_cpgs)[0]

    nvariable_pat1 = binary1[binary1.loc[:, signcpg] == 0].index
    variable_pat1 = binary1[binary1.loc[:, signcpg] == 1].index
    nvariable_pat2 = binary2[binary2.loc[:, signcpg] == 0].index
    variable_pat2 = binary2[binary2.loc[:, signcpg] == 1].index
    nvariable_pat3 = binary3[binary3.loc[:, signcpg] == 0].index
    variable_pat3 = binary3[binary3.loc[:, signcpg] == 1].index

    print(signcpg, len(variable_pat1), len(nvariable_pat1))
    ax1 = get_lineplot(
        EPIC_b=EPIC2_b,
        found_df=found_df,
        variable_pat=variable_pat1,
        nvariable_pat=nvariable_pat1,
        promoter_reg=promoter_reg,
        title=f"SWEPIC1 {gene}",
        promoter_only=promoter_only,
        background_cpgs=background_cpgs,
    )

    ax2 = get_lineplot(
        EPIC_b=EPIC3_b,
        found_df=found_df,
        variable_pat=variable_pat2,
        nvariable_pat=nvariable_pat2,
        promoter_reg=promoter_reg,
        title=f"SWEPIC2 {gene}",
        promoter_only=promoter_only,
        background_cpgs=background_cpgs,
    )

    ax3 = get_lineplot(
        EPIC_b=EPIC4_b,
        found_df=found_df,
        variable_pat=variable_pat3,
        nvariable_pat=nvariable_pat3,
        promoter_reg=promoter_reg,
        title=f"SWEPIC3 {gene}",
        promoter_only=promoter_only,
        background_cpgs=background_cpgs,
    )
    return ax1, ax2, ax3
