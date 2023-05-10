import mygene
import pandas as pd
import numpy as np
import pathlib as pl
from typing import List, Union, Dict, Tuple


def tissue_type(x):
    if "AP" in x:
        return "Adenoma"
    elif "CR" in x:
        return "Healthy"
    elif "SSA/P" in x:
        return "SSL"
    elif "UR" in x:
        return "NAC"


def query_gene_name(ens: Union[pd.Index, np.ndarray]) -> Tuple[Dict, List]:
    mg = mygene.MyGeneInfo()
    ginfo = mg.querymany(ens, scopes="ensembl.gene")

    ens_mapping = {}
    not_found = []
    for resp in ginfo:
        if "symbol" not in resp:
            not_found.append(resp["query"])
            continue
        ens_mapping[resp["query"]] = resp["symbol"]
    return ens_mapping, not_found


def download_gex_data(
    path_right: Union[str, pl.Path], path_right_cr: Union[str, pl.Path]
) -> pd.DataFrame:
    right_data = pd.read_csv(path_right)
    right_data = pd.concat(
        [
            right_data[["Ensembl_ID"]],
            right_data.loc[:, right_data.columns.str.startswith("FPKM")],
        ],
        axis=1,
    )
    right_data = right_data.set_index("Ensembl_ID").T
    right_data.index = right_data.index.str.strip("FPKM ")

    right_data_cancer = pd.read_csv(path_right_cr)
    right_data_cancer = pd.concat(
        [
            right_data_cancer[["Ensembl_ID"]],
            right_data_cancer.loc[:, right_data_cancer.columns.str.startswith("FPKM")],
        ],
        axis=1,
    )
    right_data_cancer = right_data_cancer.set_index("Ensembl_ID").T
    right_data_cancer.index = right_data_cancer.index.str.strip("FPKM ")

    tissue = pd.DataFrame(
        right_data.index, index=right_data.index, columns=["Tissue type"]
    ).applymap(tissue_type)

    right_data["type"] = tissue
    right_data_cancer["type"] = ["Cancer"] * right_data_cancer.shape[0]

    ens = right_data.columns[:-1]

    ens_mapping, not_found = query_gene_name(ens)

    right_data = right_data.drop(not_found, axis=1)
    right_data = right_data.rename(columns=ens_mapping)
    right_data_cancer = right_data_cancer.drop(
        right_data_cancer.columns.intersection(not_found), axis=1
    )
    right_data_cancer = right_data_cancer.rename(columns=ens_mapping)
    right_data_cancer = right_data_cancer.loc[
        :, ~right_data_cancer.columns.str.startswith("ENSG")
    ]

    right_data = right_data.loc[:, ~right_data.columns.duplicated()]
    right_data_cancer = right_data_cancer.loc[
        :, ~right_data_cancer.columns.duplicated()
    ]

    right_data = pd.concat([right_data, right_data_cancer])
    # For some reason a left sample slipped in here
    right_data = right_data.drop("UL-1")
    right_data = right_data.dropna(axis=1)

    return right_data
