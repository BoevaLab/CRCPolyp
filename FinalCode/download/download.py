import pandas as pd
import numpy as np
import pathlib as pl
from tqdm import tqdm

from typing import Union, List, Tuple


def get_EPIC_data(
    sample_origin: pd.DataFrame,
    EPIC_name: str,
    specimens_to_exclude_right: pd.Index,
    b_values: pd.DataFrame,
    ad_right: pd.Series,
    clinical: pd.DataFrame,
    bad_probes: np.ndarray,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Index, np.ndarray]:
    EPIC_samples = sample_origin[
        (sample_origin.batch == EPIC_name) & (sample_origin.specimen_type == "cecum")
    ]
    EPIC_samples = EPIC_samples.drop(
        EPIC_samples.index.intersection(specimens_to_exclude_right)
    ).index

    EPIC_b = b_values.loc[EPIC_samples]

    EPIC_phenotypes = ad_right.loc[EPIC_samples].ravel()

    EPIC_b = EPIC_b.loc[~EPIC_phenotypes.isna()]
    EPIC_phenotypes = EPIC_phenotypes[~EPIC_phenotypes.isna()]

    EPIC_clin = clinical.set_index("Patient ID").loc[
        sample_origin.loc[EPIC_b.index].patient_id.ravel()
    ]
    EPIC_clin.index = EPIC_b.index

    # removing potential sketchy probes
    EPIC_b = EPIC_b.drop(EPIC_b.columns.intersection(bad_probes), axis=1)

    return EPIC_b, EPIC_clin, EPIC_samples, EPIC_phenotypes


def download_EPIC(
    sample_origin_path: Union[pl.Path, str],
    base_dir: pl.Path,
    clinical_path: Union[pl.Path, str],
    target_path: Union[pl.Path, str],
    bad_probes: np.ndarray,
    EPIC4: bool,
) -> Tuple:
    sample_origin = pd.read_csv(sample_origin_path, index_col=0)
    sample_origin.index = sample_origin.index.astype(str)

    if EPIC4:
        mapping = pd.read_csv(base_dir / "sample_sheet_EPIC4.csv", index_col=0)
    else:
        mapping = pd.read_csv(base_dir / "SWEPIC_full_sample_sheet.csv", index_col=0)

    idx = (
        mapping["Sentrix_ID"].astype(str) + "_" + mapping["Sentrix_Position"]
    ).ravel()

    mapping.index = idx
    mapping = mapping["Sample_Name"].to_dict()

    b_val_dir = base_dir / "beta_values"

    b_values = []
    for f in tqdm(b_val_dir.iterdir()):
        b_values.append(pd.read_pickle(f))

    b_values = pd.concat(b_values, axis=1)
    b_values = b_values.dropna()

    b_values = b_values.rename(columns=mapping)
    b_values = b_values.T
    b_values.index = b_values.index.astype(str)

    common_specimens = sample_origin.index.intersection(b_values.index)
    sample_origin = sample_origin.loc[common_specimens]
    b_values = b_values.loc[common_specimens]

    clinical = pd.read_csv(clinical_path)

    targets = pd.read_csv(target_path, dtype=int)

    excluded_patients_right = targets[
        (targets["Polyps (right)"] != 0) & (targets["Adenoma (right)"] == 0)
    ]

    excluded_patients_right = excluded_patients_right["Patient ID"].ravel()

    specimens_to_exclude_right = sample_origin[
        sample_origin.patient_id.isin(excluded_patients_right)
    ].index

    mapping_right = targets.set_index("Patient ID")["Adenoma (right)"].to_dict()

    mapping_right[135] = np.nan
    mapping_right[875] = np.nan

    ad_right = sample_origin.patient_id.replace(mapping_right).astype("category")
    ad_right.name = "Adenoma (right)"

    if EPIC4:
        EPIC4_b, EPIC4_clin, EPIC4_samples, EPIC4_phenotypes = get_EPIC_data(
            sample_origin=sample_origin,
            EPIC_name="EPIC4",
            specimens_to_exclude_right=specimens_to_exclude_right,
            b_values=b_values,
            ad_right=ad_right,
            clinical=clinical,
            bad_probes=bad_probes,
        )

        return EPIC4_b, EPIC4_clin, EPIC4_samples, EPIC4_phenotypes
    else:
        EPIC2_b, EPIC2_clin, EPIC2_samples, EPIC2_phenotypes = get_EPIC_data(
            sample_origin=sample_origin,
            EPIC_name="EPIC2",
            specimens_to_exclude_right=specimens_to_exclude_right,
            b_values=b_values,
            ad_right=ad_right,
            clinical=clinical,
            bad_probes=bad_probes,
        )
        EPIC3_b, EPIC3_clin, EPIC3_samples, EPIC3_phenotypes = get_EPIC_data(
            sample_origin=sample_origin,
            EPIC_name="EPIC3",
            specimens_to_exclude_right=specimens_to_exclude_right,
            b_values=b_values,
            ad_right=ad_right,
            clinical=clinical,
            bad_probes=bad_probes,
        )
        return (
            EPIC2_b,
            EPIC2_clin,
            EPIC2_samples,
            EPIC2_phenotypes,
            EPIC3_b,
            EPIC3_clin,
            EPIC3_samples,
            EPIC3_phenotypes,
        )
