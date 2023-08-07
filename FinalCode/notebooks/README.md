## Description

This folder contains the notebooks used for the analysis. 

### How to reproduce the analysis

Once can reproduce the analysis presented in the paper by running the notebooks in order, after having downloaded the necessary data as described in step 0. 
To reproduce the analysis from A to Z, it is important to run the notebooks in the order presented here, as the notebook might be using some files generating in one of the previous ones. All the intermediate files generated here are however present in the Zenodo shared data, which can be used if only a part of the analysis needs to be reproduced.
In general, placeholders in the notebooks are indicated with "/add/path/here". They are also described for each notebook. 

#### General comments
In the code in general, you will see the notation adVMP - this is equivalent to aDVMC. You will also see the notations EPIC2, EPIC3, EPIC4 - these are equivalent to SWEPIC1, SWEPIC2, SWEPIC3.

#### 0. Clone the repository and download the data

First, clone this repository in a folder of your choice. Navigate to the notebooks folder to run the analyses.

> **_NOTE_**: If you wish to run the notebooks from a separate folder, don't forget to change the sys.path to the folder containing the code (FinalCode/ folder) to be able to download the helper functions.

##### In house data
You can find the preprocessed data at XXXXX  (TODO)
You can find the clinical and auxiliary data at XXXX (TODO). 

ADD IN THE GENE EXPRESSION

##### External data
For external methylation datasets, we use the [methylprep package](https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/). 

Install the package following the installation guide on their website, then run the following command, replacing GEOID by the ID of the GEO database to download (i.e., GSE132804, GSE48684, GSE199057), and the path to where you want to save the data.

```
python -m methylprep -v download -i {GEOID} -d "path/to/save/methylprep_GSE132804" --batch_size 50
```

For the external gene expression data, GSE76987, navigate to [this link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76987) and download the following files

```
GSE76987_ColonCancerProcessed.csv
GSE76987_RightCancerProcessed.csv
```

Save these files in the directory of your choice.


#### 1. `GeneralDatasetCharacteristics.ipynb`

This notebook runs the analysis to reproduce Fig 1, Suppl. Table S1 and Fig S1.

##### Placeholders

- fig_dir = pl.Path("/add/path/here/"); where the figures will be automatically saved
- base_dir = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC1/2 (see 0. Download the data) is stored 
- base_dir4 = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC3 (see 0. Download the data) is stored.
- data_dir = pl.Path("/add/path/here/"); the folder where the clinical and auxiliary data (see 0. Download the data) is stored.

#### 2. `adVMP-DMR.ipynb`

This notebook computes the comparison between differential methylation and differential variability presented in Suppl. Fig S2.

##### Placeholders

- fig_dir = pl.Path("/add/path/here/"); where the figures will be automatically saved
- base_dir = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC1/2 (see 0. Download the data) is stored 
- base_dir4 = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC3 (see 0. Download the data) is stored.
- data_dir = pl.Path("/add/path/here/"); the folder where the clinical and auxiliary data (see 0. Download the data) is stored.
- result_dir = pl.Path("/add/path/here"); the folder where the results of the analysis will be saved. 

#### 3. `adVMPDeconvLink.ipynb`

This notebook uses the EpiSCORE algorithm to get an estimate of the tissue composition, used in Suppl. Table S3

> **_NOTE_**: You will have to run this notebook twice, once for SWEPIC1/2, once for SWEPIC3.

##### Placeholders

- data_dir = pl.Path("/add/path/here/"); the folder where the clinical and auxiliary data (see 0. Download the data) is stored.
- mapping_path = pl.Path("add/path/here/"); the folder where the sample sheet (of resp. SWPEPIC1/2 or SWEPIC3) is stored.
- base_dir = pl.Path("add/path/here/sesame_processed_EPIC"); the folder where the processed DNAmeth data for SWEPIC1/2 or SWEPIC3 (see 0. Download the data) is stored.
- resdir = pl.Path("/add/path/here"); the folder where the results of the analysis will be saved.

#### 4. `adVMPDiscovery.ipynb`

This notebook computes the aDVMCs in the full cohorts and the results for Fig 2a, 2b, 2f and Suppl. Fig S6.

##### Placeholders

- fig_dir = pl.Path("/add/path/here/"); where the figures will be automatically saved
- base_dir = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC1/2 (see 0. Download the data) is stored 
- base_dir4 = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC3 (see 0. Download the data) is stored.
- data_dir = pl.Path("/add/path/here/"); the folder where the clinical and auxiliary data (see 0. Download the data) is stored.
- result_dir = pl.Path("/add/path/here"); the folder where the results of the analysis will be saved. 
- deconv_path  = pl.Path("/add/path/here/"); the folder where the results of step 3, the deconvolution results, are stored

#### 5. `adVMP-GlobalVsLocal.ipynb`

This notebook computes most variable probes in healthy tissue, results used in Suppl. Fig S5 and used to compare against for Suppl. Fig S3 and S7.

##### Placeholders

- fig_dir = pl.Path("/add/path/here/"); where the figures will be automatically saved
- base_dir = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC1/2 (see 0. Download the data) is stored 
- base_dir4 = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC3 (see 0. Download the data) is stored.
- data_dir = pl.Path("/add/path/here/"); the folder where the clinical and auxiliary data (see 0. Download the data) is stored.

#### 5bis. (Optional) `HorvathClockComparison.ipynb`

> **_NOTE_**: This notebook necessitates raw data not shared in the associated repository but that can be shared on demand. 
> The reason for this is that the probes used in the Horvath clock were mostly not present after our stringent quality control. 
> We thus included the interesection of all probes detected with the Horvath clock probes, regardless of QC metrics. 
> The result of this notebook is available in the Zenodo repository under `horvath_age.csv`

This notebook computes the Horvath clock age for each patient to be used in Suppl. Fig S4.

##### Placeholders

- data_dir = pl.Path("/add/path/here/"); the folder where the clinical and auxiliary data (see 0. Download the data) is stored.
- PATH1 = "/add/path/here"; the folder where the raw data for SWEPIC1/2 is stored.
- PATH2 = "/add/path/here"; the folder where the raw data for SWEPIC3 is stored.
- sheet_dir = pl.Path("/add/path/here/"); the folder where the sample sheets are stored 

#### 6. `adVMPCrossVal.ipynb`

This notebook computes the cross-validated analysis on aDVMCs and results for Fig 2c, 2d, 2e and Suppl. Fig S3.

##### Placeholders

- fig_dir = pl.Path("/add/path/here/"); where the figures will be automatically saved
- base_dir = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC1/2 (see 0. Download the data) is stored 
- base_dir4 = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC3 (see 0. Download the data) is stored.
- data_dir = pl.Path("/add/path/here/"); the folder where the clinical and auxiliary data (see 0. Download the data) is stored.

#### 7. `adVMPExternalValidation.ipynb`

This notebook computes the analysis of aDVMC and the performance of aDVMC-based classifiers in external datasets, as shown in Fig 2g, 2h, and Suppl. Fig S7, S8, S9, S13.

##### Placeholders

- fig_dir = pl.Path("/add/path/here/"); where the figures will be automatically saved
- data_dir = pl.Path("/add/path/here/"); the folder where the clinical and auxiliary data (see 0. Download the data) is stored.
- data_dir_GSE199057 = pl.Path("/add/path/here/"); the folder where the results of methylprep (see 0. Download the data) for GSE199057 is stored.
- data_dir_GSE132804 = pl.Path("/add/path/here/"); the folder where the results of methylprep (see 0. Download the data) for GSE132804 is stored.
- data_dir_GSE48684 = pl.Path("/add/path/here/"); the folder where the results of methylprep (see 0. Download the data) for GSE48684 is stored.

#### 8. `adVMPEnrichmentRoadmaps.ipynb`

This notebook computes the enrichment of aDVMC in the [Roadmap Epigenomics Consortium](https://www.nature.com/articles/nature14248) annotated regions, shown in Fig 3a and 3c.

##### Placeholders

- fig_dir = pl.Path("/add/path/here/"); where the figures will be automatically saved
- base_dir = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC1/2 (see 0. Download the data) is stored 
- base_dir4 = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC3 (see 0. Download the data) is stored.
- data_dir = pl.Path("/add/path/here/"); the folder where the clinical and auxiliary data (see 0. Download the data) is stored.

#### 9. `adVMPWholeCGIMeth.ipynb`

This notebook computes the analysis of the local vs regional dysregulation at aDVMCs, shown in Fig 3c and Suppl. Fig S10.

##### Placeholders

- fig_dir = pl.Path("/add/path/here/"); where the figures will be automatically saved
- base_dir = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC1/2 (see 0. Download the data) is stored 
- base_dir4 = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC3 (see 0. Download the data) is stored.
- data_dir = pl.Path("/add/path/here/"); the folder where the clinical and auxiliary data (see 0. Download the data) is stored.

#### 10. `adVMPClinicalLifestyleLink.ipynb`

This notebook computes the link between aDVMC methylation and clinical and lifestyle characteristics, as shown in Fig 3d and Suppl. Fig S11.

##### Placeholders

- fig_dir = pl.Path("/add/path/here/"); where the figures will be automatically saved
- base_dir = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC1/2 (see 0. Download the data) is stored 
- base_dir4 = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC3 (see 0. Download the data) is stored.
- data_dir = pl.Path("/add/path/here/"); the folder where the clinical and auxiliary data (see 0. Download the data) is stored.

#### 11. `adVMPLinkGEX.ipynb`

This notebook computes the analysis of aDVMC-related genes in terms of gene expression in an external dataset, as shown in Fig 3e and 3f, and Suppl. Fig S12.

##### Placeholders

- fig_dir = pl.Path("/add/path/here/"); where the figures will be automatically saved
- base_dir = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC1/2 (see 0. Download the data) is stored 
- base_dir4 = pl.Path("/add/path/here/"); the folder where the processed DNAmeth data for SWEPIC3 (see 0. Download the data) is stored.
- data_dir = pl.Path("/add/path/here/"); the folder where the clinical and auxiliary data (see 0. Download the data) is stored.
- gex_dir = pl.Path("/add/path/here"); the folder where the data from GSE76987 (see 0. Download the data) is stored.
- resdir = pl.Path("/add/path/here"); the folder where the results of the analysis will be saved.