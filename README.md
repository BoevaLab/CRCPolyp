# CRCPolyp
Supporting code for "DNA methylation variability in cecal mucosa of patients with tubular adenomas as a marker of early field cancerization"

## Installing packages to reproduce analysis

To reproduce the analysis, one needs to have the right packages installed. 
An easy way to do so it through [conda](https://www.anaconda.com/). 
With conda installed (see their website if not installed), just run 

```
conda create --name crcpolypenv --file requirements.txt
```

to create the env crcpolypenv that will contain all necessary packages to reproduce this analysis.

## Description

### FinalCode

For more details on folders, cf the respective READMEs.

##### aDVMC 

This folder contains supporting for the analysis of adenoma-related differentially variable and differentially methylated CpG sites (aDVMCs). 

##### download

This folder contains helper functions for downloading and preprocessing data.

##### utils

This folder contains helper functions for plotting 

##### notebooks

This folder contains all notebooks necessary for reproducing the analysis as well as a [detailed explanation](https://github.com/BoevaLab/CRCPolyp/blob/main/FinalCode/notebooks/README.md) of how to reproduce the figures and analysis presented in the manuscript.

## Questions

Any questions on the analysis or the code can be addressed to josephine.yates@inf.ethz.ch

## References

TBD




