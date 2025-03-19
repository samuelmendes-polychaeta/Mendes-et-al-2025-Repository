
# Overview

This repository contains all data, script and functions used to calculate the Standardized Effect Sizes (SES) of MPD and MNTD indices, as well as the RLQ and fourthcorner analysis outputs from Mendes et al. (2025). Importantly, we present all fuzzy coded trait information, with the references that justify the scores assigned to each trait modality.

We demonstrated that South Atlantic Central Water (SACW) upwelling and environmental filtering are important factors influencing the functional trait dispersion of macrobenthic polychaete assemblages.

# Repository structure

## Folders 

- Data: The folder presents a single excel file called `data.xlsx` that contains all information needed to reproduce the results from Mendes et al. (2025). This file is organized in six sheets, but only three of them should be imported to R. They are: `R_trait_df`, `abund`, and `env_variables`. The information of each sheet is detailed below.
    + `R_trait_df` contains all fuzzy coded traits. The trait modalities are abbreviated for facilitating the RLQ output visualization. Each abbreviation can be consulted in `trait_metadata` sheet. 
    + `refs` contains all references that we consulted to justify the scores assigned to each trait modality. 
    + `trait_data` contains the fuzzy trait modalities, without abbreviations, and their respective categories.
    + `trait_metadata` summary of all trait information, including their abbreviations and identification number. 
    + `abund` contains polychaete genera incidences. 
    + `env_variables` contains the abiotic variables. 
 
      
- R: this folder contains the `Script` for calculation of SES MPD, MNTD, RLQ and fourthcorner analysis 


  
# Downloading the repository

The user can download this repository to a local folder in your computer 

## downloading all files

```{r eval=FALSE, echo=TRUE}
download.file(url = "https://github.com/samuelmendes-polychaeta/Mendes-et-al-2025-Repository/archive/refs/heads/main.zip", destfile = "Mendes_et_al_2025.zip")
```


to unzip the .zip file in your computer type

```{r eval=FALSE,echo=TRUE}
unzip(zipfile = "Mendes_et_al_2025.zip")
```


# Authors

Samuel Lucas da Silva Delgado Mendes, Paulo Cesar de Paiva, Rodolfo Leandro Nascimento

