# *Code and Data for: Fire-Retardant Increases Post-Fire Soil N and P with Variable Effects on Post-Fire Vegetation One Year After Fire*

## Description of the project

This repository contains all code and data for *Drought may initiate western
spruce budworm outbreaks, but multi-year periods of increased moisture
availability promote widespread defoliation*.

## Organization of the project

The project has the following structure:

-   *.gitignore*: a plain text file that specifies intentionally
    untracked files and directories that Git should ignore

-   scripts: This subdirectory contains all code written for this project.
    In order for the code to work, files should be run sequentially
    (i.e., *1-Clean Environmental Data.R* then *2-Main Fire Retardant Analysis.R*).

    -   *1-Clean Environmental Data.R*: This code compiles environmental
        data.

    -   *2-Main Fire Retardant Analysis.R*: The code for reproducing the main
        results presented therein.

-   data: This subdirectory contains field and laboratory data used in this project.

    -   *CarterFR_analysis.csv*: Soil nutrient extractions from field collected soils
        with methods described in the manuscript. Units for all soil nutrient concentrations is
        provided in mg/g. The raw data (unformatted for R) are provided the .xlsx file with the 
        same name. 

    -   Variables

        -   Field ID: field plot name
        -   Location: field site name
        -   Burn Severity: dNBR burn severity
        -   Treatment: presence/absence of fire-retardant
        -   pH: soil pH
        -   DOC: Dissolved organic carbon
        -   TDN: Total disolved nitrogen
        -   Na: soil Na concentration (mg/g)
        -   NH4: soil NH4 concentration (mg/g)
        -   K: soil K concentration (mg/g)
        -   Mg: soil Mg concentration (mg/g)
        -   Ca: soil Ca concentration (mg/g)
        -   Cl: soil Cl concentration (mg/g)
        -   NO3: soil NO3 concentration (mg/g)
        -   PO4: soil PO4 concentration (mg/g)
        -   SO4: soil SO4 concentration (mg/g)

    -   *ENV_Data.csv*: Environmental data derived from script 1 and composed of
        Quarry_ENV.csv and StoneCanyon_ENV.csv with the same variable names. 

    -   Variables

        -   plot: field plot name
        -   site: field site name
        -   sev: dNBR burn severity
        -   trt: presence/absence of fire-retardant
        -   elev: elevation (m)
        -   slope: slope (degrees)
        -   tpi: topographic position index
        -   aspect: aspect

    -   *FieldData(comm).csv*: Raw data for vegetative cover. 

    -   Variables

        -   plot: field plot name
        -   transect: azimuth of transect within plot (0, 120, 240)
        -   subplot: location of Daubenmire quadrat along each transect
        -   code: six letter species code used to identify each species
        -   cov: percent cover of observation within the Duabenmire quadrat

    -   *FieldData(rich).csv*: Raw data for vegetative richness. 

    -   Variables

        -   plot: field plot name
        -   code: six letter species code used to identify each species
        -   species: latin name that corresponds to each code

    -   *IntroducedStatusKey.csv*: Supplementary vegetation data from USDA Plants. 

    -   Variables

        -   code: six letter species code used to identify each species
        -   species: latin name that corresponds to each code
        -   duration: life history (annual or perennial)
        -   functional group: functinoal group (shrub, forb, gram)
        -   status: native vs. introduced status (N, I)

-   Documents

-   *README.md:* this file

