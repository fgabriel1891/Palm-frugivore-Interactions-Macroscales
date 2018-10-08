---
output:
  pdf_document: default
  html_document: default
---

# Source code to replicate the analisis and figures presented in the article:

## A synthesis of animal-mediated seed dispersal of palms reveals distinct biogeographic differences in species interactions  *The Journal of Biogeography. 2018.(in revision)* 

### Gabriel Muñoz<sup>1</sup>,*, Kristian Trøjelsgaard<sup>2</sup> & W. Daniel Kissling <sup>1</sup>

#### 1: Institute for Biodiversity and Ecosystem Dynamics (IBED), University of Amsterdam, P.O. Box 94248, 1090 GE Amsterdam, The Netherlands

#### 2: Faculty of Engineering and Science, Department of Chemistry and Bioscience, Section of Biology and Environmental Science, University of Aalborg. Denmark

----------

### 1) Structure of this repository 

```
/JBI_Palms_REPO
|-- /AppendixTables
|   |-- MetadataAppendix1.pdf 
|   |-- MetadataAppendix2.pdf
|   |-- Appendix1-DataSources.xls # References to datasources for palm-frugivore interactions included in this study
|   |-- Appendix2-IndividualPalmSC.xls # Table of Sampling completeness per palm species. 
|-- /DATA
|   |-- PalmFrugDatasetOCT2018.csv # Palm-frugivore Interactions dataset
|   |-- datasetMetatada.pdf # Metadata of PalmFrugDatasetOCT2018.csv
|   |-- world.checklist.csv # Check list of palm occurrences by TDWG level 3 botanical countries
|   |-- BotanicalCountries # Shapefile of TDWG level 3 botanical countries
|   |   |-- TDWG_level3_Coordinates.dbf
|   |   |-- TDWG_level3_Coordinates.prj
|   |   |-- TDWG_level3_Coordinates.sbn
|   |   |-- TDWG_level3_Coordinates.sbx
|   |   |-- TDWG_level3_Coordinates.shp
|   |   |-- TDWG_level3_Coordinates.shp.xml
|   |   |-- TDWG_level3_Coordinates.shx
|-- /figs # Raw figures created from scripts
|   |-- ChordDiagramAfr.eps
|   |-- ChordDiagramNeo.eps
|   |-- InteractionsVsStudies.eps
|   |-- SamplingCompletenessA.eps
|   |-- SamplingCompletenessB.eps
|   |-- TraitMatching.eps
|-- /Scripts_R # Scripts to replicate the contents of the article
|   |-- Beckett2016LPA_wb_plus.R # Beckett's Modularity function
|   |-- CustomFunctions.R # Custom functions written for this article
|   |-- ReplicateAnalisis.R # Script to replicate analisis 
|   |-- ReplicateFigures.R # Script to produce the figures in /figs
|   |-- shiny.R # Shiny app to dynamically visualize sampling accumulation curves per species 
|-- JBI_Palms_REPO.Rproj # R container folder
|-- README.md
|-- README.html
|-- README.pdf

```

### 2) Instructions 

1) Open the R container folder 

2) Run ReplicateAnalisis.R and ReplicateFigures.R scripts (**other scripts are sourced within**)

**NOTE**:  Modularity replicates are TIME CONSUMING, currently those are commented in ReplicateAnalisis.R // If you want to re-run them and see the output uncomment those first. 

-------------






