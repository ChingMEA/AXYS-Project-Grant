# AXYS-Project-Grant
This repository contains the full R analysis pipeline used to examine how persistent organic pollutant (POP) concentrations (dioxins, PCBs, OCPs) in human pancreas and peripancreatic adipose tissue relate to markers of islet and β-cell function. Results generated from this pipeline are reported in the following publication from the Bruin Lab:
“Persistent organic pollutant concentrations in human pancreas correlate with markers of beta cell dysfunction”

### Analyses included in this repository
- **Within-tissue POP–POP correlations (pancreas, adipose)**
- **POP–GSIS correlations across multiple functional readouts:**
  - Static GSIS
  - Glucose perifusion
  - Leucine perifusion
  - Fatty-acid perifusion
- **Sensitivity analyses:**
  - Imputation of missing POP concentrations
  - Removal of influential donors (singly and as a group)
- **Correlations with pollutant class sums (summed dioxins, summed PCBs)**
- **Inter-lab variability**

### Abbreviations:
- **POP**: Persistent organic pollutants
- **GSIS**: Glucose-stimulated insulin secretion
- **PCB**: Polychlorinated biphenyls
- **OCP**: Organochlorine pesticides
- **LOD**: Limit of detection

### Main AXYS flags:
- **ND**: Analyte not detected (signal below the detection limit).
- **NDR**: Analyte detected but not reliably quantifiable; signal is present but does not meet laboratory criteria for accurate quantification.

### Pollutant concentration data:
The POP concentration datasets shared in this repository are processed data. Values were blank-corrected, adjusted relative to the analyte-specific limit of detection (LOD), and filtered to remove analytes meeting the exclusion criteria described below. LOD adjustment was performed by replacing blank-corrected values ≤ 0 with 1/2 the LOD for the corresponding analyte.

1. **2024-08-24_pops_data_outliers_excluded_for_single_analyte_analyses.csv**
Analytes with %ND ≥ 50% and %NDR ≥ 40% were excluded.
2. **2024-08-24_pops_data_outliers_excluded_for_summed_analyses.csv**
Analytes with and %NDR ≥ 40% were excluded.

### Descriptions of each analysis file:
1. **2025-05-22_AXYS_correlation_analyses_imputing_missing_values.Rmd**  
Imputation of missing POP concentration values.
2. **2025-05-22_AXYS_sensitivity_analyses.Rmd**  
Sensitivity analyses comparing imputed vs. non-imputed datasets and inclusion vs. exclusion of influential donor(s).
3. **2025-07-30_AXYS_correlation_analyses_GSIS_parameters.Rmd**  
Correlations between POP concentrations and GSIS parameters.
4. **2025-07-30_AXYS_correlation_analyses_donor_characteristics.Rmd**  
Correlations between POP concentrations and donor metadata.
5. **2025-08-05_AXYS_correlation_analyses_PCA.Rmd**  
Principal component analyses (PCA) and identification of influential donors.
6. **2025-08-25_AXYS_partial_correlations.Rmd**  
Partial correlations adjusting for age and/or BMI.
7. **2025-11-21_AXYS_correlation_analyses_inter-lab_variability.Rmd**  
Quantifying inter-laboratory variability (Bruin Lab glucose perifusion vs. ADI IsletCore static GSIS).

**Note:** All input files required to run the pipeline are located in the **/data** directory.
