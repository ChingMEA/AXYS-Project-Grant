# AXYS-Project-Grant
This repository contains the full R analysis pipeline used to examine how persistent organic pollutant (POP) concentrations (dioxins, PCBs, OCPs) in human pancreas and peripancreatic adipose tissues relate to markers of islet and β-cell function. The results generated from this pipeline are reported in the following publication from the Bruin Lab:  
*“Persistent organic pollutant concentrations in human pancreas and peripancreatic adipose tissues correlate with markers of beta cell dysfunction”*

### Analyses included in this repository
- **Within-tissue POP–POP correlations (pancreas, adipose)**
- **POP–GSIS correlations across multiple functional readouts:**
  - Static GSIS
  - Glucose perifusion
  - Leucine perifusion
  - Fatty acid perifusion
- **Sensitivity analyses:**
  - Imputation of missing POP concentrations
  - Removal of influential donors (singly and as a group)
- **Correlations with pollutant class sums (summed dioxins, summed PCBs)**

### Abbreviations:
- **POP**: Persistent organic pollutants
- **GSIS**: Glucose-stimulated insulin secretion
- **PCB**: Polychlorinated biphenyls
- **OCP**: Organochlorine pesticides

### File descriptions:
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

**Note:** All files references in the pipeline are found in **/data**.
