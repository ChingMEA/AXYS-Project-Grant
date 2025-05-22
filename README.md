# AXYS-Project-Grant
R code used to analyze the data in "Persistent organic pollutant concentrations in human pancreas and peripancreatic adipose tissues correlate with markers of beta cell dysfunction".

Abbreviations:
- POP: persistent organic pollutants
- GSIS: glucose-stimulated insulin secretion

File descriptions:
1. "2025-05-22_AXYS_correlation_analyses_PCA.Rmd": Code used to run principal component analyses (PCA) and identify influential donors.
2. "2025-05-22_AXYS_correlation_analyses_imputing_missing_values.nb.Rmd": Code used to impute missing POP concentration values.
3. "2025-05-22_AXYS_correlation_analyses_donor_characteristics.nb.Rmd": Code used to correlate POP concentrations with donor characteristics.
4. "2025-05-22_AXYS_correlation_analyses_GSIS_parameters.nb.Rmd": Code used to correlate POP concentrations with GSIS parameters.
5. "2025-05-22_AXYS_sensitivity_analyses.nb.Rmd": Code used to compare correlation values after using (1) imputed vs non-imputed POPs datasets, (2) including all donors vs removing influencial donor(s).
6. "2025-05-22_AXYS_partial_correlations.Rmd": Code used to perform partial correlations, adjusting for age and/or BMI.
