Phase C — Machine Learning Prediction Model

This folder contains data and scripts for the Phase C 
machine learning analysis (CN-to-MCI conversion prediction).

 Files
 ml_baseline_data.csv: Baseline data for 217 cognitively normal (CN) participants, including neuroanatomical measures, 
  biomarkers, and conversion outcome (converted to MCI/Dementia at any point during follow-up).

 Key Details
 N = 217 CN participants at baseline
 22 converted to MCI/Dementia during follow-up 
 Outcome variable: 'converted' (1 = converted, 0 = did not convert)

Known Limitation
With only 22 positive events, this analysis is exploratory/
hypothesis-generating, not a validated predictive model. 
Cross-validation should use 3 folds given the small 
number of events. Results should be framed accordingly in any 
manuscript.

 Analysis Pipeline
See PhaseC_Pipeline_Python.md in the project root for the complete Python-based ML pipeline: 
data loading, missing data checks, Model A (neuroanatomy + 
demographics) vs Model B (+ biomarkers), cross-validation, 
and SHAP feature importance.

