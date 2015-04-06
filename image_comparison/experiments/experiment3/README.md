### Steps to do experiment 3

- First summarize hcp data (summarize_hcp.py)
- Prepare meta data files and logs for analysis (prep_make_group_maps.py)
- Generate group maps (make_group_maps.py,run_make_group_maps.py)
- Generate file list of group maps (prep_test_thresholding.py)
- Calculate qa for the maps (run_qa.py)
- Generate pickle results for spearman and pearson at different thresholds (test_thresholding.py,run_test_thresholding.py)
- Compile results (compile_threshold_results.py)
- Analysis in R (threshold_analysis.R, helper_functions.R)

### Extras
- find_missing.py will identify thresholds and directions to re-run
- filter_contrasts.py will remove "neg" contrasts, or inverse of the same (eg, faces-shapes == shapes-faces)
