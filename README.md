# POMI

Title: Using multiple imputation to classify potential outcomes subgroups

Authors: Yun Li, Irina Bondarenko, Michael R. Elliott, Timothy P. Hofer and Jeremy M.G. Taylor  

Abstract: With medical tests becoming increasingly available, concerns about over-testing, over-treatment and health care cost dramatically increase. Hence, it is important to understand the influence of testing on treatment selection in general practice. Most statistical methods focus on average effects of testing on treatment decisions. However, this may be ill-advised, particularly for patient subgroups that tend not to benefit from such tests. Furthermore, missing data are common, representing large and often unaddressed threats to the validity of most statistical methods. Finally, it is often desirable to conduct analyses that can be interpreted causally. Using the Rubin Causal Model framework, we propose to classify patients into four potential outcomes subgroups, defined by whether or not a patient's treatment selection is changed by the test result and by the direction of how the test result changes treatment selection. This subgroup classification naturally captures the differential influence of medical testing on treatment selections for different patients, which can suggest targets to improve the utilization of medical tests. We can then examine patient characteristics associated with patient potential outcomes subgroup memberships. We used multiple imputation methods to simultaneously impute the missing potential outcomes as well as regular missing values. This approach can also provide estimates of many traditional causal quantities of interest. We find that explicitly incorporating causal inference assumptions into the multiple imputation process can improve the precision for some causal estimates of interest. We also find that bias can occur when the potential outcomes conditional independence assumption is violated; sensitivity analyses are proposed to assess the impact of this violation. We applied the proposed methods to examine the influence of 21-gene assay, the most commonly used genomic test in the United States, on chemotherapy selection among breast cancer patients.

----------------------------------------------------------------------------------------------------------------------------------------
Our methods can be implemented through IVEware in SAS, and through the MICE (Multivariate Imputation by Chained Equations) package in R.

For example codes in SAS, see these files in the folder: imputation_example_annotated_YL, ex_stage0, ex_chemo0, ex_chemo1, ex_rest 

For example code in R, see the file in the folder: pomi_ind_YL
