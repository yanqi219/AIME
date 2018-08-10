# AIME
## Workflow
1) Feature extration/deconvolution: apLCMS/xMSanalyzer/XCMS <br />
   Batch effect correction: ComBat
2) Data import: HILIC_import_Qi
3) Split data into groups: HILIC_split_Qi
   - split into expo/unexpo or case/control groups, contains all features, save as "/panda_input/HILIC_classlabels_GROUP.txt" and "HILIC_ftrsmzcalib_combat_GROUP.txt".
   - Also need a combined dataset, contains all features and all comparison groups (i.e., unexposed case/control and exposed case/control), save as "HILIC_ftrsmzcalib_combat_ordered_all.txt". And the residuals of all features, save as "HILIC_all_residual.RData"--- the reason is to compare the intensity of metabolites between different groups (exposure vs. disease), no need if this comparison is not required.
4) Annotate features: annotation_cluster
   - Use "/Panda_HILIC/HILIC_ftrsmzcalib_combat_ordered_all.txt" and "/Panda_HILIC/HILIC_classlabels_for_panda_all.txt" as input. Summarize the triplicate first, no filter, save as "HILIC_annotation_input.RData". And then use it to do annotation.
   - Match "HILIC_annotation_input.RData" which is just a list of all features (mz/time) with in-house library for the verification. Save as "/HILIC_annotation/HILIC_annotation_verified.txt"
   - Almost no difference if we use different dataset (all subjects, air pollution subjects, subgroups) as long as all features are included. But theoretically should we use all subjects?
5) Adjust covariates: HILIC_CovAdjust_Comp3_qy
   - Need to consider different set of covariates
   - save as "HILIC_GROUP_residual_nonorm_WGCNA.RData" (residual) and "HILIC_GROUP_classification_nonorm.RData" (non residual)
6) Feature selection: HILIC_classification_qy
   - Will remove some part in the future such as "IV.Annotation" which is out of date
   - save two outputs: one is the feature selection results, another one is the mummichog input
   - Then if needed, get the correlated metabolites using MetabNet and save as another mummichog input
   - Can also draw box plot for pathways/modules --- not needed anymore
7) More feature selection: MetBoruta
   - Use RF to select features, can be compared with PLSDA using getVenn
7) Pathway analysis: mummichog
   - Depends on which version of mummichog, have different format of output
8) Visualization: HILIC_classification_qy/pathway_Plotly/Cytoscape/pathview
   - Find this function in "C18_classification_Comp1_qy", called "Activity network for cytoscape/KEGG mapper"
   - for pathview: right now it only plot significant features annotated by mummichog, may want to change to all features annotated by mummichog
   - Cytoscape can use mummichog files in "/sif/" as input
9) More visualization: Draw Boxplot for Poster
   - Used to draw barplot for pathways
10) Create final table: Create_final_table
    - Combine results from feature selection, xMSannotator, mummichog, and in-house library verification together
