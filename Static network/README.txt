NETWORK ANALYSIS OF GRANULARITY, SUMMARY OF MATLAB SCRIPTS
UPDATED 3 APRIL 2019

assign_subject_communities_loop.m: For every subject in the specified data set, this script assigns an initial community affiliation vector based on valence (i.e., 2 communities) and determines modularity values at 3 pre-specified connection weight thresholds. Across all subjects in the data set, this script generates a spreadsheet of the modularity at each threshold, as well as histograms of the distributions.

check_connection_weights.m: This script generates sub-plotted histograms of raw connection weights (i.e., correlation values) for every subject in the given data set. Subplot titles (subject IDs) that appear in red indicate that at least one node was removed from the correlation matrix due to lack of variance.

check_construct_validity_alexithymia.m: This script can be used to check construct validity for network metrics against TAS-20 scores.

check_construct_validity_complexity.m: This script can be used to check construct validity for network metrics against RDEES scores.

check_construct_validity_emodiversity.m: This script can be used to check construct validity for network metrics against emodiversity.

check_construct_validity_granularity.m: This script can be used to check construct validity for network metrics against the traditional ICC-based measure of granularity.

check_criterion_validity_18.m: This script can be used to check criterion validity for network metrics on the 18-term dataset.

check_criterion_validity_39.m: This script can be used to check criterion validity for network metrics on the 39-term dataset.

check_subject_clustering_loop.m: This script generates, for every subject in the specified data set, a spreadsheet of community assignment vectors and modularity values at 3 pre-specified connection weight thresholds. Across all subjects in the data set, this script generates a spreadsheet of the number of communities and modularity at each threshold, as well as histograms of the distributions for these at each threshold.

check_subject_communities.m: This script generates, for a given subject in a given data set, the community assignment vector at each of 3 pre-specified connection weight thresholds.

check_subject_communities_loop.m: This script generates, for every subject in the specified data set, a spreadsheet of community assignment vectors at 3 pre-specified connection weight thresholds. Across all subjects in the data set, this script generates a spreadsheet of the number of communities at each threshold, as well as histograms of the distribution at each threshold.

check_subject_data.m: This script generates, for a given subject in a given data set, a histogram of raw connection weights (i.e., correlation values), a heatmap of the correlation matrix, and sub-plotted histograms of raw intensity ratings for every emotion term. Plot titles (subject IDs) that appear in red indicate that at least one node was removed from the correlation matrix due to lack of variance.

check_subject_data_loop.m: This script generates, for every subject in the specified data set, a histogram of raw connection weights (i.e., correlation values), a heatmap of the correlation matrix, and sub-plotted histograms of raw intensity ratings for every emotion term. Plot titles (subject IDs) that appear in red indicate that at least one node was removed from the correlation matrix due to lack of variance.

check_subject_density_loop.m: This script generates, across all subjects in the data set, weighted and unweighted network density at each of 3 pre-specified connection weight thresholds, as well as histograms of the distribution at each threshold.

check_subject_participation_loop.m: This script generates, for every subject in the specified data set, a spreadsheet for 3 participation coefficients per node, each at 3 pre-specified connection weight thresholds. Across all subjects in the data set, this script generates 3 spreadsheets of the participation coefficients at each threshold, as well as histograms of the distribution for each coefficient at each threshold.

emodiversity_calculation.m: This script calculates emodiversity for all subjects in the specified data set.

gamma_tuning_crossvalidation.m: This script runs gamma tuning with nFold cross-validation for the specified data set, and generates plots for evaluation. Optional parameters: exclude negative connection weights, binarize connection weights, threshold connection weights, and change symmetrical treatment of negative v positive connection weights.

gamma_tuning_no_partitioning.m: This script runs gamma tuning WITHOUT crossvalidation for the specified data set, and generates plots for evaluation. NOTE: This script is no longer being actively maintained.

generate_subject_data.m: This script generates, for a given subject in a given data set, the correlation matrix and community assignment vector for use in visualization (e.g., in Gephi). Optional parameters: binarize or threshold connection weights, etc.

granularity_calculations.m: This script calculates traditional, ICC-based granularity measures for all subjects in the specified data set.




