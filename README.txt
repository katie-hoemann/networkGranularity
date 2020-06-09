ATTENTION!
> dataSet 18 (and any similar file names) refers to SAMPLE 1
> dataSet 39 (and any similar file names) refers to SAMPLE 2

ARCHIVE
> Scripts and files that are either outdated or unused in final analyses; stored here as back-up

FIGURES
> Networks: files necessary for creating network visualizations in Gephi, including style guide
> Scatter plots: MATLAB script and files necessary for creating scatter plots of results

FUNCTIONS
> MATLAB functions that are necessary for running scripts
> Brain Connectivity Toolbox (BCT) functions can be downloaded here: https://sites.google.com/site/bctnet/

GENERAL
> MATLAB scripts used to compare output data or compute data set statistics
check_network_measure_correlations.m = compares network measures; generates correlation matrix, heat map, scatter plots
check_self-report_measure_correlations_18.m = compares self-report and other non-network derived measures for 18-term data set; generates correlation matrix, heat map, histograms, scatter plots
check_self-report_measure_correlations_39.m = compares self-report and other non-network derived measures for 39-term data set; generates correlation matrix, heat map, histograms, scatter plots
dataset39term_stats.m = calculates number of days and prompts statistics for 39-term data set
granularity_calculations_daily.m = calculates estimates of granularity for each day in a given data set

INPUT
> Data files that are necessary for running scripts
18ARIEOD.xlsx = raw data without day variables for 18-term data set; used in static networks
18ARIEOD_daily_filtered.xlsx = raw data including day, time, and time of day variables for 18-term data set; 'No Late Surveys' worksheet used in time-varying networks (note: 'filtered' refers to variables, not cases)
18ARIVariations_QuestionnaireData.xlsx = self-report questionnaire data for 18-term data set
39ESdata.xlsx = raw data without day variables for 39-term data set; used in static networks
39ESdata_daily_filtered.xlsx = raw data including day, time, and time of day variables for 39-term data set; used in time-varying networks (note: 'filtered' refers to variables not cases)
39ESdata_QuestionnaireData.xlsx = self-report questionnaire data for 39-term data set
words18.csv = list of words for 18-term data set, with valence and arousal ratings from Warriner et al (2014)
words39.csv = list of words for 39-term data set, with valence and arousal ratings from Warriner et al (2014)

REGRESSION ANALYSES
> SPSS databases and output files used to conduct multiple regressions of network, derived, and self-report measures
> Excel spreadsheets used for conducting mini-meta-analyses of regression results

STATIC NETWORK
> MATLAB scripts used to compute measures from all experience sampling data (i.e., overall networks)
derived_measures_calculations.m = calculates non-network measures on all data
generate_subject_data.m = generates, for a given subject in a given data set, the correlation matrix and community assignment vector for use in visualization software
network_measures_calculations_Pearson_FINAL.m = calculates measures for networks estimated using Pearson correlations (final version)

TIME-VARYING NETWORK
> MATLAB scripts used to compute measures from 3-day overlapping windows of experience sampling data (i.e., time-varying networks)
derived_measures_calculations_windows.m = calculates non-network measures on windowed data
generate_subject_data_windows.m = generates, for a given subject in a given data set, correlation matrices and community assignment vectors for use in visualization software
network_measures_calculations_windows_Pearson.m = calculates means and standard deviations of measures for time-varying networks estimated using Pearson correlations (final version)


