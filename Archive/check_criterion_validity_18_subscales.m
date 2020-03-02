% This script can be used to check criterion validity for network metrics
% on the 18-term dataset

%% Instructions:
% run the following scripts for the 18-term dataset:
% - check_subject_communities_loop (set subjectFiles to 0)
% - check_subject_clustering_loop (set subjectFiles to 0)
% - check_subject_participation_loop (set subjectFiles to 0)
% - check_subject_density_loop
% - granularity_calculations
% - emodiversity_calculation
% make sure the necessary .mat files are on your path

clear;
clc;

%% specify dataset
dataSet = 18; % must be set to 18
print = 1; % set to 1 to generate spreadsheet with results

%% load data and list of checks
filenames{1} = ['dataSet' num2str(dataSet) '_zICC.mat'];
filenames{2} = ['dataSet' num2str(dataSet) '_communities.mat'];
filenames{3} = ['dataSet' num2str(dataSet) '_clustering.mat'];
filenames{4} = ['dataSet' num2str(dataSet) '_participation.mat'];
filenames{5} = ['dataSet' num2str(dataSet) '_gateway.mat'];
filenames{6} = ['dataSet' num2str(dataSet) '_diversity.mat'];
filenames{7} = ['dataSet' num2str(dataSet) '_density.mat'];
filenames{8} = ['dataSet' num2str(dataSet) '_emodiversity.mat'];
for i_file = 1:numel(filenames)
    load(filenames{i_file});
end
checks = readtable('criterionChecks_18_subscales.csv');

%% load criterion validity measures
rawData = importdata('18ARIVariations_QuestionnaireData.xlsx');
ASIphys = rawData.data(:,4)-7; % subtract 7 to bring in line with original rating scale
ASIcog = rawData.data(:,5)-7; % subtract 7 to bring in line with original rating scale
ASIsoc = rawData.data(:,6)-7; % subtract 7 to bring in line with original rating scale
ERSsens = rawData.data(:,8);
ERSarous = rawData.data(:,9);
ERSpers = rawData.data(:,10);
TASddf = rawData.data(:,15); %37 for session 2
TASdif = rawData.data(:,16); %38 for session 2
   
%% correlations between ASI physical and network metrics
[corr18(1,1),corr18(1,2)] = corr(ASIphys(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(2,1),corr18(2,2)] = corr(ASIphys(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(3,1),corr18(3,2)] = corr(ASIphys(:),networkDensity(:,3),'rows','complete');
[corr18(4,1),corr18(4,2)] = corr(ASIphys(:),networkDensity(:,4),'rows','complete');
[corr18(5,1),corr18(5,2)] = corr(ASIphys(:),participation(:,2),'rows','complete');
[corr18(6,1),corr18(6,2)] = corr(ASIphys(:),participation(:,3),'rows','complete');
[corr18(7,1),corr18(7,2)] = corr(ASIphys(:),participation(:,5),'rows','complete');
[corr18(8,1),corr18(8,2)] = corr(ASIphys(:),participation(:,6),'rows','complete');
[corr18(9,1),corr18(9,2)] = corr(ASIphys(:),diversity(:,2),'rows','complete');
[corr18(10,1),corr18(10,2)] = corr(ASIphys(:),diversity(:,3),'rows','complete');
[corr18(11,1),corr18(11,2)] = corr(ASIphys(:),diversity(:,5),'rows','complete');
[corr18(12,1),corr18(12,2)] = corr(ASIphys(:),diversity(:,6),'rows','complete');

%% correlations between ASI cognitive and network metrics
[corr18(13,1),corr18(13,2)] = corr(ASIcog(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(14,1),corr18(14,2)] = corr(ASIcog(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(15,1),corr18(15,2)] = corr(ASIcog(:),networkDensity(:,3),'rows','complete');
[corr18(16,1),corr18(16,2)] = corr(ASIcog(:),networkDensity(:,4),'rows','complete');
[corr18(17,1),corr18(17,2)] = corr(ASIcog(:),participation(:,2),'rows','complete');
[corr18(18,1),corr18(18,2)] = corr(ASIcog(:),participation(:,3),'rows','complete');
[corr18(19,1),corr18(19,2)] = corr(ASIcog(:),participation(:,5),'rows','complete');
[corr18(20,1),corr18(20,2)] = corr(ASIcog(:),participation(:,6),'rows','complete');
[corr18(21,1),corr18(21,2)] = corr(ASIcog(:),diversity(:,2),'rows','complete');
[corr18(22,1),corr18(22,2)] = corr(ASIcog(:),diversity(:,3),'rows','complete');
[corr18(23,1),corr18(23,2)] = corr(ASIcog(:),diversity(:,5),'rows','complete');
[corr18(24,1),corr18(24,2)] = corr(ASIcog(:),diversity(:,6),'rows','complete');

%% correlations between ASI social and network metrics
[corr18(25,1),corr18(25,2)] = corr(ASIsoc(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(26,1),corr18(26,2)] = corr(ASIsoc(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(27,1),corr18(27,2)] = corr(ASIsoc(:),networkDensity(:,3),'rows','complete');
[corr18(28,1),corr18(28,2)] = corr(ASIsoc(:),networkDensity(:,4),'rows','complete');
[corr18(29,1),corr18(29,2)] = corr(ASIsoc(:),participation(:,2),'rows','complete');
[corr18(30,1),corr18(30,2)] = corr(ASIsoc(:),participation(:,3),'rows','complete');
[corr18(31,1),corr18(31,2)] = corr(ASIsoc(:),participation(:,5),'rows','complete');
[corr18(32,1),corr18(32,2)] = corr(ASIsoc(:),participation(:,6),'rows','complete');
[corr18(33,1),corr18(33,2)] = corr(ASIsoc(:),diversity(:,2),'rows','complete');
[corr18(34,1),corr18(34,2)] = corr(ASIsoc(:),diversity(:,3),'rows','complete');
[corr18(35,1),corr18(35,2)] = corr(ASIsoc(:),diversity(:,5),'rows','complete');
[corr18(36,1),corr18(36,2)] = corr(ASIsoc(:),diversity(:,6),'rows','complete');

%% correlations between ERS sensitivity and network metrics
[corr18(37,1),corr18(37,2)] = corr(ERSsens(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(38,1),corr18(38,2)] = corr(ERSsens(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(39,1),corr18(39,2)] = corr(ERSsens(:),networkDensity(:,3),'rows','complete');
[corr18(40,1),corr18(40,2)] = corr(ERSsens(:),networkDensity(:,4),'rows','complete');
[corr18(41,1),corr18(41,2)] = corr(ERSsens(:),participation(:,2),'rows','complete');
[corr18(42,1),corr18(42,2)] = corr(ERSsens(:),participation(:,3),'rows','complete');
[corr18(43,1),corr18(43,2)] = corr(ERSsens(:),participation(:,5),'rows','complete');
[corr18(44,1),corr18(44,2)] = corr(ERSsens(:),participation(:,6),'rows','complete');
[corr18(45,1),corr18(45,2)] = corr(ERSsens(:),diversity(:,2),'rows','complete');
[corr18(46,1),corr18(46,2)] = corr(ERSsens(:),diversity(:,3),'rows','complete');
[corr18(47,1),corr18(47,2)] = corr(ERSsens(:),diversity(:,5),'rows','complete');
[corr18(48,1),corr18(48,2)] = corr(ERSsens(:),diversity(:,6),'rows','complete');

%% correlations between ERS arousal and network metrics
[corr18(49,1),corr18(49,2)] = corr(ERSarous(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(50,1),corr18(50,2)] = corr(ERSarous(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(51,1),corr18(51,2)] = corr(ERSarous(:),networkDensity(:,3),'rows','complete');
[corr18(52,1),corr18(52,2)] = corr(ERSarous(:),networkDensity(:,4),'rows','complete');
[corr18(53,1),corr18(53,2)] = corr(ERSarous(:),participation(:,2),'rows','complete');
[corr18(54,1),corr18(54,2)] = corr(ERSarous(:),participation(:,3),'rows','complete');
[corr18(55,1),corr18(55,2)] = corr(ERSarous(:),participation(:,5),'rows','complete');
[corr18(56,1),corr18(56,2)] = corr(ERSarous(:),participation(:,6),'rows','complete');
[corr18(57,1),corr18(57,2)] = corr(ERSarous(:),diversity(:,2),'rows','complete');
[corr18(58,1),corr18(58,2)] = corr(ERSarous(:),diversity(:,3),'rows','complete');
[corr18(59,1),corr18(59,2)] = corr(ERSarous(:),diversity(:,5),'rows','complete');
[corr18(60,1),corr18(60,2)] = corr(ERSarous(:),diversity(:,6),'rows','complete');

%% correlations between ERS persistence and network metrics
[corr18(61,1),corr18(61,2)] = corr(ERSpers(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(62,1),corr18(62,2)] = corr(ERSpers(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(63,1),corr18(63,2)] = corr(ERSpers(:),networkDensity(:,3),'rows','complete');
[corr18(64,1),corr18(64,2)] = corr(ERSpers(:),networkDensity(:,4),'rows','complete');
[corr18(65,1),corr18(65,2)] = corr(ERSpers(:),participation(:,2),'rows','complete');
[corr18(66,1),corr18(66,2)] = corr(ERSpers(:),participation(:,3),'rows','complete');
[corr18(67,1),corr18(67,2)] = corr(ERSpers(:),participation(:,5),'rows','complete');
[corr18(68,1),corr18(68,2)] = corr(ERSpers(:),participation(:,6),'rows','complete');
[corr18(69,1),corr18(69,2)] = corr(ERSpers(:),diversity(:,2),'rows','complete');
[corr18(70,1),corr18(70,2)] = corr(ERSpers(:),diversity(:,3),'rows','complete');
[corr18(71,1),corr18(71,2)] = corr(ERSpers(:),diversity(:,5),'rows','complete');
[corr18(72,1),corr18(72,2)] = corr(ERSpers(:),diversity(:,6),'rows','complete');

%% correlations between TAS DDF and network metrics
[corr18(73,1),corr18(73,2)] = corr(TASddf(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(74,1),corr18(74,2)] = corr(TASddf(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(75,1),corr18(75,2)] = corr(TASddf(:),networkDensity(:,3),'rows','complete');
[corr18(76,1),corr18(76,2)] = corr(TASddf(:),networkDensity(:,4),'rows','complete');
[corr18(77,1),corr18(77,2)] = corr(TASddf(:),participation(:,2),'rows','complete');
[corr18(78,1),corr18(78,2)] = corr(TASddf(:),participation(:,3),'rows','complete');
[corr18(79,1),corr18(79,2)] = corr(TASddf(:),participation(:,5),'rows','complete');
[corr18(80,1),corr18(80,2)] = corr(TASddf(:),participation(:,6),'rows','complete');
[corr18(81,1),corr18(81,2)] = corr(TASddf(:),diversity(:,2),'rows','complete');
[corr18(82,1),corr18(82,2)] = corr(TASddf(:),diversity(:,3),'rows','complete');
[corr18(83,1),corr18(83,2)] = corr(TASddf(:),diversity(:,5),'rows','complete');
[corr18(84,1),corr18(84,2)] = corr(TASddf(:),diversity(:,6),'rows','complete');

%% correlations between TAS DIF and network metrics
[corr18(85,1),corr18(85,2)] = corr(TASdif(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(86,1),corr18(86,2)] = corr(TASdif(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(87,1),corr18(87,2)] = corr(TASdif(:),networkDensity(:,3),'rows','complete');
[corr18(88,1),corr18(88,2)] = corr(TASdif(:),networkDensity(:,4),'rows','complete');
[corr18(89,1),corr18(89,2)] = corr(TASdif(:),participation(:,2),'rows','complete');
[corr18(90,1),corr18(90,2)] = corr(TASdif(:),participation(:,3),'rows','complete');
[corr18(91,1),corr18(91,2)] = corr(TASdif(:),participation(:,5),'rows','complete');
[corr18(92,1),corr18(92,2)] = corr(TASdif(:),participation(:,6),'rows','complete');
[corr18(93,1),corr18(93,2)] = corr(TASdif(:),diversity(:,2),'rows','complete');
[corr18(94,1),corr18(94,2)] = corr(TASdif(:),diversity(:,3),'rows','complete');
[corr18(95,1),corr18(95,2)] = corr(TASdif(:),diversity(:,5),'rows','complete');
[corr18(96,1),corr18(96,2)] = corr(TASdif(:),diversity(:,6),'rows','complete');

%% organize and save output
correlations = array2table(corr18);
correlations = horzcat(checks,correlations);
correlations.Properties.VariableNames = {'v1','v2','gamma','threshold','r','p'};

if print == 1
    filename = ['dataSet' num2str(dataSet) '_criterion_validity_checks_copy.xlsx'];
    writetable(correlations,filename);
end