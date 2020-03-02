% This script can be used to check criterion validity for network metrics
% on the 39-term dataset

%% Instructions:
% run the following scripts for the 39-term dataset:
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
dataSet = 39; % must be set to 39
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
checks = readtable('criterionChecks_39_subscales.csv');

%% load criterion validity measures
rawData = importdata('39ESdata_QuestionnaireData.xlsx');
BAIm = rawData.data(:,11); 
BAIb = rawData.data(:,10);
TASddf = rawData.data(:,3);
TASdif = rawData.data(:,4);

%% remove subject 31 (row 29) for all calculations
rtoZ(29,:) = [];
numCommunities(29,:) = [];
clusterCoefGlobal(29,:) = [];
participation(29,:) = [];
gateway(29,:) = [];
diversity(29,:) = [];
networkDensity(29,:) = [];
commRadius(29,:) = [];
emodiversity(29,:) = [];
BAIm(29,:) = [];
BAIb(29,:) = [];
TASddf(29,:) = [];
TASdif(29,:) = [];
   
%% correlations between BAI m subscale and network metrics
[corr39(1,1),corr39(1,2)] = corr(BAIm(:),clusterCoefGlobal(:,3),'rows','complete');
[corr39(2,1),corr39(2,2)] = corr(BAIm(:),clusterCoefGlobal(:,4),'rows','complete');
[corr39(3,1),corr39(3,2)] = corr(BAIm(:),networkDensity(:,3),'rows','complete');
[corr39(4,1),corr39(4,2)] = corr(BAIm(:),networkDensity(:,4),'rows','complete');
[corr39(5,1),corr39(5,2)] = corr(BAIm(:),participation(:,2),'rows','complete');
[corr39(6,1),corr39(6,2)] = corr(BAIm(:),participation(:,3),'rows','complete');
[corr39(7,1),corr39(7,2)] = corr(BAIm(:),participation(:,5),'rows','complete');
[corr39(8,1),corr39(8,2)] = corr(BAIm(:),participation(:,6),'rows','complete');
[corr39(9,1),corr39(9,2)] = corr(BAIm(:),diversity(:,2),'rows','complete');
[corr39(10,1),corr39(10,2)] = corr(BAIm(:),diversity(:,3),'rows','complete');
[corr39(11,1),corr39(11,2)] = corr(BAIm(:),diversity(:,5),'rows','complete');
[corr39(12,1),corr39(12,2)] = corr(BAIm(:),diversity(:,6),'rows','complete');

%% correlations between BAI b subscale and network metrics
[corr39(13,1),corr39(13,2)] = corr(BAIb(:),clusterCoefGlobal(:,3),'rows','complete');
[corr39(14,1),corr39(14,2)] = corr(BAIb(:),clusterCoefGlobal(:,4),'rows','complete');
[corr39(15,1),corr39(15,2)] = corr(BAIb(:),networkDensity(:,3),'rows','complete');
[corr39(16,1),corr39(16,2)] = corr(BAIb(:),networkDensity(:,4),'rows','complete');
[corr39(17,1),corr39(17,2)] = corr(BAIb(:),participation(:,2),'rows','complete');
[corr39(18,1),corr39(18,2)] = corr(BAIb(:),participation(:,3),'rows','complete');
[corr39(19,1),corr39(19,2)] = corr(BAIb(:),participation(:,5),'rows','complete');
[corr39(20,1),corr39(20,2)] = corr(BAIb(:),participation(:,6),'rows','complete');
[corr39(21,1),corr39(21,2)] = corr(BAIb(:),diversity(:,2),'rows','complete');
[corr39(22,1),corr39(22,2)] = corr(BAIb(:),diversity(:,3),'rows','complete');
[corr39(23,1),corr39(23,2)] = corr(BAIb(:),diversity(:,5),'rows','complete');
[corr39(24,1),corr39(24,2)] = corr(BAIb(:),diversity(:,6),'rows','complete');

%% correlations between TAS DDF and network metrics
[corr39(25,1),corr39(25,2)] = corr(TASddf(:),clusterCoefGlobal(:,3),'rows','complete');
[corr39(26,1),corr39(26,2)] = corr(TASddf(:),clusterCoefGlobal(:,4),'rows','complete');
[corr39(27,1),corr39(27,2)] = corr(TASddf(:),networkDensity(:,3),'rows','complete');
[corr39(28,1),corr39(28,2)] = corr(TASddf(:),networkDensity(:,4),'rows','complete');
[corr39(29,1),corr39(29,2)] = corr(TASddf(:),participation(:,2),'rows','complete');
[corr39(30,1),corr39(30,2)] = corr(TASddf(:),participation(:,3),'rows','complete');
[corr39(31,1),corr39(31,2)] = corr(TASddf(:),participation(:,5),'rows','complete');
[corr39(32,1),corr39(32,2)] = corr(TASddf(:),participation(:,6),'rows','complete');
[corr39(33,1),corr39(33,2)] = corr(TASddf(:),diversity(:,2),'rows','complete');
[corr39(34,1),corr39(34,2)] = corr(TASddf(:),diversity(:,3),'rows','complete');
[corr39(35,1),corr39(35,2)] = corr(TASddf(:),diversity(:,5),'rows','complete');
[corr39(36,1),corr39(36,2)] = corr(TASddf(:),diversity(:,6),'rows','complete');

%% correlations between TAS DIF and network metrics
[corr39(37,1),corr39(37,2)] = corr(TASdif(:),clusterCoefGlobal(:,3),'rows','complete');
[corr39(38,1),corr39(38,2)] = corr(TASdif(:),clusterCoefGlobal(:,4),'rows','complete');
[corr39(39,1),corr39(39,2)] = corr(TASdif(:),networkDensity(:,3),'rows','complete');
[corr39(40,1),corr39(40,2)] = corr(TASdif(:),networkDensity(:,4),'rows','complete');
[corr39(41,1),corr39(41,2)] = corr(TASdif(:),participation(:,2),'rows','complete');
[corr39(42,1),corr39(42,2)] = corr(TASdif(:),participation(:,3),'rows','complete');
[corr39(43,1),corr39(43,2)] = corr(TASdif(:),participation(:,5),'rows','complete');
[corr39(44,1),corr39(44,2)] = corr(TASdif(:),participation(:,6),'rows','complete');
[corr39(45,1),corr39(45,2)] = corr(TASdif(:),diversity(:,2),'rows','complete');
[corr39(46,1),corr39(46,2)] = corr(TASdif(:),diversity(:,3),'rows','complete');
[corr39(47,1),corr39(47,2)] = corr(TASdif(:),diversity(:,5),'rows','complete');
[corr39(48,1),corr39(48,2)] = corr(TASdif(:),diversity(:,6),'rows','complete');

%% organize and save output
correlations = array2table(corr39);
correlations = horzcat(checks,correlations);
correlations.Properties.VariableNames = {'v1','v2','gamma','threshold','r','p'};

if print == 1
    filename = ['dataSet' num2str(dataSet) '_criterion_validity_checks_copy.xlsx'];
    writetable(correlations,filename);
end