% This script can be used to check criterion validity for network metrics
% on the 39-term dataset

%% Instructions:
% run the following scripts for the 39-term dataset:
% - assign_subject_communities_loop (set subjectFiles to 0)
% - check_subject_communities_loop (set subjectFiles to 0)
% - check_subject_clustering_loop (set subjectFiles to 0)
% - check_subject_participation_loop (set subjectFiles to 0)
% - check_subject_density_loop
% - granularity_calculations
% - emodiversity_calculation
% - affect_calculation
% make sure the necessary .mat files are on your path

clear;
clc;

%% specify dataset
dataSet = 39; % must be set to 39
print = 1; % set to 1 to generate spreadsheet with results

%% load data and list of checks
filenames{1} = ['dataSet' num2str(dataSet) '_zICC.mat'];
filenames{2} = ['dataSet' num2str(dataSet) '_communities.mat'];
filenames{3} = ['dataSet' num2str(dataSet) '_modularity.mat'];
filenames{4} = ['dataSet' num2str(dataSet) '_modularity_VA.mat'];
filenames{5} = ['dataSet' num2str(dataSet) '_clustering.mat'];
filenames{6} = ['dataSet' num2str(dataSet) '_participation.mat'];
filenames{7} = ['dataSet' num2str(dataSet) '_gateway.mat'];
filenames{8} = ['dataSet' num2str(dataSet) '_diversity.mat'];
filenames{9} = ['dataSet' num2str(dataSet) '_density.mat'];
filenames{10} = ['dataSet' num2str(dataSet) '_comm_radius.mat'];
filenames{11} = ['dataSet' num2str(dataSet) '_emodiversity.mat'];
filenames{12} = ['dataSet' num2str(dataSet) '_affect.mat'];
for i_file = 1:numel(filenames)
    load(filenames{i_file});
end
checks = readtable('criterionChecks_39.csv');

%% load criterion validity measures
rawData = importdata('39ESdata_QuestionnaireData.xlsx');
BAItotal = rawData.data(:,9); 
BDItotal = rawData.data(:,12);
SWLStotal = rawData.data(:,13);

%% remove subject 31 (row 29) for all calculations
rtoZ(29,:) = [];
numCommunities(29,:) = [];
modularity(29,:) = [];
modularity_VA(29,:) = [];
clusterCoefGlobal(29,:) = [];
participation(29,:) = [];
gateway(29,:) = [];
diversity(29,:) = [];
networkDensity(29,:) = [];
commRadius(29,:) = [];
emodiversity(29,:) = [];
affect(29,:) = [];
BAItotal(29,:) = [];
BDItotal(29,:) = [];
SWLStotal(29,:) = [];
   
%% correlations between BAI and network metrics
[corr39(1,1),corr39(1,2)] = corr(BAItotal(:),numCommunities(:,2),'rows','complete');
[corr39(2,1),corr39(2,2)] = corr(BAItotal(:),numCommunities(:,3),'rows','complete');
[corr39(3,1),corr39(3,2)] = corr(BAItotal(:),numCommunities(:,4),'rows','complete');
[corr39(61,1),corr39(61,2)] = corr(BAItotal(:),modularity(:,2),'rows','complete');
[corr39(62,1),corr39(62,2)] = corr(BAItotal(:),modularity(:,3),'rows','complete');
[corr39(63,1),corr39(63,2)] = corr(BAItotal(:),modularity(:,4),'rows','complete');
[corr39(64,1),corr39(64,2)] = corr(BAItotal(:),modularity_VA(:,2),'rows','complete');
[corr39(65,1),corr39(65,2)] = corr(BAItotal(:),modularity_VA(:,3),'rows','complete');
[corr39(66,1),corr39(66,2)] = corr(BAItotal(:),modularity_VA(:,4),'rows','complete');
[corr39(4,1),corr39(4,2)] = corr(BAItotal(:),clusterCoefGlobal(:,2),'rows','complete');
[corr39(5,1),corr39(5,2)] = corr(BAItotal(:),clusterCoefGlobal(:,3),'rows','complete');
[corr39(6,1),corr39(6,2)] = corr(BAItotal(:),clusterCoefGlobal(:,4),'rows','complete');
[corr39(7,1),corr39(7,2)] = corr(BAItotal(:),networkDensity(:,3),'rows','complete');
[corr39(8,1),corr39(8,2)] = corr(BAItotal(:),networkDensity(:,4),'rows','complete');
[corr39(67,1),corr39(67,2)] = corr(BAItotal(:),networkDensity(:,5),'rows','complete');
[corr39(68,1),corr39(68,2)] = corr(BAItotal(:),networkDensity(:,6),'rows','complete');
[corr39(69,1),corr39(69,2)] = corr(BAItotal(:),networkDensity(:,7),'rows','complete');
[corr39(9,1),corr39(9,2)] = corr(BAItotal(:),participation(:,2),'rows','complete');
[corr39(10,1),corr39(10,2)] = corr(BAItotal(:),participation(:,3),'rows','complete');
[corr39(11,1),corr39(11,2)] = corr(BAItotal(:),participation(:,5),'rows','complete');
[corr39(12,1),corr39(12,2)] = corr(BAItotal(:),participation(:,6),'rows','complete');
[corr39(13,1),corr39(13,2)] = corr(BAItotal(:),gateway(:,2),'rows','complete'); % node strength
[corr39(14,1),corr39(14,2)] = corr(BAItotal(:),gateway(:,3),'rows','complete'); % node strength
[corr39(15,1),corr39(15,2)] = corr(BAItotal(:),gateway(:,5),'rows','complete'); % node strength
[corr39(16,1),corr39(16,2)] = corr(BAItotal(:),gateway(:,6),'rows','complete'); % node strength
[corr39(17,1),corr39(17,2)] = corr(BAItotal(:),diversity(:,2),'rows','complete');
[corr39(18,1),corr39(18,2)] = corr(BAItotal(:),diversity(:,3),'rows','complete');
[corr39(19,1),corr39(19,2)] = corr(BAItotal(:),diversity(:,5),'rows','complete');
[corr39(20,1),corr39(20,2)] = corr(BAItotal(:),diversity(:,6),'rows','complete');
[corr39(103,1),corr39(103,2)] = corr(BAItotal(:),commRadius(:,2),'rows','complete');
[corr39(104,1),corr39(104,2)] = corr(BAItotal(:),commRadius(:,3),'rows','complete');
[corr39(105,1),corr39(105,2)] = corr(BAItotal(:),commRadius(:,4),'rows','complete');

%% correlations between BDI and network metrics
[corr39(21,1),corr39(21,2)] = corr(BDItotal(:),numCommunities(:,2),'rows','complete');
[corr39(22,1),corr39(22,2)] = corr(BDItotal(:),numCommunities(:,3),'rows','complete');
[corr39(23,1),corr39(23,2)] = corr(BDItotal(:),numCommunities(:,4),'rows','complete');
[corr39(70,1),corr39(70,2)] = corr(BDItotal(:),modularity(:,2),'rows','complete');
[corr39(71,1),corr39(71,2)] = corr(BDItotal(:),modularity(:,3),'rows','complete');
[corr39(72,1),corr39(72,2)] = corr(BDItotal(:),modularity(:,4),'rows','complete');
[corr39(73,1),corr39(73,2)] = corr(BDItotal(:),modularity_VA(:,2),'rows','complete');
[corr39(74,1),corr39(74,2)] = corr(BDItotal(:),modularity_VA(:,3),'rows','complete');
[corr39(75,1),corr39(75,2)] = corr(BDItotal(:),modularity_VA(:,4),'rows','complete');
[corr39(24,1),corr39(24,2)] = corr(BDItotal(:),clusterCoefGlobal(:,2),'rows','complete');
[corr39(25,1),corr39(25,2)] = corr(BDItotal(:),clusterCoefGlobal(:,3),'rows','complete');
[corr39(26,1),corr39(26,2)] = corr(BDItotal(:),clusterCoefGlobal(:,4),'rows','complete');
[corr39(27,1),corr39(27,2)] = corr(BDItotal(:),networkDensity(:,3),'rows','complete');
[corr39(28,1),corr39(28,2)] = corr(BDItotal(:),networkDensity(:,4),'rows','complete');
[corr39(76,1),corr39(76,2)] = corr(BDItotal(:),networkDensity(:,5),'rows','complete');
[corr39(77,1),corr39(77,2)] = corr(BDItotal(:),networkDensity(:,6),'rows','complete');
[corr39(78,1),corr39(78,2)] = corr(BDItotal(:),networkDensity(:,7),'rows','complete');
[corr39(29,1),corr39(29,2)] = corr(BDItotal(:),participation(:,2),'rows','complete');
[corr39(30,1),corr39(30,2)] = corr(BDItotal(:),participation(:,3),'rows','complete');
[corr39(31,1),corr39(31,2)] = corr(BDItotal(:),participation(:,5),'rows','complete');
[corr39(32,1),corr39(32,2)] = corr(BDItotal(:),participation(:,6),'rows','complete');
[corr39(33,1),corr39(33,2)] = corr(BDItotal(:),gateway(:,2),'rows','complete'); % node strength
[corr39(34,1),corr39(34,2)] = corr(BDItotal(:),gateway(:,3),'rows','complete'); % node strength
[corr39(35,1),corr39(35,2)] = corr(BDItotal(:),gateway(:,5),'rows','complete'); % node strength
[corr39(36,1),corr39(36,2)] = corr(BDItotal(:),gateway(:,6),'rows','complete'); % node strength
[corr39(37,1),corr39(37,2)] = corr(BDItotal(:),diversity(:,2),'rows','complete');
[corr39(38,1),corr39(38,2)] = corr(BDItotal(:),diversity(:,3),'rows','complete');
[corr39(39,1),corr39(39,2)] = corr(BDItotal(:),diversity(:,5),'rows','complete');
[corr39(40,1),corr39(40,2)] = corr(BDItotal(:),diversity(:,6),'rows','complete');
[corr39(106,1),corr39(106,2)] = corr(BDItotal(:),commRadius(:,2),'rows','complete');
[corr39(107,1),corr39(107,2)] = corr(BDItotal(:),commRadius(:,3),'rows','complete');
[corr39(108,1),corr39(108,2)] = corr(BDItotal(:),commRadius(:,4),'rows','complete');

%% correlations between SWLS and network metrics
[corr39(41,1),corr39(41,2)] = corr(SWLStotal(:),numCommunities(:,2),'rows','complete');
[corr39(42,1),corr39(42,2)] = corr(SWLStotal(:),numCommunities(:,3),'rows','complete');
[corr39(43,1),corr39(43,2)] = corr(SWLStotal(:),numCommunities(:,4),'rows','complete');
[corr39(79,1),corr39(79,2)] = corr(SWLStotal(:),modularity(:,2),'rows','complete');
[corr39(80,1),corr39(80,2)] = corr(SWLStotal(:),modularity(:,3),'rows','complete');
[corr39(81,1),corr39(81,2)] = corr(SWLStotal(:),modularity(:,4),'rows','complete');
[corr39(82,1),corr39(82,2)] = corr(SWLStotal(:),modularity_VA(:,2),'rows','complete');
[corr39(83,1),corr39(83,2)] = corr(SWLStotal(:),modularity_VA(:,3),'rows','complete');
[corr39(84,1),corr39(84,2)] = corr(SWLStotal(:),modularity_VA(:,4),'rows','complete');
[corr39(44,1),corr39(44,2)] = corr(SWLStotal(:),clusterCoefGlobal(:,2),'rows','complete');
[corr39(45,1),corr39(45,2)] = corr(SWLStotal(:),clusterCoefGlobal(:,3),'rows','complete');
[corr39(46,1),corr39(46,2)] = corr(SWLStotal(:),clusterCoefGlobal(:,4),'rows','complete');
[corr39(47,1),corr39(47,2)] = corr(SWLStotal(:),networkDensity(:,3),'rows','complete');
[corr39(48,1),corr39(48,2)] = corr(SWLStotal(:),networkDensity(:,4),'rows','complete');
[corr39(85,1),corr39(85,2)] = corr(SWLStotal(:),networkDensity(:,5),'rows','complete');
[corr39(86,1),corr39(86,2)] = corr(SWLStotal(:),networkDensity(:,6),'rows','complete');
[corr39(87,1),corr39(87,2)] = corr(SWLStotal(:),networkDensity(:,7),'rows','complete');
[corr39(49,1),corr39(49,2)] = corr(SWLStotal(:),participation(:,2),'rows','complete');
[corr39(50,1),corr39(50,2)] = corr(SWLStotal(:),participation(:,3),'rows','complete');
[corr39(51,1),corr39(51,2)] = corr(SWLStotal(:),participation(:,5),'rows','complete');
[corr39(52,1),corr39(52,2)] = corr(SWLStotal(:),participation(:,6),'rows','complete');
[corr39(53,1),corr39(53,2)] = corr(SWLStotal(:),gateway(:,2),'rows','complete'); % node strength
[corr39(54,1),corr39(54,2)] = corr(SWLStotal(:),gateway(:,3),'rows','complete'); % node strength
[corr39(55,1),corr39(55,2)] = corr(SWLStotal(:),gateway(:,5),'rows','complete'); % node strength
[corr39(56,1),corr39(56,2)] = corr(SWLStotal(:),gateway(:,6),'rows','complete'); % node strength
[corr39(57,1),corr39(57,2)] = corr(SWLStotal(:),diversity(:,2),'rows','complete');
[corr39(58,1),corr39(58,2)] = corr(SWLStotal(:),diversity(:,3),'rows','complete');
[corr39(59,1),corr39(59,2)] = corr(SWLStotal(:),diversity(:,5),'rows','complete');
[corr39(60,1),corr39(60,2)] = corr(SWLStotal(:),diversity(:,6),'rows','complete');
[corr39(109,1),corr39(109,2)] = corr(SWLStotal(:),commRadius(:,2),'rows','complete');
[corr39(110,1),corr39(110,2)] = corr(SWLStotal(:),commRadius(:,3),'rows','complete');
[corr39(111,1),corr39(111,2)] = corr(SWLStotal(:),commRadius(:,4),'rows','complete');

%% correlations between zICC and psychological measures
[corr39(88,1),corr39(88,2)] = corr(BAItotal(:),rtoZ(:,3),'rows','complete');
[corr39(89,1),corr39(89,2)] = corr(BDItotal(:),rtoZ(:,3),'rows','complete');
[corr39(90,1),corr39(90,2)] = corr(SWLStotal(:),rtoZ(:,3),'rows','complete');

%% correlations between affect and psychological measures
[corr39(91,1),corr39(91,2)] = corr(BAItotal(:),affect(:,1),'rows','complete');
[corr39(92,1),corr39(92,2)] = corr(BAItotal(:),affect(:,2),'rows','complete');
[corr39(93,1),corr39(93,2)] = corr(BAItotal(:),affect(:,3),'rows','complete');
[corr39(94,1),corr39(94,2)] = corr(BAItotal(:),affect(:,4),'rows','complete');
[corr39(95,1),corr39(95,2)] = corr(BDItotal(:),affect(:,1),'rows','complete');
[corr39(96,1),corr39(96,2)] = corr(BDItotal(:),affect(:,2),'rows','complete');
[corr39(97,1),corr39(97,2)] = corr(BDItotal(:),affect(:,3),'rows','complete');
[corr39(98,1),corr39(98,2)] = corr(BDItotal(:),affect(:,4),'rows','complete');
[corr39(99,1),corr39(99,2)] = corr(SWLStotal(:),affect(:,1),'rows','complete');
[corr39(100,1),corr39(100,2)] = corr(SWLStotal(:),affect(:,2),'rows','complete');
[corr39(101,1),corr39(101,2)] = corr(SWLStotal(:),affect(:,3),'rows','complete');
[corr39(102,1),corr39(102,2)] = corr(SWLStotal(:),affect(:,4),'rows','complete');

%% organize and save output
correlations = array2table(corr39);
correlations = horzcat(checks,correlations);
correlations.Properties.VariableNames = {'v1','v2','gamma','threshold','r','p'};

if print == 1
    filename = ['dataSet' num2str(dataSet) '_criterion_validity_checks_poly.xlsx'];
    writetable(correlations,filename);
end