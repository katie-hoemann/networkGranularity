% This script can be used to check construct validity for network metrics
% against emodiversity

%% Instructions:
% run the following scripts:
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
dataSet = 39; % set to 18, 39, or 88 and the rest will take care of itself!
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
checks = readtable('validityChecks_emodiversity.csv');
if dataSet == 18
    rawData = importdata('18ARIVariations_QuestionnaireData.xlsx');
    TAStotal = rawData.data(:,17); %39 for session 2
    TASddf = rawData.data(:,15); %37 for session 2
    TASdif = rawData.data(:,16); %38 for session 2
    TASeot = rawData.data(:,14); %36 for session 2
    range = rawData.data(:,19); %44 for session 2
    diff = rawData.data(:,18); %43 for session 2
    RDEES = rawData.data(:,20); %45 for session 2
elseif dataSet == 39
    rawData = importdata('39ESdata_QuestionnaireData.xlsx');
    TAStotal = rawData.data(:,2);
    TASddf = rawData.data(:,3);
    TASdif = rawData.data(:,4);
    TASeot = rawData.data(:,5);
    range = rawData.data(:,8);
    diff = rawData.data(:,7);
    RDEES = rawData.data(:,6);
end 

%% specify subjects
if dataSet == 88
    clusterCoefGlobal = clusterCoefGlobal(1:length(rtoZ),:);
    diversity = diversity(1:length(rtoZ),:);
    gateway = gateway(1:length(rtoZ),:);
    networkDensity = networkDensity(1:length(rtoZ),:);
    numCommunities = numCommunities(1:length(rtoZ),:);
    participation = participation(1:length(rtoZ),:);
end
    
%% correlations between emodiversity and number of communities at thresholds
[emoDivCorr(1,1),emoDivCorr(1,2)] = corr(emodiversity(:),numCommunities(:,2),'rows','complete');
[emoDivCorr(2,1),emoDivCorr(2,2)] = corr(emodiversity(:),numCommunities(:,3),'rows','complete');
[emoDivCorr(3,1),emoDivCorr(3,2)] = corr(emodiversity(:),numCommunities(:,4),'rows','complete');

%% correlations between emodiversity and global clustering coefficients at thresholds
[emoDivCorr(4,1),emoDivCorr(4,2)] = corr(emodiversity(:),clusterCoefGlobal(:,2),'rows','complete');
[emoDivCorr(5,1),emoDivCorr(5,2)] = corr(emodiversity(:),clusterCoefGlobal(:,3),'rows','complete');
[emoDivCorr(6,1),emoDivCorr(6,2)] = corr(emodiversity(:),clusterCoefGlobal(:,4),'rows','complete');

%% correlations between emodiversity and network density at thresholds
[emoDivCorr(7,1),emoDivCorr(7,2)] = corr(emodiversity(:),networkDensity(:,2),'rows','complete');
[emoDivCorr(8,1),emoDivCorr(8,2)] = corr(emodiversity(:),networkDensity(:,3),'rows','complete');
[emoDivCorr(9,1),emoDivCorr(9,2)] = corr(emodiversity(:),networkDensity(:,4),'rows','complete');

%% correlations between emodiversity and participation coefficients at thresholds
[emoDivCorr(10,1),emoDivCorr(10,2)] = corr(emodiversity(:),participation(:,2),'rows','complete');
[emoDivCorr(11,1),emoDivCorr(11,2)] = corr(emodiversity(:),participation(:,3),'rows','complete');
[emoDivCorr(12,1),emoDivCorr(12,2)] = corr(emodiversity(:),participation(:,4),'rows','complete');

[emoDivCorr(13,1),emoDivCorr(13,2)] = corr(emodiversity(:),participation(:,5),'rows','complete');
[emoDivCorr(14,1),emoDivCorr(14,2)] = corr(emodiversity(:),participation(:,6),'rows','complete');
[emoDivCorr(15,1),emoDivCorr(15,2)] = corr(emodiversity(:),participation(:,7),'rows','complete');

%% correlations between emodiversity and gateway coefficients at thresholds
[emoDivCorr(16,1),emoDivCorr(16,2)] = corr(emodiversity(:),gateway(:,2),'rows','complete'); % based on node strength
[emoDivCorr(17,1),emoDivCorr(17,2)] = corr(emodiversity(:),gateway(:,3),'rows','complete'); % based on node strength
[emoDivCorr(18,1),emoDivCorr(18,2)] = corr(emodiversity(:),gateway(:,4),'rows','complete'); % based on node strength

[emoDivCorr(19,1),emoDivCorr(19,2)] = corr(emodiversity(:),gateway(:,5),'rows','complete'); % based on node strength
[emoDivCorr(20,1),emoDivCorr(20,2)] = corr(emodiversity(:),gateway(:,6),'rows','complete'); % based on node strength
[emoDivCorr(21,1),emoDivCorr(21,2)] = corr(emodiversity(:),gateway(:,7),'rows','complete'); % based on node strength

[emoDivCorr(22,1),emoDivCorr(22,2)] = corr(emodiversity(:),gateway(:,8),'rows','complete'); % based on betweenness centrality
[emoDivCorr(23,1),emoDivCorr(23,2)] = corr(emodiversity(:),gateway(:,9),'rows','complete'); % based on betweenness centrality
[emoDivCorr(24,1),emoDivCorr(24,2)] = corr(emodiversity(:),gateway(:,10),'rows','complete'); % based on betweenness centrality

[emoDivCorr(25,1),emoDivCorr(25,2)] = corr(emodiversity(:),gateway(:,11),'rows','complete'); % based on betweenness centrality
[emoDivCorr(26,1),emoDivCorr(26,2)] = corr(emodiversity(:),gateway(:,12),'rows','complete'); % based on betweenness centrality
[emoDivCorr(27,1),emoDivCorr(27,2)] = corr(emodiversity(:),gateway(:,13),'rows','complete'); % based on betweenness centrality

%% correlations between emodiversity and diversity coefficients at thresholds
[emoDivCorr(28,1),emoDivCorr(28,2)] = corr(emodiversity(:),diversity(:,2),'rows','complete');
[emoDivCorr(29,1),emoDivCorr(29,2)] = corr(emodiversity(:),diversity(:,3),'rows','complete');
[emoDivCorr(30,1),emoDivCorr(30,2)] = corr(emodiversity(:),diversity(:,4),'rows','complete');

[emoDivCorr(31,1),emoDivCorr(31,2)] = corr(emodiversity(:),diversity(:,5),'rows','complete');
[emoDivCorr(32,1),emoDivCorr(32,2)] = corr(emodiversity(:),diversity(:,6),'rows','complete');
[emoDivCorr(33,1),emoDivCorr(33,2)] = corr(emodiversity(:),diversity(:,7),'rows','complete');

%% correlation between emodiversity and zICC, alexithymia, range, and differentiation
[emoDivCorr(34,1),emoDivCorr(34,2)] = corr(emodiversity(:),rtoZ(:,3),'rows','complete');
[emoDivCorr(35,1),emoDivCorr(35,2)] = corr(emodiversity(:),TAStotal(:),'rows','complete');
[emoDivCorr(36,1),emoDivCorr(36,2)] = corr(emodiversity(:),TASddf(:),'rows','complete');
[emoDivCorr(37,1),emoDivCorr(37,2)] = corr(emodiversity(:),TASdif(:),'rows','complete');
[emoDivCorr(38,1),emoDivCorr(38,2)] = corr(emodiversity(:),TASeot(:),'rows','complete');
[emoDivCorr(39,1),emoDivCorr(39,2)] = corr(emodiversity(:),range(:),'rows','complete');
[emoDivCorr(40,1),emoDivCorr(40,2)] = corr(emodiversity(:),diff(:),'rows','complete');
[emoDivCorr(41,1),emoDivCorr(41,2)] = corr(emodiversity(:),RDEES(:),'rows','complete');

%% organize and save output
correlations = array2table(emoDivCorr);
correlations = horzcat(checks,correlations);
correlations.Properties.VariableNames = {'CheckName','v1','v2','gamma','threshold','r','p'};

if print == 1
    filename = ['dataSet' num2str(dataSet) '_construct_validity_checks_emodiversity.xlsx'];
    writetable(correlations,filename);
end

%% still to do
% could calculate separate emodiversity for pos and neg
% could grab SDs for participation and clustering coefs
% gini coefficient has not been integrated