% This script can be used to check construct validity for network metrics
% against TAS-20 scores

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
checks = readtable('validityChecks_alexithymia.csv');
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
    
%% correlations between TAS total and number of communities at thresholds
[TAScorr(1,1),TAScorr(1,2)] = corr(TAStotal(:),numCommunities(:,2),'rows','complete');
[TAScorr(2,1),TAScorr(2,2)] = corr(TAStotal(:),numCommunities(:,3),'rows','complete');
[TAScorr(3,1),TAScorr(3,2)] = corr(TAStotal(:),numCommunities(:,4),'rows','complete');

%% correlations between TAS total and global clustering coefficients at thresholds
[TAScorr(4,1),TAScorr(4,2)] = corr(TAStotal(:),clusterCoefGlobal(:,2),'rows','complete');
[TAScorr(5,1),TAScorr(5,2)] = corr(TAStotal(:),clusterCoefGlobal(:,3),'rows','complete');
[TAScorr(6,1),TAScorr(6,2)] = corr(TAStotal(:),clusterCoefGlobal(:,4),'rows','complete');

%% correlations between TAS total and network density at thresholds
[TAScorr(7,1),TAScorr(7,2)] = corr(TAStotal(:),networkDensity(:,2),'rows','complete');
[TAScorr(8,1),TAScorr(8,2)] = corr(TAStotal(:),networkDensity(:,3),'rows','complete');
[TAScorr(9,1),TAScorr(9,2)] = corr(TAStotal(:),networkDensity(:,4),'rows','complete');

%% correlations between TAS total and participation coefficients at thresholds
[TAScorr(10,1),TAScorr(10,2)] = corr(TAStotal(:),participation(:,2),'rows','complete');
[TAScorr(11,1),TAScorr(11,2)] = corr(TAStotal(:),participation(:,3),'rows','complete');
[TAScorr(12,1),TAScorr(12,2)] = corr(TAStotal(:),participation(:,4),'rows','complete');

[TAScorr(13,1),TAScorr(13,2)] = corr(TAStotal(:),participation(:,5),'rows','complete');
[TAScorr(14,1),TAScorr(14,2)] = corr(TAStotal(:),participation(:,6),'rows','complete');
[TAScorr(15,1),TAScorr(15,2)] = corr(TAStotal(:),participation(:,7),'rows','complete');

%% correlations between TAS total and gateway coefficients at thresholds
[TAScorr(16,1),TAScorr(16,2)] = corr(TAStotal(:),gateway(:,2),'rows','complete'); % based on node strength
[TAScorr(17,1),TAScorr(17,2)] = corr(TAStotal(:),gateway(:,3),'rows','complete'); % based on node strength
[TAScorr(18,1),TAScorr(18,2)] = corr(TAStotal(:),gateway(:,4),'rows','complete'); % based on node strength

[TAScorr(19,1),TAScorr(19,2)] = corr(TAStotal(:),gateway(:,5),'rows','complete'); % based on node strength
[TAScorr(20,1),TAScorr(20,2)] = corr(TAStotal(:),gateway(:,6),'rows','complete'); % based on node strength
[TAScorr(21,1),TAScorr(21,2)] = corr(TAStotal(:),gateway(:,7),'rows','complete'); % based on node strength

[TAScorr(22,1),TAScorr(22,2)] = corr(TAStotal(:),gateway(:,8),'rows','complete'); % based on betweenness centrality
[TAScorr(23,1),TAScorr(23,2)] = corr(TAStotal(:),gateway(:,9),'rows','complete'); % based on betweenness centrality
[TAScorr(24,1),TAScorr(24,2)] = corr(TAStotal(:),gateway(:,10),'rows','complete'); % based on betweenness centrality

[TAScorr(25,1),TAScorr(25,2)] = corr(TAStotal(:),gateway(:,11),'rows','complete'); % based on betweenness centrality
[TAScorr(26,1),TAScorr(26,2)] = corr(TAStotal(:),gateway(:,12),'rows','complete'); % based on betweenness centrality
[TAScorr(27,1),TAScorr(27,2)] = corr(TAStotal(:),gateway(:,13),'rows','complete'); % based on betweenness centrality

%% correlations between TAS total and diversity coefficients at thresholds
[TAScorr(28,1),TAScorr(28,2)] = corr(TAStotal(:),diversity(:,2),'rows','complete');
[TAScorr(29,1),TAScorr(29,2)] = corr(TAStotal(:),diversity(:,3),'rows','complete');
[TAScorr(30,1),TAScorr(30,2)] = corr(TAStotal(:),diversity(:,4),'rows','complete');

[TAScorr(31,1),TAScorr(31,2)] = corr(TAStotal(:),diversity(:,5),'rows','complete');
[TAScorr(32,1),TAScorr(32,2)] = corr(TAStotal(:),diversity(:,6),'rows','complete');
[TAScorr(33,1),TAScorr(33,2)] = corr(TAStotal(:),diversity(:,7),'rows','complete');

%% correlation between TAS total and zICC, emodiversity, range, and differentiation
[TAScorr(34,1),TAScorr(34,2)] = corr(TAStotal(:),rtoZ(:,3),'rows','complete');
[TAScorr(35,1),TAScorr(35,2)] = corr(TAStotal(:),emodiversity(:),'rows','complete');
[TAScorr(36,1),TAScorr(36,2)] = corr(TAStotal(:),TASddf(:),'rows','complete');
[TAScorr(37,1),TAScorr(37,2)] = corr(TAStotal(:),TASdif(:),'rows','complete');
[TAScorr(38,1),TAScorr(38,2)] = corr(TAStotal(:),TASeot(:),'rows','complete');
[TAScorr(39,1),TAScorr(39,2)] = corr(TAStotal(:),range(:),'rows','complete');
[TAScorr(40,1),TAScorr(40,2)] = corr(TAStotal(:),diff(:),'rows','complete');
[TAScorr(41,1),TAScorr(41,2)] = corr(TAStotal(:),RDEES(:),'rows','complete');

%% organize and save output
correlations = array2table(TAScorr);
correlations = horzcat(checks,correlations);
correlations.Properties.VariableNames = {'CheckName','v1','v2','gamma','threshold','r','p'};

if print == 1
    filename = ['dataSet' num2str(dataSet) '_construct_validity_checks_alexithymia.xlsx'];
    writetable(correlations,filename);
end

%% still to do
% could grab SDs for participation and clustering coefs
% gini coefficient has not been integrated