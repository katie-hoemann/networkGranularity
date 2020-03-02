% This script can be used to check construct validity for network metrics
% against RDEES scores

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
checks = readtable('validityChecks_complexity.csv');
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
    
%% correlations between RDEES total and number of communities at thresholds
[RDEEScorr(1,1),RDEEScorr(1,2)] = corr(RDEES(:),numCommunities(:,2),'rows','complete');
[RDEEScorr(2,1),RDEEScorr(2,2)] = corr(RDEES(:),numCommunities(:,3),'rows','complete');
[RDEEScorr(3,1),RDEEScorr(3,2)] = corr(RDEES(:),numCommunities(:,4),'rows','complete');

%% correlations between RDEES total and global clustering coefficients at thresholds
[RDEEScorr(4,1),RDEEScorr(4,2)] = corr(RDEES(:),clusterCoefGlobal(:,2),'rows','complete');
[RDEEScorr(5,1),RDEEScorr(5,2)] = corr(RDEES(:),clusterCoefGlobal(:,3),'rows','complete');
[RDEEScorr(6,1),RDEEScorr(6,2)] = corr(RDEES(:),clusterCoefGlobal(:,4),'rows','complete');

%% correlations between RDEES total and network density at thresholds
[RDEEScorr(7,1),RDEEScorr(7,2)] = corr(RDEES(:),networkDensity(:,2),'rows','complete');
[RDEEScorr(8,1),RDEEScorr(8,2)] = corr(RDEES(:),networkDensity(:,3),'rows','complete');
[RDEEScorr(9,1),RDEEScorr(9,2)] = corr(RDEES(:),networkDensity(:,4),'rows','complete');

%% correlations between RDEES total and participation coefficients at thresholds
[RDEEScorr(10,1),RDEEScorr(10,2)] = corr(RDEES(:),participation(:,2),'rows','complete');
[RDEEScorr(11,1),RDEEScorr(11,2)] = corr(RDEES(:),participation(:,3),'rows','complete');
[RDEEScorr(12,1),RDEEScorr(12,2)] = corr(RDEES(:),participation(:,4),'rows','complete');

[RDEEScorr(13,1),RDEEScorr(13,2)] = corr(RDEES(:),participation(:,5),'rows','complete');
[RDEEScorr(14,1),RDEEScorr(14,2)] = corr(RDEES(:),participation(:,6),'rows','complete');
[RDEEScorr(15,1),RDEEScorr(15,2)] = corr(RDEES(:),participation(:,7),'rows','complete');

%% correlations between RDEES total and gateway coefficients at thresholds
[RDEEScorr(16,1),RDEEScorr(16,2)] = corr(RDEES(:),gateway(:,2),'rows','complete'); % based on node strength
[RDEEScorr(17,1),RDEEScorr(17,2)] = corr(RDEES(:),gateway(:,3),'rows','complete'); % based on node strength
[RDEEScorr(18,1),RDEEScorr(18,2)] = corr(RDEES(:),gateway(:,4),'rows','complete'); % based on node strength

[RDEEScorr(19,1),RDEEScorr(19,2)] = corr(RDEES(:),gateway(:,5),'rows','complete'); % based on node strength
[RDEEScorr(20,1),RDEEScorr(20,2)] = corr(RDEES(:),gateway(:,6),'rows','complete'); % based on node strength
[RDEEScorr(21,1),RDEEScorr(21,2)] = corr(RDEES(:),gateway(:,7),'rows','complete'); % based on node strength

[RDEEScorr(22,1),RDEEScorr(22,2)] = corr(RDEES(:),gateway(:,8),'rows','complete'); % based on betweenness centrality
[RDEEScorr(23,1),RDEEScorr(23,2)] = corr(RDEES(:),gateway(:,9),'rows','complete'); % based on betweenness centrality
[RDEEScorr(24,1),RDEEScorr(24,2)] = corr(RDEES(:),gateway(:,10),'rows','complete'); % based on betweenness centrality

[RDEEScorr(25,1),RDEEScorr(25,2)] = corr(RDEES(:),gateway(:,11),'rows','complete'); % based on betweenness centrality
[RDEEScorr(26,1),RDEEScorr(26,2)] = corr(RDEES(:),gateway(:,12),'rows','complete'); % based on betweenness centrality
[RDEEScorr(27,1),RDEEScorr(27,2)] = corr(RDEES(:),gateway(:,13),'rows','complete'); % based on betweenness centrality

%% correlations between RDEES total and diversity coefficients at thresholds
[RDEEScorr(28,1),RDEEScorr(28,2)] = corr(RDEES(:),diversity(:,2),'rows','complete');
[RDEEScorr(29,1),RDEEScorr(29,2)] = corr(RDEES(:),diversity(:,3),'rows','complete');
[RDEEScorr(30,1),RDEEScorr(30,2)] = corr(RDEES(:),diversity(:,4),'rows','complete');

[RDEEScorr(31,1),RDEEScorr(31,2)] = corr(RDEES(:),diversity(:,5),'rows','complete');
[RDEEScorr(32,1),RDEEScorr(32,2)] = corr(RDEES(:),diversity(:,6),'rows','complete');
[RDEEScorr(33,1),RDEEScorr(33,2)] = corr(RDEES(:),diversity(:,7),'rows','complete');

%% correlation between RDEES total and zICC, emodiversity, alexithymia, range, and differentiation
[RDEEScorr(34,1),RDEEScorr(34,2)] = corr(RDEES(:),rtoZ(:,3),'rows','complete');
[RDEEScorr(35,1),RDEEScorr(35,2)] = corr(RDEES(:),emodiversity(:),'rows','complete');
[RDEEScorr(36,1),RDEEScorr(36,2)] = corr(RDEES(:),TAStotal(:),'rows','complete');
[RDEEScorr(37,1),RDEEScorr(37,2)] = corr(RDEES(:),TASddf(:),'rows','complete');
[RDEEScorr(38,1),RDEEScorr(38,2)] = corr(RDEES(:),TASdif(:),'rows','complete');
[RDEEScorr(39,1),RDEEScorr(39,2)] = corr(RDEES(:),TASeot(:),'rows','complete');
[RDEEScorr(40,1),RDEEScorr(40,2)] = corr(RDEES(:),range(:),'rows','complete');
[RDEEScorr(41,1),RDEEScorr(41,2)] = corr(RDEES(:),diff(:),'rows','complete');


%% organize and save output
correlations = array2table(RDEEScorr);
correlations = horzcat(checks,correlations);
correlations.Properties.VariableNames = {'CheckName','v1','v2','gamma','threshold','r','p'};

if print == 1
    filename = ['dataSet' num2str(dataSet) '_construct_validity_checks_complexity.xlsx'];
    writetable(correlations,filename);
end

%% still to do
% could grab SDs for participation and clustering coefs
% gini coefficient has not been integrated