% This script can be used to check construct validity for network metrics
% against the traditional ICC-based measure of granularity

%% Instructions:
% run the following scripts:
% - assign_subject_communities_loop (set subjectFiles to 0)
% - check_subject_communities_loop (set subjectFiles to 0)
% - check_subject_clustering_loop (set subjectFiles to 0)
% - check_subject_participation_loop (set subjectFiles to 0)
% - check_subject_density_loop
% - check_subject_spatial_metrics_loop
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
filenames{3} = ['dataSet' num2str(dataSet) '_modularity.mat'];
filenames{4} = ['dataSet' num2str(dataSet) '_modularity_VA.mat'];
filenames{5} = ['dataSet' num2str(dataSet) '_clustering.mat'];
filenames{6} = ['dataSet' num2str(dataSet) '_participation.mat'];
filenames{7} = ['dataSet' num2str(dataSet) '_gateway.mat'];
filenames{8} = ['dataSet' num2str(dataSet) '_diversity.mat'];
filenames{9} = ['dataSet' num2str(dataSet) '_density.mat'];
filenames{10} = ['dataSet' num2str(dataSet) '_comm_radius.mat'];
filenames{11} = ['dataSet' num2str(dataSet) '_emodiversity.mat'];
for i_file = 1:numel(filenames)
    load(filenames{i_file});
end
checks = readtable('validityChecks_granularity.csv'); 
if dataSet == 18
    rawData = importdata('18ARIVariations_QuestionnaireData.xlsx');
    TAStotal = rawData.data(:,17); %39 for session 2
    TASddf = rawData.data(:,15); %37 for session 2
    TASdif = rawData.data(:,16); %38 for session 2
    TASeot = rawData.data(:,14); %36 for session 2
    range = rawData.data(:,19)/7; %44 for session 2 % divide by 7 to get mean per item
    diff = rawData.data(:,18)/7; %43 for session 2 % divide by 7 to get mean per item
    RDEES = rawData.data(:,20)/14; %45 for session 2 % divide by 14 to get mean per item
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
if dataSet == 39 % remove subject 31 (row 29) for all calculations
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
    TAStotal(29,:) = [];
    TASddf(29,:) = [];
    TASdif(29,:) = [];
    TASeot(29,:) = [];
    range(29,:) = [];
    diff(29,:) = [];
    RDEES(29,:) = [];
elseif dataSet == 88 % exclude subject ID in 500 range from all calculations
    clusterCoefGlobal = clusterCoefGlobal(1:length(rtoZ),:);
    diversity = diversity(1:length(rtoZ),:);
    gateway = gateway(1:length(rtoZ),:);
    networkDensity = networkDensity(1:length(rtoZ),:);
    numCommunities = numCommunities(1:length(rtoZ),:);
    modularity = modularity(1:length(rtoZ),:);
    modularity_VA = modularity_VA(1:length(rtoZ),:);
    participation = participation(1:length(rtoZ),:);
    commRadius = commRadius(1:length(rtoZ),:);
end
    
%% correlations between zICC and number of communities at thresholds
[zICCcorr(1,1),zICCcorr(1,2)] = corr(rtoZ(:,3),numCommunities(:,2),'rows','complete');
[zICCcorr(2,1),zICCcorr(2,2)] = corr(rtoZ(:,3),numCommunities(:,3),'rows','complete');
[zICCcorr(3,1),zICCcorr(3,2)] = corr(rtoZ(:,3),numCommunities(:,4),'rows','complete');

%% correlation between zICC and modularity values when communities are detected or assigned
[zICCcorr(45,1),zICCcorr(45,2)] = corr(rtoZ(:,3),modularity(:,2),'rows','complete');
[zICCcorr(46,1),zICCcorr(46,2)] = corr(rtoZ(:,3),modularity(:,3),'rows','complete');
[zICCcorr(47,1),zICCcorr(47,2)] = corr(rtoZ(:,3),modularity(:,4),'rows','complete');

[zICCcorr(48,1),zICCcorr(48,2)] = corr(rtoZ(:,3),modularity_VA(:,2),'rows','complete');
[zICCcorr(49,1),zICCcorr(49,2)] = corr(rtoZ(:,3),modularity_VA(:,3),'rows','complete');
[zICCcorr(50,1),zICCcorr(50,2)] = corr(rtoZ(:,3),modularity_VA(:,4),'rows','complete');

%% correlations between zICC and global clustering coefficients at thresholds
[zICCcorr(4,1),zICCcorr(4,2)] = corr(rtoZ(:,3),clusterCoefGlobal(:,2),'rows','complete');
[zICCcorr(5,1),zICCcorr(5,2)] = corr(rtoZ(:,3),clusterCoefGlobal(:,3),'rows','complete');
[zICCcorr(6,1),zICCcorr(6,2)] = corr(rtoZ(:,3),clusterCoefGlobal(:,4),'rows','complete');

%% correlations between zICC and network density at thresholds
[zICCcorr(7,1),zICCcorr(7,2)] = corr(rtoZ(:,3),networkDensity(:,2),'rows','complete');
[zICCcorr(8,1),zICCcorr(8,2)] = corr(rtoZ(:,3),networkDensity(:,3),'rows','complete');
[zICCcorr(9,1),zICCcorr(9,2)] = corr(rtoZ(:,3),networkDensity(:,4),'rows','complete');

[zICCcorr(42,1),zICCcorr(42,2)] = corr(rtoZ(:,3),networkDensity(:,5),'rows','complete');
[zICCcorr(43,1),zICCcorr(43,2)] = corr(rtoZ(:,3),networkDensity(:,6),'rows','complete');
[zICCcorr(44,1),zICCcorr(44,2)] = corr(rtoZ(:,3),networkDensity(:,7),'rows','complete');

%% correlations between zICC and participation coefficients at thresholds
[zICCcorr(10,1),zICCcorr(10,2)] = corr(rtoZ(:,3),participation(:,2),'rows','complete');
[zICCcorr(11,1),zICCcorr(11,2)] = corr(rtoZ(:,3),participation(:,3),'rows','complete');
[zICCcorr(12,1),zICCcorr(12,2)] = corr(rtoZ(:,3),participation(:,4),'rows','complete');

[zICCcorr(13,1),zICCcorr(13,2)] = corr(rtoZ(:,3),participation(:,5),'rows','complete');
[zICCcorr(14,1),zICCcorr(14,2)] = corr(rtoZ(:,3),participation(:,6),'rows','complete');
[zICCcorr(15,1),zICCcorr(15,2)] = corr(rtoZ(:,3),participation(:,7),'rows','complete');

%% correlations between zICC and gateway coefficients at thresholds
[zICCcorr(16,1),zICCcorr(16,2)] = corr(rtoZ(:,3),gateway(:,2),'rows','complete'); % based on node strength
[zICCcorr(17,1),zICCcorr(17,2)] = corr(rtoZ(:,3),gateway(:,3),'rows','complete'); % based on node strength
[zICCcorr(18,1),zICCcorr(18,2)] = corr(rtoZ(:,3),gateway(:,4),'rows','complete'); % based on node strength

[zICCcorr(19,1),zICCcorr(19,2)] = corr(rtoZ(:,3),gateway(:,5),'rows','complete'); % based on node strength
[zICCcorr(20,1),zICCcorr(20,2)] = corr(rtoZ(:,3),gateway(:,6),'rows','complete'); % based on node strength
[zICCcorr(21,1),zICCcorr(21,2)] = corr(rtoZ(:,3),gateway(:,7),'rows','complete'); % based on node strength

[zICCcorr(22,1),zICCcorr(22,2)] = corr(rtoZ(:,3),gateway(:,8),'rows','complete'); % based on betweenness centrality
[zICCcorr(23,1),zICCcorr(23,2)] = corr(rtoZ(:,3),gateway(:,9),'rows','complete'); % based on betweenness centrality
[zICCcorr(24,1),zICCcorr(24,2)] = corr(rtoZ(:,3),gateway(:,10),'rows','complete'); % based on betweenness centrality

[zICCcorr(25,1),zICCcorr(25,2)] = corr(rtoZ(:,3),gateway(:,11),'rows','complete'); % based on betweenness centrality
[zICCcorr(26,1),zICCcorr(26,2)] = corr(rtoZ(:,3),gateway(:,12),'rows','complete'); % based on betweenness centrality
[zICCcorr(27,1),zICCcorr(27,2)] = corr(rtoZ(:,3),gateway(:,13),'rows','complete'); % based on betweenness centrality

%% correlations between zICC and community radius at thresholds
[zICCcorr(51,1),zICCcorr(51,2)] = corr(rtoZ(:,3),commRadius(:,2),'rows','complete');
[zICCcorr(52,1),zICCcorr(52,2)] = corr(rtoZ(:,3),commRadius(:,3),'rows','complete');
[zICCcorr(53,1),zICCcorr(53,2)] = corr(rtoZ(:,3),commRadius(:,4),'rows','complete');

%% correlations between zICC and diversity coefficients at thresholds
[zICCcorr(28,1),zICCcorr(28,2)] = corr(rtoZ(:,3),diversity(:,2),'rows','complete');
[zICCcorr(29,1),zICCcorr(29,2)] = corr(rtoZ(:,3),diversity(:,3),'rows','complete');
[zICCcorr(30,1),zICCcorr(30,2)] = corr(rtoZ(:,3),diversity(:,4),'rows','complete');

[zICCcorr(31,1),zICCcorr(31,2)] = corr(rtoZ(:,3),diversity(:,5),'rows','complete');
[zICCcorr(32,1),zICCcorr(32,2)] = corr(rtoZ(:,3),diversity(:,6),'rows','complete');
[zICCcorr(33,1),zICCcorr(33,2)] = corr(rtoZ(:,3),diversity(:,7),'rows','complete');

%% correlation between zICC and emodiversity, alexithymia, range, and differentiation
[zICCcorr(34,1),zICCcorr(34,2)] = corr(rtoZ(:,3),emodiversity(:),'rows','complete');
if dataSet == 18 | dataSet == 39
    [zICCcorr(35,1),zICCcorr(35,2)] = corr(rtoZ(:,3),TAStotal(:),'rows','complete');
    [zICCcorr(36,1),zICCcorr(36,2)] = corr(rtoZ(:,3),TASddf(:),'rows','complete');
    [zICCcorr(37,1),zICCcorr(37,2)] = corr(rtoZ(:,3),TASdif(:),'rows','complete');
    [zICCcorr(38,1),zICCcorr(38,2)] = corr(rtoZ(:,3),TASeot(:),'rows','complete');
    [zICCcorr(39,1),zICCcorr(39,2)] = corr(rtoZ(:,3),range(:),'rows','complete');
    [zICCcorr(40,1),zICCcorr(40,2)] = corr(rtoZ(:,3),diff(:),'rows','complete');
    [zICCcorr(41,1),zICCcorr(41,2)] = corr(rtoZ(:,3),RDEES(:),'rows','complete');
end

%% organize and save output
correlations = array2table(zICCcorr);
correlations = horzcat(checks,correlations);
correlations.Properties.VariableNames = {'CheckName','v1','v2','gamma','threshold','r','p'};

if print == 1
    filename = ['dataSet' num2str(dataSet) '_construct_validity_checks_granularity_poly.xlsx'];
    writetable(correlations,filename);
end

%% Predicted relationships:
% more communities = higher granularity (lower ICC); negative correlation
% greater clustering = lower granularity (higher ICC); positive correlation

%% still to do
% could grab SDs for participation and clustering coefs
% gini coefficient has not been integrated