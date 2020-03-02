%% Instructions:
% run the following scripts for the 18-term dataset:
% - assign_subject_communities_loop (set subjectFiles to 0)
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
threshold = 2; % set to 0, 1, or 2

%% load network measures
filenames{1} = ['dataSet' num2str(dataSet) '_zICC.mat'];
filenames{2} = ['dataSet' num2str(dataSet) '_communities.mat'];
filenames{3} = ['dataSet' num2str(dataSet) '_modularity.mat'];
filenames{4} = ['dataSet' num2str(dataSet) '_modularity_VA.mat'];
filenames{5} = ['dataSet' num2str(dataSet) '_clustering.mat'];
filenames{6} = ['dataSet' num2str(dataSet) '_participation.mat'];
filenames{7} = ['dataSet' num2str(dataSet) '_gateway.mat'];
filenames{8} = ['dataSet' num2str(dataSet) '_diversity.mat'];
filenames{9} = ['dataSet' num2str(dataSet) '_density.mat'];
filenames{10} = ['dataSet' num2str(dataSet) '_emodiversity.mat'];
for i_file = 1:numel(filenames)
    load(filenames{i_file});
end

%% load psychological measures
rawData = importdata('18ARIVariations_QuestionnaireData.xlsx');
ASItotal = rawData.data(:,7);
ERStotal = rawData.data(:,11); 
GADtotal = rawData.data(:,12); 
PSStotal = rawData.data(:,21); %46 for session 2
PHQ15total = rawData.data(:,22); %47 for session 2
PHQ8total = rawData.data(:,23); %48 for session 2
neuroticism = rawData.data(:,49);
TAStotal = rawData.data(:,17); %39 for session 2
TASddf = rawData.data(:,15); %37 for session 2
TASdif = rawData.data(:,16); %38 for session 2
TASeot = rawData.data(:,14); %36 for session 2
range = rawData.data(:,19); %44 for session 2
diff = rawData.data(:,18); %43 for session 2
RDEES = rawData.data(:,20); %45 for session 2

%% define variable sets
if threshold == 0
    %variableSetA = [numCommunities(:,2) modularity_VA(:,2) clusterCoefGlobal(:,2) networkDensity(:,5) participation(:,2) participation(:,5) diversity(:,2) diversity(:,5)];
    variableSetA = [numCommunities(:,2) modularity_VA(:,2) clusterCoefGlobal(:,2) networkDensity(:,5) diversity(:,2) diversity(:,5)];
elseif threshold == 1
    %variableSetA = [numCommunities(:,3) modularity_VA(:,3) clusterCoefGlobal(:,3) networkDensity(:,6) participation(:,3) participation(:,6) diversity(:,3) diversity(:,6)];
    variableSetA = [numCommunities(:,3) modularity_VA(:,3) clusterCoefGlobal(:,3) networkDensity(:,6) diversity(:,3) diversity(:,6)];
elseif threshold == 2
    variableSetA = [numCommunities(:,4) modularity_VA(:,4) clusterCoefGlobal(:,4) networkDensity(:,7)];
end
variableSetB = [ASItotal ERStotal GADtotal PSStotal PHQ15total PHQ8total TAStotal RDEES];
%% set empty values to 0
variableSetA(isnan(variableSetA)) = 0;
variableSetB(isnan(variableSetB)) = 0;

%% run canonical correlation
[A,B,r,U,V,stats] = canoncorr(variableSetA,variableSetB);

%% plot canonical variables scores for each function
for i = 1:length(A)
    figure;
    plot(U(:,i),V(:,i),'.');
end

%% generate table of canonical coefficients for each variable set
% canonCoeffA = array2table(A);
% rowNamesA = cell2table(rawDataA.colheaders(3:end)');
% canonCoeffA = horzcat(rowNamesA,canonCoeffA);
% canonCoeffA = sortrows(canonCoeffA,2,'descend');
% 
% canonCoeffB = array2table(B);
% rowNamesB = cell2table(rawDataB.colheaders(3:end)');
% canonCoeffB = horzcat(rowNamesB,canonCoeffB);
% canonCoeffB = sortrows(canonCoeffB,2,'descend');

