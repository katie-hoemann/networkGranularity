% This script runs gamma tuning WITHOUT crossvalidation for the
% specified data set, and generates plots for evaluation
% NOTE: This script is no longer being actively maintained

clear;
clc;

%% specify dataset
dataSet = 18; % set to 18, 39, or 88 and the rest will take care of itself!

%% load data file, along with word file that includes raw norms
if dataSet == 18
    dataFile = '18ARIEOD.xlsx';
    wordFile = 'words18.csv'; 
elseif dataSet == 39
    dataFile = '39ESdata.xlsx';
    wordFile = 'words39.csv';
else
    dataFile = '88PANASdata.xlsx';
    wordFile = 'words88.csv';
end
rawData = importdata(dataFile);
subjectIDlist = unique(rawData.data(:,1)); % grab subject IDs from first column of data file
words = readtable(wordFile); 
wordList = rawData.colheaders(2:end)';  % grab sampled words from top row of data file

%% set parameters
minimumNetworkSize = (3/4)*length(wordList);
gMin = 0.05; % set starting value for tuning parameter gamma; set > 0
gMax = 5; % set ending value for tuning parameter gamma
gStepSize = .05;  % set step size for tuning parameter gamma
if dataSet == 39
    startRatingsat1 = 1;
    maximumRating = 5; % set to highest valid scale value
elseif dataSet == 88
    startRatingsat1 = 0;
    maximumRating = 6;
    skipSubjectRange = 500; % set to invalid subject ID range
    subjectIDlist(subjectIDlist >= skipSubjectRange,:) = [];
else
    startRatingsat1 = 0;
    maximumRating = 6;
end

%% initialize variables
rawICC = zeros(length(subjectIDlist),3);
rtoZ = zeros(length(subjectIDlist),3);
corrMatrix = zeros(length(wordList),length(wordList),length(subjectIDlist));
deletedCorr = nan(length(subjectIDlist),length(wordList));
deletedNodes = cell(length(subjectIDlist),1);
gArray = [gMin:gStepSize:gMax];
nGamma = length(gArray);
nodes = nan(length(subjectIDlist),1);
community = nan(length(subjectIDlist),nGamma);
modularity = nan(length(subjectIDlist),nGamma);

%% set valence and arousal categories for sampled words
for i_word = 1:height(words) % define valence categories
if words.Valence(i_word) > 5 % derived based on the database mean for Warriner et al (2013)
    valCat(i_word) = {'Positive'};
    positive(i_word) = 1;
else
    valCat(i_word) = {'Negative'};
    positive(i_word) = 0;
end
end 

for i_word = 1:height(words) % define arousal categories
if words.Arousal(i_word) > 4.6 % derived based on the sample mean for 88 PANAS-X terms in Warriner et al (2013)
    aroCat(i_word) = {'High'};
    high(i_word) = 1;
else
    aroCat(i_word) = {'Low'};
    high(i_word) = 0;
end
end 

words = [words valCat' aroCat']; % append table with category assignments
words.Properties.VariableNames(5:end) = {'ValCat' 'AroCat'}; % label new variables
labels = [positive' high']; % create matrix for logical indexing in ICC commands

%% grab data for each subject and run through calculations
for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index = find(rawData.data(:,1)==subjectID);
    subjectData = rawData.data(index,2:end);
    %% remove invalid values and missing data
    subjectData(subjectData > maximumRating) = NaN;
    missingData = isnan(subjectData); 
    missingData2 = any(missingData,2); 
    subjectData = subjectData(~missingData2,:); 
    %% if necessary, rescale data to start at 0
    if startRatingsat1 == 1
        subjectData = subjectData-1;
    end
    %% compute ICCs - positive, negative, valence average
    rawICC(i_subject,1) = ICC(subjectData(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
    rawICC(i_subject,2) = ICC(subjectData(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
    rawICC(i_subject,3) = (rawICC(i_subject,1)+rawICC(i_subject,2))/2; % valence mean ICC
    %% Fisher transform ICCs (noted as zICCs)
    rtoZ(i_subject,1) = 0.5*log((1+rawICC(i_subject,1))/(1-rawICC(i_subject,1)));
    rtoZ(i_subject,2) = 0.5*log((1+rawICC(i_subject,2))/(1-rawICC(i_subject,2)));
    rtoZ(i_subject,3) = 0.5*log((1+rawICC(i_subject,3))/(1-rawICC(i_subject,3)));
    %% create correlation matrix
    subjectMatrix = corrcoef(subjectData); % calculate n-by-n correlation matrix
    subjectMatrix(logical(eye(size(subjectMatrix)))) = 0; % set diagonal to 0 for BCT functions
    %% save off correlation matrix before continuing
    corrMatrix(:,:,i_subject) = subjectMatrix;
    %% delete NaN values from correlation matrix
    missingCorr = find(isnan(subjectMatrix(1,:)));
    deletedCorr(i_subject,:) = isnan(subjectMatrix(1,:)); % save off which row/col were removed
    deletedNodes{i_subject} = wordList(deletedCorr(i_subject,:)==1); % note which nodes are not in the network
    subjectMatrix(missingCorr,:) = []; % delete missing row
    subjectMatrix(:,missingCorr) = []; % delete missing column
    %% calculate clustering coefficient
    [C,C_null,Ctot,Ctot_null] = clustering_coef_wu_sign(subjectMatrix,3);
    clusterCoefGlobal(i_subject) = Ctot; % save off global value
    clusterCoefNode{i_subject} = C; % save off node level vector
    %% community detection finetuning loop
    for i_gamma = 1:nGamma
            W = subjectMatrix;  % set subject matrix
            n = size(W,1);      % number of nodes
            M = 1:n;            % initial community affiliations
            Q0 = -1; Q1 = 0;    % initialize modularity values
                while Q1 - Q0 > 1e-5    % while modularity increases
                    Q0 = Q1;            % perform community detection
                    [M, Q1] = community_louvain(W,gArray(i_gamma),[],'negative_asym');
                end
            nodes(i_subject) = n;       % save number of nodes
            communityAssignment{i_subject,i_gamma} = M; % save community assignment vector
            community(i_subject,i_gamma) = max(M);   % save number of communities
            modularity(i_subject,i_gamma) = Q1;         % save modularity (q value)
    end
end

%% remove subjects that have a large number of deleted nodes
missingNodes = sum(deletedCorr,2); % sum how many nodes missing per subject
sparseSubjects = (length(wordList) - missingNodes) < minimumNetworkSize; % find which subjects are below threshold
removedSubjects = subjectIDlist(sparseSubjects); % save off which subjects were removed
modularity(sparseSubjects,:) = []; % delete missing row from modularity matrix
community(sparseSubjects,:) = []; % delete missing row from community matrix
rtoZclean = rtoZ; % create clean copy of rtoZ table
rtoZclean(sparseSubjects,:) = []; % delete missing row from clean rtoZ table

%% calculate distributional statistics on the modularity values
meanModularity = mean(modularity,1);
sdModularity = std(modularity,1);
varModularity = var(modularity,1);
varModMin = min(varModularity);                     % find minimum variance of q
gammaOptimalbyModVariance = gArray(varModularity==varModMin);    % find gamma at which variance of q is minimum

%% calculate distributional statistics on the number of communities
meanCommunity = mean(community,1);
sdCommunity = std(community,1);
varCommunity = var(community,1);
minCommunity = min(community);
maxCommunity = max(community);
rangeCommunity = maxCommunity - minCommunity;

%% run normality tests on distributions of number of communities
communitySWtest = nan(1,nGamma); % Shapiro-Wilk test of normality
for g = 1:nGamma
    communitySWtest(g) = swtest(community(:,g));
end

%% compare modularity and number of communities against zICC for valence mean, determine optimal values of gamma
ICCtoModularity = corr(rtoZclean(:,3),modularity);
ICCtoCommunity = corr(rtoZclean(:,3),community);

corrModMax = max(ICCtoModularity);                  % find maximum positive correlation between zICC and modularity (expected direction)
corrModMin = min(ICCtoModularity);                  % find maximum negative correlation between zICC and modularity
gammaOptimalbyModCorrPos = gArray(ICCtoModularity==corrModMax);    % find gamma at which positive correlation is highest
gammaOptimalbyModCorrNeg = gArray(ICCtoModularity==corrModMin);    % find gamma at which negative correlation is highest
corrComMax = max(ICCtoCommunity);                   % find maximum positive correlation between zICC and number of communities
corrComMin = min(ICCtoCommunity);                   % find maximum negative correlation between zICC and number of communities (expected direction)
gammaOptimalbyComCorrPos = gArray(ICCtoCommunity==corrComMax);    % find gamma at which positive correlation is highest
gammaOptimalbyComCorrNeg = gArray(ICCtoCommunity==corrComMin);    % find gamma at which negative correlation is highest

diffICCtoCom = diff(ICCtoCommunity);                % take differential of correlation between zICC and number of communities
elbow = max(diffICCtoCom);                          % find maximum rate of increase (elbow)
gammaOptimalbyElbow = gArray(diffICCtoCom==elbow)+gStepSize;    % find gamma at which elbow occurs

%% create tables of results
gammaTuning = table(gArray', meanCommunity', sdCommunity', rangeCommunity', communitySWtest', meanModularity', sdModularity', varModularity');
gammaTuning.Properties.VariableNames = {'Gamma' 'Mean_Num_Com' 'SD_Num_Com' 'Range_Num_Com' 'SW_Num_Com' 'Mean_Mod' 'SD_Mod' 'Var_Mod'};

%% generate plots of results
% plot for mean modularity with SD error bars
figure;
errorbar(gArray,meanModularity,sdModularity); 
title ('Mean modularity (+/- SD) by gamma value')
xlabel ('gamma')
ylabel ('q')

% plot for mean number of communities with SD error bars
figure;
errorbar(gArray,meanCommunity,sdCommunity); 
title ('Mean number of communities (+/- SD) by gamma value')
xlabel ('gamma')
ylabel ('Number of communities')

% plot of correlation between modularity and zICC by gamma value
figure;
plot(gArray,ICCtoModularity);
title ('Correlation between modularity and zICC by gamma value')
xlabel ('gamma')
ylabel ('r')
line ([gammaOptimalbyModCorrNeg gammaOptimalbyModCorrNeg], get(gca, 'ylim'), 'Color','b');
line ([gammaOptimalbyModCorrPos gammaOptimalbyModCorrPos], get(gca, 'ylim'), 'Color','b');
line ([gammaOptimalbyElbow gammaOptimalbyElbow], get(gca, 'ylim'), 'Color','r');

% plot of correlation between community and zICC by gamma value
figure;
plot(gArray,ICCtoCommunity);
title ('Correlation between number of communities and zICC by gamma value')
xlabel ('gamma')
ylabel ('r')
line ([gammaOptimalbyComCorrNeg gammaOptimalbyComCorrNeg], get(gca, 'ylim'), 'Color','b');
line ([gammaOptimalbyComCorrPos gammaOptimalbyComCorrPos], get(gca, 'ylim'), 'Color','b');
line ([gammaOptimalbyElbow gammaOptimalbyElbow], get(gca, 'ylim'), 'Color','r');

% plot of correlation between zICC and community, against range in community
figure;
scatter(ICCtoCommunity,rangeCommunity);
title ('Correlation between zICC and number of communities, against range in number of communities')
xlabel ('r')
ylabel ('Range in number of communities')

% plot relationships between zICC and community
figure;
scatter(rtoZclean(:,3),community(:,gArray==gammaOptimalbyElbow));
title ('Correlation between number of communities and zICC; expected negative')
xlabel ('zICC meanV')
ylabel ('Number of communities')

% plot relationship between zICC and modularity
figure;
scatter(rtoZclean(:,3),modularity(:,gArray==gammaOptimalbyElbow));
title ('Correlation between modularity and zICC; expected positive')
xlabel ('zICC meanV')
ylabel ('q')

% plot histogram of deleted nodes
figure;
histogram(sum(deletedCorr'));
title('Number of deleted nodes per participant')
xlabel('Deleted nodes')
ylabel('Total')


%% To-do EOD 27 February 2019: 
% - create separate script for gamma tuning from other calculations
% - clean all lists/matrices of removed participants?

%% Predicted relationships:
% expectations about number of communities will be based on n in dataset
% more communities = higher granularity (lower ICC); negative correlation
% more communities = lower modularity = lower ICC; positive correlation
% so we want to optimize gamma by ComCorrNeg and/or ModCorrPos
