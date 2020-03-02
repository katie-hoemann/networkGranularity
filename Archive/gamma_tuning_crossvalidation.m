% This script runs gamma tuning with nFold cross-validation for the
% specified data set, and generates plots for evaluation
% Optional parameters: exclude negative connection weights, binarize
% connection weights, threshold connection weights, and change symmetrical
% treatment of negative v positive connection weights

clear;
clc;
seed = 1;
rng(seed);

%% specify dataset
dataSet = 18; % set to 18, 39, or 88 and the rest will take care of itself!

%% specify version
noNeg = 0; % sets negative weights to 0 (maintains positive weights)
binarize = 0; % sets negative weights to 0, positive weights to 1
threshold = 0.2; % clears weak correlation weights with absolute value at or below threshold
B = 'negative_asym'; % set to negative_sym, negative_asym, or modularity

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
nFold = 5; % set number of folds for cross-validation partition
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
deletedCorr = nan(length(subjectIDlist),length(wordList));
gArray = [gMin:gStepSize:gMax];
nGamma = length(gArray);
nodes = nan(length(subjectIDlist),1);
community = nan(length(subjectIDlist),nGamma);
modularity = nan(length(subjectIDlist),nGamma);
ICCtoCommunity = nan(nFold,nGamma);
diffICCtoCom = nan(nFold,nGamma-1);
elbow = nan(nFold,1);

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

%% screen data for each subject
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
    %% create correlation matrix
    subjectMatrix = corrcoef(subjectData); % calculate n-by-n correlation matrix
    subjectMatrix(logical(eye(size(subjectMatrix)))) = 0; % set diagonal to 0 for BCT functions
    %% delete NaN values from correlation matrix
    missingCorr = find(isnan(subjectMatrix(1,:)));
    deletedCorr(i_subject,:) = isnan(subjectMatrix(1,:)); % save off which row/col were removed
    subjectMatrix(missingCorr,:) = []; % delete missing row
    subjectMatrix(:,missingCorr) = []; % delete missing column
end
%% remove subjects that have a large number of deleted nodes
missingNodes = sum(deletedCorr,2); % sum how many nodes missing per subject
sparseSubjects = (length(wordList) - missingNodes) < minimumNetworkSize; % find which subjects are below threshold
removedSubjects = subjectIDlist(sparseSubjects); % save off which subjects were removed

%% initialize variables according to filtered subject list
filteredSubjectIDlist = subjectIDlist(~sparseSubjects==1);
rawICC = zeros(length(filteredSubjectIDlist),3);
rtoZ = zeros(length(filteredSubjectIDlist),3);
corrMatrix = zeros(length(wordList),length(wordList),length(filteredSubjectIDlist));
    
%% grab data for each included subject, compute ICC values, generate correlation matrix
for i_subject = 1:length(filteredSubjectIDlist)
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
    %% delete NaN values from correlation matrix
    missingCorr = find(isnan(subjectMatrix(1,:)));
    deletedCorr(i_subject,:) = isnan(subjectMatrix(1,:)); % save off which row/col were removed
    subjectMatrix(missingCorr,:) = []; % delete missing row
    subjectMatrix(:,missingCorr) = []; % delete missing column
    %% set negative weights to 0
    if noNeg == 1
        for i_row = 1:length(subjectMatrix)
            for i_col = 1:length(subjectMatrix)
                if subjectMatrix(i_row,i_col) < 0
                    subjectMatrix(i_row,i_col) = 0;
                end
            end
        end
    end
    %% binarize network
    if binarize == 1
        subjectMatrix = subjectMatrix > 0;
    end
    %% threshold network
    if threshold > 0
        for i_row = 1:length(subjectMatrix)
            for i_col = 1:length(subjectMatrix)
                if abs(subjectMatrix(i_row,i_col)) < threshold
                    subjectMatrix(i_row,i_col) = 0;
                end
            end
        end
    end
    %% save off correlation matrix before continuing
    corrMatrix(:,:,i_subject) = subjectMatrix;
end

%% partition data and run through cross-validation
CV=cvpartition(filteredSubjectIDlist,'KFold',nFold); % divide dataset into n training/test folds

for i_fold = 1:CV.NumTestSets
    trainingSet = find(CV.training(i_fold)==1); % index training set for this fold
    testSet = find(CV.test(i_fold)==1); % index test set for this fold
    %% find optimal gamma value for each training set based on correlation w ICC
    for i_subject = 1:length(trainingSet)
        subjectMatrixTrain = corrMatrix(:,:,trainingSet(i_subject)); 
        %% delete NaN values from correlation matrix
        missingCorr = find(isnan(subjectMatrixTrain(1,:)));
        subjectMatrixTrain(missingCorr,:) = []; % delete missing row
        subjectMatrixTrain(:,missingCorr) = []; % delete missing column
        %% community detection finetuning loop
        for i_gamma = 1:nGamma
            W = subjectMatrixTrain;  % set subject matrix
            n = size(W,1);      % number of nodes
            M = 1:n;            % initial community affiliations
            Q0 = -1; Q1 = 0;    % initialize modularity values
                while Q1 - Q0 > 1e-5    % while modularity increases
                    Q0 = Q1;            % perform community detection
                    [M, Q1] = community_louvain(W,gArray(i_gamma),[],B);
                end
            nodes(i_subject) = n;       % save number of nodes
            comAssignTrain{i_subject,i_gamma} = M; % save community assignment vector
            comTrain(i_subject,i_gamma) = max(M);   % save number of communities
            modTrain(i_subject,i_gamma) = Q1;         % save modularity (q value)
        end 
        clear subjectMatrixTrain
        clear missingCorr
    end
       
    %% calculate distributional statistics on the modularity values
    meanModTrain = mean(modTrain,1);
    sdModTrain = std(modTrain,1);
    varModTrain = var(modTrain,1);

    %% calculate distributional statistics on the number of communities
    meanComTrain = mean(comTrain,1);
    sdComTrain = std(comTrain,1);
    varComTrain = var(comTrain,1);
    minComTrain = min(comTrain);
    maxComTrain = max(comTrain);
    rangeComTrain = maxComTrain - minComTrain;

    %% compare modularity and number of communities against ICC
    ICCtoModularity = corr(rtoZ(find(CV.training(i_fold)==1),3),modTrain);
    ICCtoCommunity(i_fold,:) = corr(rtoZ(find(CV.training(i_fold)==1),3),comTrain);
      
    clear modTrain
    clear comTrain

    corrModMax = max(ICCtoModularity);                  % find maximum positive correlation between zICC and modularity (expected direction)
    corrModMin = min(ICCtoModularity);                  % find maximum negative correlation between zICC and modularity
    nValues = numel(gArray(ICCtoModularity==corrModMax));
    gammaOptimalbyModCorrPos(i_fold,1:nValues) = gArray(ICCtoModularity==corrModMax);    % find gamma at which positive correlation is highest
    nValues = numel(gArray(ICCtoModularity==corrModMin));
    gammaOptimalbyModCorrNeg(i_fold,1:nValues)= gArray(ICCtoModularity==corrModMin);    % find gamma at which negative correlation is highest
    corrComMax = max(ICCtoCommunity(i_fold,:));                   % find maximum positive correlation between zICC and number of communities
    corrComMin = min(ICCtoCommunity(i_fold,:));                   % find maximum negative correlation between zICC and number of communities (expected direction)
    nValues = numel(gArray(ICCtoCommunity(i_fold,:)==corrComMax));
    gammaOptimalbyComCorrPos(i_fold,1:nValues) = gArray(ICCtoCommunity(i_fold,:)==corrComMax);    % find gamma at which positive correlation is highest
    nValues = numel(gArray(ICCtoCommunity(i_fold,:)==corrComMin));
    gammaOptimalbyComCorrNeg(i_fold,1:nValues) = gArray(ICCtoCommunity(i_fold,:)==corrComMin);    % find gamma at which negative correlation is highest
        
    diffICCtoCom(i_fold,:) = diff(ICCtoCommunity(i_fold,:));              % take differential of correlation between zICC and number of communities
    elbow(i_fold) = max(diffICCtoCom(i_fold,:));                          % find maximum rate of increase (elbow)
    nValues = numel(gArray(diffICCtoCom(i_fold,:)==elbow(i_fold)));
    gammaOptimalbyElbow(i_fold,1:nValues) = gArray(diffICCtoCom(i_fold,:)==elbow(i_fold))+gStepSize;    % find gamma at which elbow occurs
    
    % plot correlation between community and zICC by gamma for this set
    figure;
    plot(gArray,ICCtoCommunity(i_fold,:));
    title (sprintf('Fold %d correlation between number of communities and zICC', i_fold))
    xlabel ('gamma')
    ylabel ('r')
    line ([gammaOptimalbyComCorrNeg(i_fold,1:nValues) gammaOptimalbyComCorrNeg(i_fold,1:nValues)], get(gca, 'ylim'), 'Color','b');
    line ([gammaOptimalbyComCorrPos(i_fold,1:nValues) gammaOptimalbyComCorrPos(i_fold,1:nValues)], get(gca, 'ylim'), 'Color','b');
    line ([gammaOptimalbyElbow(i_fold,1:nValues) gammaOptimalbyElbow(i_fold,1:nValues)], get(gca, 'ylim'), 'Color','r');
    
    %% run community detection on test set using optimal gamma identified in training set
    for i_subject = 1:length(testSet)
        subjectMatrixTest = corrMatrix(:,:,testSet(i_subject)); 
        %% delete NaN values from correlation matrix
        missingCorr = find(isnan(subjectMatrixTest(1,:)));
        subjectMatrixTest(missingCorr,:) = [];
        subjectMatrixTest(:,missingCorr) = [];
        %% run community detection algorithm using gamma value identified in training
            W = subjectMatrixTest;  % set subject matrix
            n = size(W,1);      % number of nodes
            M = 1:n;            % initial community affiliations
            Q0 = -1; Q1 = 0;    % initialize modularity values
                while Q1 - Q0 > 1e-5    % while modularity increases
                    Q0 = Q1;            % perform community detection
                    [M, Q1] = community_louvain(W,mean(gammaOptimalbyElbow(i_fold,:)),[],B);
                end
            comTest(i_subject) = max(M);   % save number of communities
            modTest(i_subject) = Q1;         % save modularity (q value)
   
        clear subjectMatrixTest
        clear missingCorr
    end
    
        %% calculate distributional statistics on the modularity values
        meanModTest(i_fold) = mean(modTest,2);
        sdModTest(i_fold) = std(modTest);
        varModTest(i_fold) = var(modTest);

        %% calculate distributional statistics on the number of communities
        meanComTest(i_fold) = mean(comTest,2);
        sdComTest(i_fold) = std(comTest);
        varComTest(i_fold) = var(comTest);
        minComTest(i_fold) = min(comTest);
        maxComTest(i_fold) = max(comTest);
        rangeComTest(i_fold) = maxComTest(i_fold) - minComTest(i_fold);
end

%% determine optimal values of gamma to minimize modularity variance
varModMinTest = min(varModTest);                     % find minimum variance of q
gammaOptimalbyModVar = gammaOptimalbyElbow(varModTest==varModMinTest,:);
gammaOptimalbyModVar = mean(gammaOptimalbyModVar);    % find gamma at which variance of q is minimum

%% generate plots
gammaToPlot = gammaOptimalbyElbow; % fixing the way I generated the variable
for i_row = 1:length(gammaToPlot);
    if numel(gammaToPlot(i_row,:)) > 1;
        gammaToPlot(i_row,1) = mean(gammaToPlot(i_row,:)); % only want one gamma value per fold for plotting
    end
end

meanICCtoCommunity = mean(ICCtoCommunity,1);
meanDiffICCtoCom = mean(diffICCtoCom,1);
meanElbow = max(meanDiffICCtoCom);                      
gammaOptimalbyMeanElbow = gArray(meanDiffICCtoCom==meanElbow)+gStepSize;

% plot correlation between community and zICC by gamma for all train sets
figure;
plot(gArray,ICCtoCommunity);
hold on
plot(gArray,meanICCtoCommunity,'Color','k','LineWidth',2);
title ('Correlation between number of communities and zICC')
xlabel ('gamma')
ylabel ('r')
line ([gammaOptimalbyMeanElbow gammaOptimalbyMeanElbow], get(gca, 'ylim'), 'Color','r');
hold off

figure;
errorbar(gammaToPlot(:,1),meanModTest,varModTest);
title ('Mean Modularity at Test (+/- Variance) by Optimized Gamma Value')
xlabel ('Gamma Value')
ylabel ('Modularity')
axis([(min(gammaToPlot(:,1))-.1) (max(gammaToPlot(:,1))+.1) (min(meanModTest)-.1) (max(meanModTest)+.1)])
line ([gammaOptimalbyModVar gammaOptimalbyModVar], get(gca, 'ylim'), 'Color','r');

figure;
errorbar(1:nFold,gammaToPlot(:,1),varModTest);
title ('Optimized Gamma Value (+/- Variance in Modularity) per Fold')
xlabel ('Fold')
ylabel ('Gamma Value')
axis([0 6 (min(gammaToPlot(:,1))-.1) (max(gammaToPlot(:,1))+.1)])
line (get(gca,'xlim'),[gammaOptimalbyModVar gammaOptimalbyModVar], 'Color','r');
