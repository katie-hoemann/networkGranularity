clear;
clc;

%% load data file, set tuning parameter, initialize variables
datafile = 'ARIEOD_n17_040718.xlsx';
rawdata = importdata(datafile);
subjectIDlist = unique(rawdata.data(:,1)); % grab subject IDs from first column
wordlist = rawdata.colheaders(2:end)'; % grab sampled words from top row
corrMatrix = zeros(length(wordlist),length(wordlist),length(subjectIDlist));
rawICC_mat = zeros(length(subjectIDlist),3);
rtoZ_mat = zeros(length(subjectIDlist),3);
gMin = 0; % set starting value for tuning parameter gamma
gMax = 5; % set ending value for tuning parameter gamma
gStepSize = .05; % set step size for tuning parameter gamma
gArray = [gMin:gStepSize:gMax];
nGamma = length(gArray);
nFold = 5; % set number of folds for cross-validation partition

%% set valence and arousal categories for sampled words
wordfile = 'words18.csv'; %load word list that includes raw norms
words = readtable(wordfile); 

for i_word=1:height(words) % define valence categories
if words.Valence(i_word)>5
    valCat(i_word)={'Positive'};
    positive(i_word)=1;
else
    valCat(i_word)={'Negative'};
    positive(i_word)=0;
end
end 

for i_word=1:height(words) % define arousal categories
if words.Arousal(i_word)>4.6
    aroCat(i_word)={'High'};
    high(i_word)=1;
else
    aroCat(i_word)={'Low'};
    high(i_word)=0;
end
end 

words=[words valCat' aroCat']; % append table with category assignments
words.Properties.VariableNames(5:end)={'ValCat' 'AroCat'}; % label new variables
labels=[positive' high']; % create matrix for logical indexing in ICC commands

%% grab data for each subject, generate correlation matrix and ICC values
for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index = find(rawdata.data(:,1)==subjectID);
    subjectData = rawdata.data(index,2:end);
    %% remove missing data
    missing_data = isnan(subjectData); 
    missing_data2 = any(missing_data,2); 
    subjectData = subjectData(~missing_data2,:); 
    %% compute ICCs - positive, negative, and valence average
    rawICC_mat(i_subject,1) = ICC(subjectData(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
    rawICC_mat(i_subject,2) = ICC(subjectData(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
    rawICC_mat(i_subject,3) = (rawICC_mat(i_subject,1)+rawICC_mat(i_subject,2))/2; % valence average ICC
    %% Fisher transform ICCs
    rtoZ_mat(i_subject,1) = 0.5*log((1+rawICC_mat(i_subject,1))/(1-rawICC_mat(i_subject,1)));
    rtoZ_mat(i_subject,2) = 0.5*log((1+rawICC_mat(i_subject,2))/(1-rawICC_mat(i_subject,2)));
    rtoZ_mat(i_subject,3) = 0.5*log((1+rawICC_mat(i_subject,3))/(1-rawICC_mat(i_subject,3)));
    %% create correlation matrix
    subjectMatrix = corrcoef(subjectData); % calculate n-by-n correlation matrix
    subjectMatrix(logical(eye(size(subjectMatrix)))) = 0; % set diagonal to 0 for BCT functions
    %% save off correlation matrix before continuing
    corrMatrix(:,:,i_subject) = subjectMatrix;
end

%% partition data and run through cross-validation
CV=cvpartition(subjectIDlist,'KFold',nFold); % divide dataset into n training/test folds

for i_fold = 1:CV.NumTestSets
    trainingSet = find(CV.training(i_fold)==1); % index training set for this fold
    saveTrainingSets{i_fold} = trainingSet; % save off training set cases for this fold
    testSet = find(CV.test(i_fold)==1); % index test set for this fold
    saveTestSets{i_fold} = testSet; % save off test set cases for this fold
    %% find optimal gamma value for each training set based on correlation w ICC
    for i_subject = 1:length(trainingSet)
        subjectMatrixTrain = corrMatrix(:,:,trainingSet(i_subject)); 
        %% delete NaN values from correlation matrix
        missing_corr = find(isnan(subjectMatrixTrain(1,:)));
        subjectMatrixTrain(missing_corr,:) = []; % delete missing row
        subjectMatrixTrain(:,missing_corr) = []; % delete missing column
        %% community detection finetuning loop
        for i_gamma = 1:nGamma
            W = subjectMatrixTrain;  % set subject matrix
            n = size(W,1);      % number of nodes
            M = 1:n;            % initial community affiliations
            Q0 = -1; Q1 = 0;    % initialize modularity values
                while Q1 - Q0 > 1e-5    % while modularity increases
                    Q0 = Q1;            % perform community detection
                    [M, Q1] = community_louvain(W,gArray(i_gamma),[],'negative_sym');
                end
            nodes(i_subject) = n;       % save number of nodes
            comTrain(i_subject,i_gamma) = max(M);   % save number of communities
            modTrain(i_subject,i_gamma) = Q1;         % save modularity (q value)
            if max(M) == n              % stop if number of communities equals number of nodes
                break;
            end
        end 
        clear subjectMatrixTrain
        clear missing_corr
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
    ICCtoModularity = corr(rtoZ_mat(find(CV.training(i_fold)==1),3),modTrain);
    ICCtoCommunity = corr(rtoZ_mat(find(CV.training(i_fold)==1),3),comTrain);
    
    clear modTrain
    clear comTrain

    corrModMax = max(ICCtoModularity);                  % find maximum positive correlation between ICC and modularity (expected direction)
    corrModMin = min(ICCtoModularity);                  % find maximum negative correlation between ICC and modularity
    nValues = numel(gArray(ICCtoModularity==corrModMax));
    gammaOptimalbyModCorrPos(i_fold,1:nValues) = gArray(ICCtoModularity==corrModMax);    % find gamma at which positive correlation is highest
    nValues = numel(gArray(ICCtoModularity==corrModMin));
    gammaOptimalbyModCorrNeg(i_fold,1:nValues)= gArray(ICCtoModularity==corrModMin);    % find gamma at which negative correlation is highest
    corrComMax = max(ICCtoCommunity);                   % find maximum positive correlation between ICC and number of communities
    corrComMin = min(ICCtoCommunity);                   % find maximum negative correlation between ICC and number of communities (expected direction)
    nValues = numel(gArray(ICCtoCommunity==corrComMax));
    gammaOptimalbyComCorrPos(i_fold,1:nValues) = gArray(ICCtoCommunity==corrComMax);    % find gamma at which positive correlation is highest
    nValues = numel(gArray(ICCtoCommunity==corrComMin));
    gammaOptimalbyComCorrNeg(i_fold,1:nValues) = gArray(ICCtoCommunity==corrComMin);    % find gamma at which negative correlation is highest
    
    %% run community detection on test set using optimal gamma identified in training set
    for i_subject = 1:length(testSet)
        subjectMatrixTest = corrMatrix(:,:,testSet(i_subject)); 
        %% delete NaN values from correlation matrix
        missing_corr = find(isnan(subjectMatrixTest(1,:)));
        subjectMatrixTest(missing_corr,:) = [];
        subjectMatrixTest(:,missing_corr) = [];
        %% run community detection algorithm using gamma value identified in training
            %for i_gTest = 1:numel(gammaOptimalbyComCorrNeg(i_fold,:)); % for loop in case > 1 gamma
                W = subjectMatrixTest;  % set subject matrix
                n = size(W,1);      % number of nodes
                M = 1:n;            % initial community affiliations
                Q0 = -1; Q1 = 0;    % initialize modularity values
                    while Q1 - Q0 > 1e-5    % while modularity increases
                        Q0 = Q1;            % perform community detection
                        [M, Q1] = community_louvain(W,mean(gammaOptimalbyComCorrNeg(i_fold)),[],'negative_sym');
                        %[M, Q1] = community_louvain(W,gammaOptimalbyComCorrNeg(i_fold,i_gTest),[],'negative_sym');
                    end
                %comTest(i_subject,i_gTest) = max(M);   % save number of communities
                %modTest(i_subject,i_gTest) = Q1;         % save modularity (q value)
                comTest(i_subject) = max(M);   % save number of communities
                modTest(i_subject) = Q1;         % save modularity (q value)
            %end
   
        clear subjectMatrixTest
        clear missing_corr
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

%% determine value of gamma where variance of modularity is minimized across test sets
varModMinTest = min(varModTest);                     % find minimum variance of q
gammaOptimalbyModVar = gammaOptimalbyComCorrNeg(varModTest==varModMinTest,:);
gammaOptimalbyModVar = mean(gammaOptimalbyModVar);    % find gamma at which variance of q is minimum

%% generate plots for mean modularity with variance error bars
gammaToPlot = gammaOptimalbyComCorrNeg; % fixing the way I generated the variable
for i_row = 1:length(gammaToPlot);
    if numel(gammaToPlot(i_row,:)) > 1;
        gammaToPlot(i_row,1) = mean(gammaToPlot(i_row,:)); % only want one gamma value per fold for plotting
    end
end
errorbar(gammaToPlot(:,1),meanModTest,varModTest); 
title ('Mean Modularity at Test (+/- Variance) by Optimized Gamma Value')
xlabel ('Gamma Value')
ylabel ('Modularity')
line ([gammaOptimalbyModVar gammaOptimalbyModVar], get(gca, 'ylim'), 'Color','r');