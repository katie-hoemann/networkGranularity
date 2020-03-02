clear;
clc;

%% specify dataset
dataSet = 18; % set to 18, 39, or 88 and the rest will take care of itself!

%% load data file, along with word file that includes raw norms
if dataSet == 18
    dataFile = '18ARIEOD.xlsx';
    wordFile = 'words18.csv';
elseif dataSet == 39
    dataFile = '39ES_fixed.xlsx';
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
    
    %% Compute the network metrics based on spatial location
    gammaIndex = find(gArray == 1);
    subjectIndex = i_subject;
    keepNode = ~deletedCorr(subjectIndex,:);
    locations = table2array(words(keepNode,[3 4]));
    partitions = communityAssignment{subjectIndex, gammaIndex};
    
    comm_spatial_diameter_array = comm_spatial_diameter(partitions,locations);
    comm_spatial_diameter_value(subjectIndex) = comm_spatial_diameter_array(end,end);
    comm_radius_array = comm_radius(partitions,locations);
    comm_radius_value(subjectIndex) = comm_radius_array(end,end);
    comm_ave_pairwise_spatial_dist_array = comm_ave_pairwise_spatial_dist(partitions,locations);
    comm_ave_pairwise_spatial_dist_value(subjectIndex) = comm_ave_pairwise_spatial_dist_array(end,end);
    
    clearvars comm_spatial_diameter_array comm_radius_array comm_ave_pairwise_spatial_dist_array
end

%% remove subjects that have a large number of deleted nodes
missingNodes = sum(deletedCorr,2); % sum how many nodes missing per subject
sparseSubjects = (length(wordList) - missingNodes) < minimumNetworkSize; % find which subjects are below threshold
removedSubjects = subjectIDlist(sparseSubjects); % save off which subjects were removed
modularity(sparseSubjects,:) = []; % delete missing row from modularity matrix
community(sparseSubjects,:) = []; % delete missing row from community matrix


