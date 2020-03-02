% This script generates, for a given subject in a given data set, the
% community assignment vector at each of 3 pre-specified connection weight thresholds

clear;
clc;

%% specify dataset and subject number, whether to print results to file
dataSet = 18; % set to 18, 39, or 88
subjectID = 3;
weightInstanceIntensity = 0; % set to 1 to weight instances by total intensity
weightEmotionTotal = 0; % set to 1 to weight nodes by number of endorsements
weightInstanceTotal = 0; % set to 1 to weight instances by number of emotions endorsed
weightEmotionIntensity = 0; % set to 1 to weight nodes by maximum intensity
polychoric = 0; % set to 0 to use standard Pearson correlation matrix
defaultGamma = 1; % set to 1 to use default gamma value of 1
B = 'negative_asym'; % set to negative_sym, negative_asym, or modularity
print = 0; % set to 1 to print subject matrix and community assignment

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
if dataSet == 39
    startRatingsat1 = 1;
    maximumRating = 5; % set to highest valid scale value
else
    startRatingsat1 = 0;
    maximumRating = 6;
end

%% generate matrix for subject
subjectData = [];
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
%% if desired, weight instances by total intensity
if weightInstanceIntensity == 1
    for i_row = 1:length(subjectData)
        rowIntensity(i_row) = sum(subjectData(i_row,:));
        subjectData(i_row,:) = subjectData(i_row,:)./rowIntensity(i_row);
    end
end
%% if desired, weight nodes by number of endorsements
if weightEmotionTotal == 1
    for i_col = 1:size(subjectData,2)
        colTotal(i_col) = sum(subjectData(:,i_col)>0);
        subjectData(:,i_col) = subjectData(:,i_col).*colTotal(i_col);
    end
end
%% if desired, weight instances by number of emotions endorsed
if weightInstanceTotal == 1
    for i_row = 1:length(subjectData)
        rowTotal(i_row) = sum(subjectData(i_row,:)>0);
        subjectData(i_row,:) = subjectData(i_row,:)./rowTotal(i_row);
    end
end
%% if desired, weight nodes by maximum intensity
if weightEmotionIntensity == 1
    for i_col = 1:size(subjectData,2)
        colIntensity(i_col) = max(subjectData(:,i_col));
        subjectData(:,i_col) = subjectData(:,i_col)./colIntensity(i_col);
    end
end
%% create correlation matrix
if polychoric == 1
    %% delete zero-variance rows/columns from data
    colVar = var(subjectData,0,1); % find columns (nodes) with zero variance
    rowVar = var(subjectData,0,2); % find rows (instances) with zero variance
    missingColVar = find(colVar==0); % index zero-variance nodes
    missingRowVar = find(rowVar==0); % index zero-variance instances
    subjectData(:,missingColVar) = []; % remove zero-variance nodes
    subjectData(missingRowVar,:) = []; % remove zero-variance instances
    deletedNodes = wordList(missingColVar); % note which nodes are not in the network
    remainingNodes = wordList;
    remainingNodes(missingColVar) = [];
    subjectMatrix = polychoric_proc_missing(subjectData,NaN); % calculate n-by-n polychoric correlation matrix
else
    subjectMatrix = corrcoef(subjectData); % calculate n-by-n correlation matrix
    %% delete NaN values from correlation matrix
    missingCorr = find(isnan(subjectMatrix(1,:)));
    deletedCorr = isnan(subjectMatrix(1,:)); % save off which row/col were removed
    deletedNodes = wordList(deletedCorr==1); % note which nodes are not in the network
    remainingNodes = wordList(~deletedCorr==1);
    subjectMatrix(missingCorr,:) = []; % delete missing row
    subjectMatrix(:,missingCorr) = []; % delete missing column
end
subjectMatrix(logical(eye(size(subjectMatrix)))) = 0; % set diagonal to 0 for BCT functions

%% community detection on original connection weights
if defaultGamma == 1
    gamma = 1;
elseif dataSet == 18 % based on cross-validated data 27.Mar.19
    gamma = 1.9;
elseif dataSet == 39 % based on cross-validated data 20.Mar.19
    gamma = 1.1;
elseif dataSet == 88 % based on cross-validated data 6.Mar.19
    gamma = 1.45; 
else
    gamma = 1; % default value
end
W = subjectMatrix;  % set subject matrix
n = size(W,1);      % number of nodes
M = 1:n;            % initial community affiliations
Q0 = -1; Q = 0;    % initialize modularity values
    while Q - Q0 > 1e-5    % while modularity increases
        Q0 = Q;            % perform community detection
        [M, Q] = community_louvain(W,gamma,[],B); 
    end
communityAssignment(:,1) = M; % save community assignment vector    
numCommunities(1) = max(M);   % save number of communities

%% community detection on thresholded connection weights, first test value
threshold = 0.1; % removes correlations below given absolute value
subjectMatrixT1 = subjectMatrix;
for i_row = 1:length(subjectMatrixT1)
    for i_col = 1:length(subjectMatrixT1)
        if abs(subjectMatrixT1(i_row,i_col)) < threshold
            subjectMatrixT1(i_row,i_col) = 0;
        end
    end
end
if defaultGamma == 1
    gamma = 1;
elseif dataSet == 18 % based on cross-validated data 27.Mar.19
    gamma = 1.9;
elseif dataSet == 39 % based on cross-validated data 20.Mar.19
    gamma = 1.1;
elseif dataSet == 88 % based on cross-validated data 6.Mar.19
    gamma = 1.45; 
else
    gamma = 1; % default value
end
W = subjectMatrixT1;  % set subject matrix
n = size(W,1);      % number of nodes
M = 1:n;            % initial community affiliations
Q0 = -1; Q = 0;    % initialize modularity values
    while Q - Q0 > 1e-5    % while modularity increases
        Q0 = Q;            % perform community detection
        [M, Q] = community_louvain(W,gamma,[],B); 
    end
communityAssignment(:,2) = M; % save community assignment vector       
numCommunities(2) = max(M);   % save number of communities

%% community detection on thresholded connection weights, second test value
threshold = 0.2; % removes correlations below given absolute value
subjectMatrixT2 = subjectMatrix;
for i_row = 1:length(subjectMatrixT2)
    for i_col = 1:length(subjectMatrixT2)
        if abs(subjectMatrixT2(i_row,i_col)) < threshold
            subjectMatrixT2(i_row,i_col) = 0;
        end
    end
end
if defaultGamma == 1
    gamma = 1;
elseif dataSet == 18 % based on cross-validated data 27.Mar.19
    gamma = 1.9;
elseif dataSet == 39 % based on cross-validated data 20.Mar.19
    gamma = 1.3;
elseif dataSet == 88 % based on cross-validated data 6.Mar.19
    gamma = 2.1; 
else
    gamma = 1; % default value
end
W = subjectMatrixT2;  % set subject matrix
n = size(W,1);      % number of nodes
M = 1:n;            % initial community affiliations
Q0 = -1; Q = 0;    % initialize modularity values
    while Q - Q0 > 1e-5    % while modularity increases
        Q0 = Q;            % perform community detection
        [M, Q] = community_louvain(W,gamma,[],B); 
    end
communityAssignment(:,3) = M; % save community assignment vector       
numCommunities(3) = max(M);   % save number of communities

%% write results
rowNames = remainingNodes;
colNames = remainingNodes;
matrix = array2table(subjectMatrix,'RowNames',rowNames,'VariableNames',colNames);

communities = array2table(communityAssignment);
rowNames = cell2table(remainingNodes);
communities = horzcat(rowNames,communities);
communities.Properties.VariableNames = {'Emotion','OriginalData','Threshold1','Threshold2'};
communities = sortrows(communities,2);
    
if print == 1
%     writetable(matrix,'subject_correlation_matrix.xlsx');
    writetable(communities,'subject_community_assignment.xlsx');
end
