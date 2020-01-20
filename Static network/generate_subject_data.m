% This script generates, for a given subject in a given data set, the
% correlation matrix and community assignment vector for use in
% visualization (e.g., in Gephi)
% Optional parameters: binarize or threshold connection weights, etc.

clear;
clc;

%% specify dataset and subject number, whether to print results to file, version parameters
dataSet = 18; % set to 18, 39, or 88
subjectID = 38;
print = 0; % set to 1 to print subject matrix and community assignment
noNeg = 0; % sets negative weights to 0 (maintains positive weights)
binarize = 0; % sets negative weights to 0, positive weights to 1
polychoric = 0; % set to 0 to use standard Pearson correlation matrix
threshold = 0; % set > 0 to remove weak correlations below given absolute value
defaultGamma = 0; % set to 1 to use default gamma value of 1
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
if noNeg == 1
    for i_row = 1:length(subjectMatrix)
        for i_col = 1:length(subjectMatrix)
            if subjectMatrix(i_row,i_col) < 0
                subjectMatrix(i_row,i_col) = 0;
            end
        end
    end
end
if binarize == 1
    subjectMatrix = subjectMatrix > 0;
end
if threshold > 0
    for i_row = 1:length(subjectMatrix)
        for i_col = 1:length(subjectMatrix)
            if abs(subjectMatrix(i_row,i_col)) < threshold
                subjectMatrix(i_row,i_col) = 0;
            end
        end
    end
end

%% calculate clustering coefficient
[C,C_null,Ctot,Ctot_null] = clustering_coef_wu_sign(subjectMatrix,3);

%% community detection
if defaultGamma == 1
    gamma = 1;
elseif dataSet == 18 % based on cross-validated data 21.Nov.18
    if noNeg == 1
        gamma = 1.7;
    elseif binarize == 1
        gamma = 1.65;
%     elseif max(B == 'negative_sym') == 1
%         gamma = 1
    elseif threshold == 0.1
        gamma = 1.9; % based on cross-validated data 27.Mar.19
    elseif threshold == 0.2
        gamma = 1.9; % based on cross-validated data 27.Mar.19
    else
        gamma = 1.9; % based on cross-validated data 27.Mar.19
    end
elseif dataSet == 39 % based on cross-validated data 20.Mar.19
    if threshold == 0.1
        gamma = 1.1; % based on cross-validated data 20.Mar.19
    elseif threshold == 0.2
        gamma = 1.3; % based on cross-validated data 20.Mar.19
    else
        gamma = 1.1;
    end
elseif dataSet == 88
    if threshold == 0.1
        gamma = 1.45; % based on cross-validated data 6.Mar.19
    elseif threshold == 0.2
        gamma = 2.1; % based on cross-validated data 6.Mar.19
    else
        gamma = 1.45; % based on cross-validated data 6.Mar.19
    end
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
numCommunities = max(M);   % save number of communities

%% write results
rowNames = remainingNodes;
colNames = remainingNodes;
matrix = array2table(subjectMatrix,'RowNames',rowNames,'VariableNames',colNames);
communities = array2table(M,'RowNames',rowNames);
communities = sortrows(communities,1);

if print == 1
    writetable(matrix,'subject_correlation_matrix.xlsx');
    writetable(communities,'subject_community_assignment.xlsx');
end
