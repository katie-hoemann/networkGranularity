% This script generates, for a given subject in a given data set, the
% correlation matrix and community assignment vector for use in
% visualization (e.g., in Gephi)
% Optional parameters: binarize or threshold connection weights, etc.

clear;
clc;

%% specify dataset and subject number, whether to print results to file, version parameters
dataSet = 39; % set to 18, 39, or 88
subjectID = 38; %in DS18: 19 and 21 have high gran/4 com; 55 and 56 have low gran/2 com; in DS39: 5, 7, and 17 have low/2; 2, 32, and 40 have high/4; 38 has 5 com
print = 1; % set to 1 to print subject matrix and community assignment
noNeg = 0; % sets negative weights to 0 (maintains positive weights)
binarize = 0; % sets negative weights to 0, positive weights to 1
polychoric = 0; % set to 0 to use standard Pearson correlation matrix
threshold = 0; % set > 0 to remove weak correlations below given absolute value
defaultGamma = 1; % set to 1 to use default gamma value of 1
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
% run community detection, iterating 1000x to maximize Q
for i_iter = 1:1000
    W = subjectMatrix;  % set subject matrix
    n = size(W,1);      % number of nodes
    M = 1:n;            % initial community affiliations
    Q0 = -1; Q = 0;    % initialize modularity values
        while Q - Q0 > 1e-5    % while modularity increases
            Q0 = Q;            % perform community detection
            [M, Q] = community_louvain(W,gamma,[],B); 
        end
    communityAssignment_Iter(:,i_iter) = M; % save community assignment vector for given iteration
    modularity_Iter(i_iter) = Q; % save modularity value for given iteration
end
maxModularity_Index = find(modularity_Iter==max(modularity_Iter));
if numel(maxModularity_Index) > 1
    maxModularity_Index = maxModularity_Index(1);
end
communityAssignment = communityAssignment_Iter(:,maxModularity_Index); % save community assignment vector    
numCommunities = max(communityAssignment);   % save number of communities

%% prep data for use in Gephi
% prune edges for display using backbone
[nRow,nCol] = size(subjectMatrix); % get size of original matrix for later use
subjectMatrixNeg = subjectMatrix<0; % capture where negative weights occurred in original matrix
subjectMatrixAbs = abs(subjectMatrix); % create positive-only matrix
[~,subjectMatrixAbs] = backbone_wu(subjectMatrixAbs,round(nCol-(nCol/2))); % create network backbone (positive edges only)
for i_cell = 1:numel(subjectMatrixAbs) % reset negative edge weights where they occurred in the original matrix
    if subjectMatrixNeg(i_cell) == 1
        subjectMatrixBackbone(i_cell) = subjectMatrixAbs(i_cell)*-1;
    else
        subjectMatrixBackbone(i_cell) = subjectMatrixAbs(i_cell);
    end
end
subjectMatrixBackbone = reshape(subjectMatrixBackbone,nRow,nCol); % reshape matrix to original dimensions
G = graph(subjectMatrixBackbone,remainingNodes); % create graph object (use plot function to visualize)
negEdges = G.Edges.Weight<0; % identify negative edges
subjectEdgesBackbone = G.Edges; % return the edge list
isNeg = array2table(negEdges,'VariableNames',{'IsNeg'}); % create array for negative weight index
absWeight = abs(table2array(subjectEdgesBackbone(:,2))); % create array for absolute value weight
absWeight = array2table(absWeight,'VariableNames',{'absWeight'});
subjectEdgesBackbone = [subjectEdgesBackbone absWeight isNeg]; % add arrays to edge list
subjectEdgesBackbone(:,2) = []; % remove unnecessary variable
subjectEdgesBackbone = splitvars(subjectEdgesBackbone); % split nested edges columns
subjectEdgesBackbone.Properties.VariableNames = {'Source' 'Target' 'Weight' 'IsNeg'}; % relabel columns for Gephi import

%% write results
rowNames = remainingNodes;
colNames = remainingNodes;
matrix = array2table(subjectMatrix,'RowNames',rowNames,'VariableNames',colNames);
communities = array2table(communityAssignment);
communities = horzcat(rowNames,lower(rowNames),communities); % add columns for node ID and name
communities.Properties.VariableNames = {'Id' 'Label' 'Community'}; % relabel columns for Gephi import
if print == 1
    writetable(matrix,['subject' num2str(subjectID) '_correlation_matrix.xlsx']);
    writetable(subjectEdgesBackbone,['subject' num2str(subjectID) '_edge_list_backbone.xlsx']);
    writetable(communities,['subject' num2str(subjectID) '_community_assignment.xlsx']);
end
