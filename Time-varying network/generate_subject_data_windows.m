% This script generates, for a given subject in a given data set, the
% correlation matrix and community assignment vector for use in
% dynamic (time-varying) visualization (e.g., in Gephi)

clear;
clc;

%% specify dataset and subject number, whether to print results to file, version parameters
dataSet = 39; % set to 18 or 39
subjectID = 40; %in DS18: 19 and 21 have high gran/4 com; 55 and 56 have low gran/2 com; in DS39: 5, 7, and 17 have low/2; 2, 32, and 40 have high/4; 38 has 5 com
print = 1; % set to 1 to print subject matrix and community assignment

%% load data file, along with word file that includes raw norms
if dataSet == 18
    dataFile = '18ARIEOD_daily_filtered.xlsx';
    wordFile = 'words18.csv'; 
    rawData = importdata(dataFile);
    allData = rawData.data.NoLateSurveys;
    subjectIDlist = unique(allData(:,1)); % grab subject IDs from first column of data file
    words = readtable(wordFile); 
    wordList = rawData.colheaders.NoLateSurveys(5:end)';  % grab sampled words from top row of data file
elseif dataSet == 39
    dataFile = '39ESdata_daily_filtered.xlsx';
    wordFile = 'words39.csv';
    rawData = importdata(dataFile);
    allData = rawData.data.Sheet1;
    subjectIDlist = unique(allData(:,1)); % grab subject IDs from first column of data file
    words = readtable(wordFile); 
    wordList = rawData.colheaders.Sheet1(5:end)';  % grab sampled words from top row of data file
end

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
index1 = find(allData(:,1)==subjectID);
subjectData = allData(index1,2:end);
dayIDlist = unique(subjectData(:,1));
for i_window = 1:length(dayIDlist)-2
    focusDay = i_window+1;
    windowData = [];
    index2 = [find(subjectData(:,1)==(focusDay-1)); find(subjectData(:,1)==(focusDay)); find(subjectData(:,1)==(focusDay+1))];
    windowData = subjectData(index2,4:end);
    % remove invalid values and missing data
    windowData(windowData > maximumRating) = NaN;
    missingData = isnan(windowData); 
    missingData2 = any(missingData,2); 
    windowData = windowData(~missingData2,:); 
    if numel(windowData) == 0
        continue
    end
    % if necessary, rescale data to start at 0
    if startRatingsat1 == 1
        windowData = windowData-1;
    end
    windowMatrix = corrcoef(windowData); % calculate n-by-n correlation matrix
    % delete NaN values from correlation matrix
    i_corr = 1;
    missingCorr = find(isnan(windowMatrix(i_corr,:)));
    while numel(missingCorr) == length(wordList)
        i_corr = i_corr+1;
        missingCorr = find(isnan(windowMatrix(i_corr,:)));
    end
    deletedCorr = isnan(windowMatrix(1,:)); % save off which row/col were removed
    remainingNodes = wordList;
    remainingNodes(missingCorr) = [];
    windowMatrix(missingCorr,:) = []; % delete missing rows
    windowMatrix(:,missingCorr) = []; % delete missing columns
    windowMatrix(logical(eye(size(windowMatrix)))) = 0; % set diagonal to 0 for BCT functions
    % run community detection, iterating 1000x to maximize Q
    for i_iter = 1:1000
        gamma = 1;
        W = windowMatrix;  % set subject matrix
        B = 'negative_asym';
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
    if numel(maxModularity_Index) == 0
        maxModularity_Index = 1;
    end
    communityAssignment = communityAssignment_Iter(:,maxModularity_Index); % save community assignment vector    
    numCommunities(i_window) = max(communityAssignment);   % save number of communities

    %% prep data for use in Gephi
    % prune edges for display using backbone
    [nRow,nCol] = size(windowMatrix); % get size of original matrix for later use
    windowMatrixNeg = windowMatrix<0; % capture where negative weights occurred in original matrix
    windowMatrixAbs = abs(windowMatrix); % create positive-only matrix
    [~,windowMatrixAbs] = backbone_wu(windowMatrixAbs,round(nCol-(nCol/2))); % create network backbone (positive edges only)
    for i_cell = 1:numel(windowMatrixAbs) % reset negative edge weights where they occurred in the original matrix
        if windowMatrixNeg(i_cell) == 1
            windowMatrixBackbone(i_cell) = windowMatrixAbs(i_cell)*-1;
        else
            windowMatrixBackbone(i_cell) = windowMatrixAbs(i_cell);
        end
    end
    windowMatrixBackbone = reshape(windowMatrixBackbone,nRow,nCol); % reshape matrix to original dimensions
    G1 = graph(windowMatrixBackbone,remainingNodes,'upper'); % create graph object (use plot function to visualize)
    negEdges = G1.Edges.Weight<0; % identify negative edges
    windowEdgesBackbone = G1.Edges; % return the edge list
    isNeg = array2table(negEdges,'VariableNames',{'IsNeg'}); % create array for negative weight index
    absWeight = abs(table2array(windowEdgesBackbone(:,2))); % create array for absolute value weight
    absWeight = array2table(absWeight,'VariableNames',{'absWeight'});
    windowEdgesBackbone = [windowEdgesBackbone absWeight isNeg]; % add arrays to edge list
    windowEdgesBackbone(:,2) = []; % remove unnecessary variable
    windowEdgesBackbone = splitvars(windowEdgesBackbone); % split nested edges columns
    windowEdgesBackbone.Properties.VariableNames = {'Source' 'Target' 'Weight' 'IsNeg'}; % relabel columns for Gephi import
    
    % binarize network and remove negative edges
    windowMatrixPos = windowMatrix>0; % create binarized matrix of only positive edges
    G2 = graph(windowMatrixPos,remainingNodes); % create graph object
    windowEdgesPos = G2.Edges; % return the edge list
    windowEdgesPos = splitvars(windowEdgesPos,'EndNodes','NewVariableNames',{'Source','Target'}); % split nested edge columns
    windowTimeStamp = repmat(datetime(2020,1,i_window,0,0,0,'Format','uuuu-MM-dd''T''HH:mm:ss'),1,height(windowEdgesPos))'; % create timestamp column
    windowEdgesPos = table2timetable(windowEdgesPos,'RowTimes',windowTimeStamp); % add timestamp column to edge list
    windowEdgesPos = timetable2table(windowEdgesPos);
    windowEdgesPos = movevars(windowEdgesPos,'Time','After','Target');
    windowEdgesPos.Properties.VariableNames = {'Source' 'Target' 'Timeset'}; % relabel columns for Gephi import
    
    %% write results
    rowNames = remainingNodes;
    colNames = remainingNodes;
    matrix = array2table(windowMatrix,'RowNames',rowNames,'VariableNames',colNames);
    communities = array2table(communityAssignment);
    communities = horzcat(rowNames,lower(rowNames),communities); % add columns for node ID and name
    communities.Properties.VariableNames = {'Id' 'Label' 'Community'}; % relabel columns for Gephi import
        
    if print == 1
        %writetable(matrix,['subject' num2str(subjectID) '_window' num2str(i_window) '_correlation_matrix.xlsx']);
        writetable(windowEdgesBackbone,['subject' num2str(subjectID) '_window' num2str(i_window) '_edge_list_backbone.xlsx']);
        writetable(communities,['subject' num2str(subjectID) '_window' num2str(i_window) '_community_assignment.xlsx']);
    end

    if i_window == 1
        dynamicEdgesPos = windowEdgesPos;
    else
        dynamicEdgesPos = vertcat(dynamicEdgesPos,windowEdgesPos);
    end
    clear deletedCorr remainingNodes communityAssignment_Iter communityAssignment windowMatrixBackbone
end

if print == 1
    writetable(dynamicEdgesPos,['subject' num2str(subjectID) '_dynamic_edge_list']); % write to .txt to preserve timestamp format
end