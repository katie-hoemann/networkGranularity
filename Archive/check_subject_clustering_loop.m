% This script generates, for every subject in the specified data set, a
% spreadsheet of node-level clustering coefficients at 3 pre-specified
% connection weight thresholds
% Across all subjects in the data set, this script generates a spreadsheet
% of the global (mean, network-level) clustering coefficients at each threshold,
% as well as histograms of the distribution at each threshold
% This script also generates a spreadsheet of the standard deviations of
% clustering coefficients at each threshold (no histograms)

clear;
clc;

%% specify dataset
dataSet = 88; % set to 18, 39, or 88
polychoric = 1; % set to 0 to use standard Pearson correlation matrix
subjectFiles = 0; % set to 1 to generate subject-(node-)level files

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

for i_subject = 1:length(subjectIDlist)
    %% generate matrix for subject
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
    
    %% calculate clustering coefficient on original connection weights
    [C,C_null,Ctot,Ctot_null] = clustering_coef_wu_sign(subjectMatrix,3);
    clusterCoefGlobal(i_subject,1) = Ctot; % save off global value
    clusterCoefNode(:,1) = C; % save off node level vector
    
    %% clustering coefficient on thresholded connection weights, first test value
    threshold = 0.1; % removes correlations below given absolute value
    subjectMatrixT1 = subjectMatrix;
    for i_row = 1:length(subjectMatrixT1)
        for i_col = 1:length(subjectMatrixT1)
            if abs(subjectMatrixT1(i_row,i_col)) < threshold
                subjectMatrixT1(i_row,i_col) = 0;
            end
        end
    end
    [C,C_null,Ctot,Ctot_null] = clustering_coef_wu_sign(subjectMatrixT1,3);
    clusterCoefGlobal(i_subject,2) = Ctot; % save off global value
    clusterCoefNode(:,2) = C; % save off node level vector

    %% clustering coefficient on thresholded connection weights, second test value
    threshold = 0.2; % removes correlations below given absolute value
    subjectMatrixT2 = subjectMatrix;
    for i_row = 1:length(subjectMatrixT2)
        for i_col = 1:length(subjectMatrixT2)
            if abs(subjectMatrixT2(i_row,i_col)) < threshold
                subjectMatrixT2(i_row,i_col) = 0;
            end
        end
    end
    [C,C_null,Ctot,Ctot_null] = clustering_coef_wu_sign(subjectMatrixT2,3);
    clusterCoefGlobal(i_subject,3) = Ctot; % save off global value
    clusterCoefNode(:,3) = C; % save off node level vector

    %% save variance (standard deviation) for node-level measure
    clusterCoefGlobalSD(i_subject,:) = std(clusterCoefNode);
    
    %% write results for node-level clustering coefficient per subject
    if subjectFiles == 1
        filename = ['subject' num2str(subjectID) '_clustering_coefficient_nodes.xlsx'];
        clusteringNode = array2table(clusterCoefNode);
        rowNames = cell2table(remainingNodes);
        clusteringNode = horzcat(rowNames,clusteringNode);
        clusteringNode.Properties.VariableNames = {'Emotion','OriginalData','Threshold1','Threshold2'};
        writetable(clusteringNode,filename);
    end
    
    %% clear variables
    clear clusterCoefNode;
end

%% write results for global clustering coefficient and distribution across subjects
clusterCoefGlobal = horzcat(subjectIDlist,clusterCoefGlobal);
colNames = {'SubjectID','OriginalData','Threshold1','Threshold2'};
clusteringGlobal = array2table(clusterCoefGlobal,'VariableNames',colNames);
writetable(clusteringGlobal,'clustering_coefficient_global_distribution.xlsx');
filename = ['dataSet' num2str(dataSet) '_clustering.mat'];
save(filename,'clusterCoefGlobal');

clusterCoefGlobalSD = horzcat(subjectIDlist,clusterCoefGlobalSD);
colNames = {'SubjectID','OriginalData','Threshold1','Threshold2'};
clusteringGlobalSD = array2table(clusterCoefGlobalSD,'VariableNames',colNames);
writetable(clusteringGlobalSD,'clustering_coefficient_global_SD_distribution.xlsx');
filename = ['dataSet' num2str(dataSet) '_clustering_SD.mat'];
save(filename,'clusterCoefGlobalSD');

figure;
histogram1 = histogram(clusterCoefGlobal(:,2));
xlim([0 max(clusterCoefGlobal(:,4))+1]);
ylim([1 length(subjectIDlist)]);
title('Global clustering coefficient in original data');
saveas(histogram1,'Clustering_coefficient_global_distribution_original','tiff');

figure;
histogram2 = histogram(clusterCoefGlobal(:,3));
xlim([0 max(clusterCoefGlobal(:,4))+1]);
ylim([1 length(subjectIDlist)]);
title('Global clustering coefficient at first threshold');
saveas(histogram2,'Clustering_coefficient_global_distribution_threshold1','tiff');

figure;
histogram3 = histogram(clusterCoefGlobal(:,4));
xlim([0 max(clusterCoefGlobal(:,4))+1]);
ylim([1 length(subjectIDlist)]);
title('Global clustering coefficient at second threshold');
saveas(histogram3,'Clustering_coefficient_global_distribution_threshold2','tiff');

close all;

