% This script generates sub-plotted histograms of raw connection weights
% (i.e., correlation values) for every subject in the given data set
% Subplot titles (subject IDs) that appear in red indicate that at least
% one node was removed from the correlation matrix due to lack of variance

clear;
clc;

%% specify dataset
dataSet = 18; % set to 18, 39, or 88 and the rest will take care of itself!
polychoric = 0; % set to 0 to use standard Pearson correlation matrix

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
corrMatrix = zeros(length(wordList),length(wordList),length(subjectIDlist));
deletedCorr = nan(length(subjectIDlist),length(wordList));
deletedNodes = cell(length(subjectIDlist),1);

%% grab data for each subject
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
        %% save off correlation matrix before continuing
        corrMatrix(:,:,i_subject) = subjectMatrix;
        %% delete NaN values from correlation matrix
        missingCorr = find(isnan(subjectMatrix(1,:)));
        deletedCorr(i_subject,:) = isnan(subjectMatrix(1,:)); % save off which row/col were removed
        deletedNodes{i_subject} = wordList(deletedCorr(i_subject,:)==1); % note which nodes are not in the network
        subjectMatrix(missingCorr,:) = []; % delete missing row
        subjectMatrix(:,missingCorr) = []; % delete missing column
    end
    subjectMatrix(logical(eye(size(subjectMatrix)))) = 0; % set diagonal to 0 for BCT functions

    %% histogram connection weights (correlations)
    halfSubjectMatrix = triu(subjectMatrix);
    halfSubjectMatrix = reshape(halfSubjectMatrix,[numel(halfSubjectMatrix),1]);
    halfSubjectMatrix = halfSubjectMatrix(halfSubjectMatrix~=0);
    subplot(7,8,i_subject);
    histogram(halfSubjectMatrix);
    ylim([0 100]);
    xlim([-1 1]);
    if numel(halfSubjectMatrix) < 153
        title(subjectID,'Color','r');
    else
        title(subjectID);
    end
%     %% create heatmap
%     subplot(7,8,i_subject);
%     colormap('parula')
%     imagesc(subjectMatrix);
%     colorbar;
%     caxis([-1 1]);
%     set(gca,'XTick',1:1:18,'XTickLabel',wordList);
%     xtickangle(45)
%     set(gca,'YTick',1:1:18,'YTickLabel',wordList);
%     if numel(halfSubjectMatrix) < 153
%         title(subjectID,'Color','r');
%     else
%         title(subjectID);
%     end
end