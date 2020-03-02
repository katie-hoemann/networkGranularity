% This script generates, for every subject in the specified data set, a
% histogram of raw connection weights (i.e., correlation values), the number
% of positive, negative, and null links, a heatmap of the correlation matrix, 
% and sub-plotted histograms of raw intensity ratings for every emotion term
% Plot titles (subject IDs) that appear in red indicate that at least one
% node was removed from the correlation matrix due to lack of variance
% Generate spreadsheet of link counts for all participants in the data set

clear;
clc;

%% specify dataset and version parameters
dataSet = 18; % set to 18, 39, or 88
noNeg = 0;
binarize = 0;
polychoric = 1; % set to 0 to use standard Pearson correlation matrix
threshold = 0; % set > 0 to remove weak correlations below given absolute value

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

    %% histogram intensity ratings
    if dataSet == 18
        edges = [0 1 2 3 4 5 6 7];
        for i_emotion = 1:length(remainingNodes)
            subplot(6,3,i_emotion);
            intensityRatings = histogram(subjectData(:,i_emotion),edges);
            ylim([0 150]);
            xlim([0 7]);
            title(remainingNodes(i_emotion));
        end
    elseif dataset == 39
        edges = [1 2 3 4 5 6];
        for i_emotion = 1:length(remainingNodes)
            subplot(13,3,i_emotion);
            intensityRatings = histogram(subjectData(:,i_emotion),edges);
            ylim([0 150]);
            xlim([1 6]);
            title(remainingNodes(i_emotion));
        end
    end
    filename = ['subject' num2str(subjectID) '_histogram_intensity_ratings'];
    saveas(intensityRatings,filename,'tiff');

    %% histogram connection weights (correlations)
    halfSubjectMatrix = triu(subjectMatrix);
    halfSubjectMatrix = reshape(halfSubjectMatrix,[numel(halfSubjectMatrix),1]);
    halfSubjectMatrix = halfSubjectMatrix(halfSubjectMatrix~=0);
    figure;
    connectionWeights = histogram(halfSubjectMatrix);
    xlim([-1 1]);
    if numel(halfSubjectMatrix) < 153
        title(subjectID,'Color','r');
    else
        title(subjectID);
    end
    filename = ['subject' num2str(subjectID) '_histogram_connection_weights'];
    saveas(connectionWeights,filename,'tiff');
    
    %% count positive vs. negative links (connection weights/correlations)
    positiveLinks(i_subject) = sum(halfSubjectMatrix > 0);
    negativeLinks(i_subject) = sum(halfSubjectMatrix < 0);
    nullLinks(i_subject) = (numel(subjectMatrix)-length(subjectMatrix))/2-(positiveLinks(i_subject) + negativeLinks(i_subject));
    
    %% create heatmap
    figure;
    colormap('parula')
    heatmap = imagesc(subjectMatrix);
    colorbar;
    caxis([-1 1]);
    set(gca,'XTick',1:1:length(remainingNodes),'XTickLabel',remainingNodes);
    xtickangle(45)
    set(gca,'YTick',1:1:length(remainingNodes),'YTickLabel',remainingNodes);
    if numel(halfSubjectMatrix) < 153
        title(subjectID,'Color','r');
    else
        title(subjectID);
    end
    filename = ['subject' num2str(subjectID) '_heatmap_correlation_matrix'];
    saveas(heatmap,filename,'tiff');
    
    close all;
end

%% write file for positive vs. negative link counts
links = horzcat(subjectIDlist,positiveLinks',negativeLinks',nullLinks');
colNames = {'SubjectID','PositiveLinks','NegativeLinks','NullLinks'};
linksTable = array2table(links,'VariableNames',colNames);
writetable(linksTable,'links_counts.xlsx');
filename = ['dataSet' num2str(dataSet) '_links.mat'];
save(filename,'links');