% This script generates, across all subjects in the data set, weighted
% and unweighted network density at each of 3 pre-specified connection weight thresholds,
% as well as histograms of the distribution at each threshold

clear;
clc;

%% specify dataset
dataSet = 18; % set to 18, 39, or 88
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
    
    %% calculate network density on original connection weights
    kden(i_subject,1) = density_und(subjectMatrix);
    halfSubjectMatrix = triu(subjectMatrix);
    halfSubjectMatrix = reshape(halfSubjectMatrix,[numel(halfSubjectMatrix),1]);
    halfSubjectMatrix = halfSubjectMatrix(halfSubjectMatrix~=0);
    kdenW(i_subject,1) = mean(halfSubjectMatrix);
    
    %% network density on thresholded connection weights, first test value
    threshold = 0.1; % removes correlations below given absolute value
    subjectMatrixT1 = subjectMatrix;
    for i_row = 1:length(subjectMatrixT1)
        for i_col = 1:length(subjectMatrixT1)
            if abs(subjectMatrixT1(i_row,i_col)) < threshold
                subjectMatrixT1(i_row,i_col) = 0;
            end
        end
    end
    kden(i_subject,2) = density_und(subjectMatrixT1);
    halfSubjectMatrixT1 = triu(subjectMatrixT1);
    halfSubjectMatrixT1 = reshape(halfSubjectMatrixT1,[numel(halfSubjectMatrixT1),1]);
    halfSubjectMatrixT1 = halfSubjectMatrixT1(halfSubjectMatrixT1~=0);
    kdenW(i_subject,2) = mean(halfSubjectMatrixT1);

    %% network density on thresholded connection weights, second test value
    threshold = 0.2; % removes correlations below given absolute value
    subjectMatrixT2 = subjectMatrix;
    for i_row = 1:length(subjectMatrixT2)
        for i_col = 1:length(subjectMatrixT2)
            if abs(subjectMatrixT2(i_row,i_col)) < threshold
                subjectMatrixT2(i_row,i_col) = 0;
            end
        end
    end
    kden(i_subject,3) = density_und(subjectMatrixT2);
    halfSubjectMatrixT2 = triu(subjectMatrixT2);
    halfSubjectMatrixT2 = reshape(halfSubjectMatrixT2,[numel(halfSubjectMatrixT2),1]);
    halfSubjectMatrixT2 = halfSubjectMatrixT2(halfSubjectMatrixT2~=0);
    kdenW(i_subject,3) = mean(halfSubjectMatrixT2);
end

%% write results for network density and distribution across subjects
networkDensity = horzcat(subjectIDlist,kden,kdenW);
colNames = {'SubjectID','OriginalData','Threshold1','Threshold2','OriginalWeighted','T1Weighted','T2Weighted'};
densityTable = array2table(networkDensity,'VariableNames',colNames);
writetable(densityTable,'network_density_distribution.xlsx');
filename = ['dataSet' num2str(dataSet) '_density.mat'];
save(filename,'networkDensity');

figure;
histogram1 = histogram(networkDensity(:,2));
xlim([0 max(networkDensity(:,4))+1]);
ylim([1 length(subjectIDlist)]);
title('Network density in original data');
saveas(histogram1,'Network_density_distribution_original','tiff');

figure;
histogram2 = histogram(networkDensity(:,3));
xlim([0 max(networkDensity(:,4))+1]);
ylim([1 length(subjectIDlist)]);
title('Network density at first threshold');
saveas(histogram2,'Network_density_distribution_threshold1','tiff');

figure;
histogram3 = histogram(networkDensity(:,4));
xlim([0 max(networkDensity(:,4))+1]);
ylim([1 length(subjectIDlist)]);
title('Network density at second threshold');
saveas(histogram3,'Network_density_distribution_threshold2','tiff');

figure;
histogram4 = histogram(networkDensity(:,5));
xlim([0 max(networkDensity(:,7))]);
ylim([1 length(subjectIDlist)]);
title('Weighted network density in original data');
saveas(histogram4,'Weighted_network_density_distribution_original','tiff');

figure;
histogram5 = histogram(networkDensity(:,6));
xlim([0 max(networkDensity(:,7))]);
ylim([1 length(subjectIDlist)]);
title('Weighted network density at first threshold');
saveas(histogram5,'Weighted_network_density_distribution_threshold1','tiff');

figure;
histogram6 = histogram(networkDensity(:,7));
xlim([0 max(networkDensity(:,7))]);
ylim([1 length(subjectIDlist)]);
title('Weighted network density at second threshold');
saveas(histogram6,'Weighted_network_density_distribution_threshold2','tiff');

close all;

