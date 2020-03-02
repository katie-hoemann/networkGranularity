% For every subject in the specified data set, this script assigns an
% initial community affiliation vector based on valence (i.e., 2 communities) 
% and determines modularity values at 3 pre-specified connection weight thresholds
% Across all subjects in the data set, this script generates a spreadsheet
% of the modularity at each threshold, as well as histograms of the distributions

clear;
clc;

%% specify dataset
dataSet = 18; % set to 18, 39, or 88
polychoric = 0; % set to 0 to use standard Pearson correlation matrix
defaultGamma = 1; % set to 1 to use default gamma value of 1
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

%% determine community affiliation based on valence of sampled words
for i_word = 1:height(words) % define valence categories
if words.Valence(i_word) > 5 % derived based on the database mean for Warriner et al (2013)
    assignedM(i_word) = 1;
else
    assignedM(i_word) = 2;
end
end 

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
        subjectM0 = assignedM; % copy overall community affiliation vector
        subjectM0(missingColVar) = []; % delete missing node from community affiliation vector
    else
        subjectMatrix = corrcoef(subjectData); % calculate n-by-n correlation matrix
        %% delete NaN values from correlation matrix
        missingCorr = find(isnan(subjectMatrix(1,:)));
        deletedCorr = isnan(subjectMatrix(1,:)); % save off which row/col were removed
        deletedNodes = wordList(deletedCorr==1); % note which nodes are not in the network
        remainingNodes = wordList(~deletedCorr==1);
        subjectMatrix(missingCorr,:) = []; % delete missing row
        subjectMatrix(:,missingCorr) = []; % delete missing column
        subjectM0 = assignedM; % copy overall community affiliation vector
        subjectM0(missingCorr) = []; % delete missing node from community affiliation vector
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
    M0 = subjectM0; % set initial community affiliations vector
    B = 'negative_asym';
    n = size(W,1);      % number of nodes
    M = 1:n;            % initial community affiliations
%     Q0 = -1; Q = 0;    % initialize modularity values
    [M, Q] = community_louvain(W,gamma,M0,B);
%         while Q - Q0 > 1e-5    % while modularity increases
%             Q0 = Q;            % perform community detection
%             [M, Q] = community_louvain(W,gamma,M0,B); 
%         end
    communityAssignment_VA(:,1) = M; % save community assignment vector    
    numCommunities_VA(i_subject,1) = max(M);   % save number of communities
    modularity_VA(i_subject,1) = Q; % save modularity value

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
    M0 = subjectM0; % set initial community affiliations vector
    B = 'negative_asym';
    n = size(W,1);      % number of nodes
    M = 1:n;            % initial community affiliations
%     Q0 = -1; Q = 0;    % initialize modularity values
    [M, Q] = community_louvain(W,gamma,M0,B);
%         while Q - Q0 > 1e-5    % while modularity increases
%             Q0 = Q;            % perform community detection
%             [M, Q] = community_louvain(W,gamma,M0,B); 
%         end
    communityAssignment_VA(:,2) = M; % save community assignment vector       
    numCommunities_VA(i_subject,2) = max(M);   % save number of communities
    modularity_VA(i_subject,2) = Q; % save modularity value

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
    M0 = subjectM0; % set initial community affiliations vector
    B = 'negative_asym';
    n = size(W,1);      % number of nodes
    M = 1:n;            % initial community affiliations
%     Q0 = -1; Q = 0;    % initialize modularity values
    [M, Q] = community_louvain(W,gamma,M0,B);
%         while Q - Q0 > 1e-5    % while modularity increases
%             Q0 = Q;            % perform community detection
%             [M, Q] = community_louvain(W,gamma,M0,B); 
%         end
    communityAssignment_VA(:,3) = M; % save community assignment vector       
    numCommunities_VA(i_subject,3) = max(M);   % save number of communities
    modularity_VA(i_subject,3) = Q; % save modularity value

    %% write results for community assignment per subject
    if subjectFiles == 1
        filename = ['subject' num2str(subjectID) '_community_assignment_VA.xlsx'];
        communities_VA = array2table(communityAssignment_VA);
        rowNames = cell2table(remainingNodes);
        communities_VA = horzcat(rowNames,communities_VA);
        communities_VA.Properties.VariableNames = {'Emotion','OriginalData','Threshold1','Threshold2'};
        communities_VA = sortrows(communities_VA,2);
        writetable(communities_VA,filename);
    end
    
    %% clear variables
    clear communityAssignment_VA;
end

%% write results for number of communities and distribution across subjects
numCommunities_VA = horzcat(subjectIDlist,numCommunities_VA);
colNames = {'SubjectID','OriginalData','Threshold1','Threshold2'};
numCommunitiesVAtable = array2table(numCommunities_VA,'VariableNames',colNames);
writetable(numCommunitiesVAtable,'number_communities_distribution_VA.xlsx');
filename = ['dataSet' num2str(dataSet) '_communities_VA.mat'];
save(filename,'numCommunities_VA');

figure;
histogram1 = histogram(numCommunities_VA(:,2));
xlim([1 max(numCommunities_VA(:,4))+1]);
ylim([1 length(subjectIDlist)]);
title('Number of communities in original data');
saveas(histogram1,'Community_distribution_original_VA','tiff');

figure;
histogram2 = histogram(numCommunities_VA(:,3));
xlim([1 max(numCommunities_VA(:,4))+1]);
ylim([1 length(subjectIDlist)]);
title('Number of communities at first threshold');
saveas(histogram2,'Community_distribution_threshold1_VA','tiff');

figure;
histogram3 = histogram(numCommunities_VA(:,4));
xlim([1 max(numCommunities_VA(:,4))+1]);
ylim([1 length(subjectIDlist)]);
title('Number of communities at second threshold');
saveas(histogram3,'Community_distribution_threshold2_VA','tiff');

close all;

%% write results for modularity and distribution across subjects
modularity_VA = horzcat(subjectIDlist,modularity_VA);
colNames = {'SubjectID','OriginalData','Threshold1','Threshold2'};
modularityVAtable = array2table(modularity_VA,'VariableNames',colNames);
writetable(modularityVAtable,'modularity_distribution_VA.xlsx');
filename = ['dataSet' num2str(dataSet) '_modularity_VA.mat'];
save(filename,'modularity_VA');

figure;
histogram4 = histogram(modularity_VA(:,2));
xlim([-1 1]);
ylim([1 length(subjectIDlist)]);
title('Modularity in original data');
saveas(histogram4,'Modularity_distribution_original_VA','tiff');

figure;
histogram5 = histogram(modularity_VA(:,3));
xlim([-1 1]);
ylim([1 length(subjectIDlist)]);
title('Modularity at first threshold');
saveas(histogram5,'Modularity_distribution_threshold1_VA','tiff');

figure;
histogram6 = histogram(modularity_VA(:,4));
xlim([-1 1]);
ylim([1 length(subjectIDlist)]);
title('Modularity at second threshold');
saveas(histogram6,'Modularity_distribution_threshold2_VA','tiff');

close all;