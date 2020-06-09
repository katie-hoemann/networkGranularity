% ISSUE: Update column IDs on second pass for zero-variance 

clear;
clc;

%% specify dataset and parameters
dataSet = 39; % set to 18, 39, or 88
print = 1; 

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
questionnaires = importdata(dataFile);
allData = questionnaires.data;
subjectIDlist = unique(questionnaires.data(:,1)); % grab subject IDs from first column of data file
words = readtable(wordFile); 
wordList = questionnaires.colheaders(2:end)';  % grab sampled words from top row of data file

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

%% set valence and arousal categories for sampled words
for i_word = 1:height(words) % define valence categories
if words.Valence(i_word) > 5 % derived based on the database mean for Warriner et al (2013)
    valCat(i_word) = {'Positive'};
    positive(i_word) = 1;
    valence(i_word) = 1;
else
    valCat(i_word) = {'Negative'};
    positive(i_word) = 0;
    valence(i_word) = 2;
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

%% determine community affiliation based on valence of sampled words
for i_word = 1:height(words) % define valence categories
if words.Valence(i_word) > 5 % derived based on the database mean for Warriner et al (2013)
    assignedM(i_word) = 1;
else
    assignedM(i_word) = 2;
end
end 

%% grab data for each subject and run through calculations
for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index1 = find(allData(:,1)==subjectID);
    subjectData = allData(index1,2:end);
    % remove invalid values and missing data
    subjectData(subjectData > maximumRating) = NaN;
    missingData = isnan(subjectData); 
    missingData2 = any(missingData,2); 
    subjectData = subjectData(~missingData2,:);
    if numel(subjectData) == 0
        continue
    end
    % if necessary, rescale data to start at 0
    if startRatingsat1 == 1
        subjectData = subjectData-1;
    end
    % delete zero-variance rows/columns from data
    colVar1 = var(subjectData,0,1); % find columns (nodes) with zero variance
    rowVar1 = var(subjectData,0,2); % find rows (instances) with zero variance
    missingColVar1 = find(colVar1==0); % index zero-variance nodes
    missingRowVar1 = find(rowVar1==0); % index zero-variance instances
    subjectData(:,missingColVar1) = []; % remove zero-variance nodes
    subjectData(missingRowVar1,:) = []; % remove zero-variance instances
    % repeat in case removal resulted in new zero-variance rows/columns
    colVar2 = var(subjectData,0,1); % find columns (nodes) with zero variance
    rowVar2 = var(subjectData,0,2); % find rows (instances) with zero variance
    missingColVar2 = find(colVar2==0); % index zero-variance nodes
    missingRowVar2 = find(rowVar2==0); % index zero-variance instances
    subjectData(:,missingColVar2) = []; % remove zero-variance nodes
    subjectData(missingRowVar2,:) = []; % remove zero-variance instances
    % track missing/removed nodes
    missingColVarAll = [missingColVar1 missingColVar2]; % one vector with all missing nodes
    deletedNodes{i_subject} = wordList(missingColVarAll); % note which nodes are not in the network
    remainingNodes = wordList;
    remainingNodes(missingColVarAll) = [];
    numNodes(i_subject) = size(remainingNodes,1); % record number of remaining nodes
    % calculate correlation matrix, iterating 100x to get mean values
    for i_iter = 1:100
        subjectMatrix_Iter(:,:,i_iter) = polychoric_proc_missing(subjectData,NaN); % calculate n-by-n polychoric correlation matrix
    end
    subjectMatrix = mean(subjectMatrix_Iter,3,'omitnan');
    subjectMatrix(isnan(subjectMatrix)) = 0; % set NaN values to 0
    subjectMatrix(logical(eye(size(subjectMatrix)))) = 0; % set diagonal to 0 for BCT functions
    % calculate clustering coefficient
    [C,C_null,Ctot,Ctot_null] = clustering_coef_wu_sign(subjectMatrix,3);
    clusterCoefGlobal(i_subject) = Ctot; % save off global value
    clusterCoefGlobalSD(i_subject) = std(C); % save off SD of node-level measure
    % calculate weighted network density
    kden(i_subject) = density_und(subjectMatrix);
    halfSubjectMatrix = triu(subjectMatrix);
    halfSubjectMatrix = reshape(halfSubjectMatrix,[numel(halfSubjectMatrix),1]);
    halfSubjectMatrix = halfSubjectMatrix(halfSubjectMatrix~=0);
    kdenW(i_subject) = mean(halfSubjectMatrix); 
    % run community detection, iterating 1000x to maximize Q
    for i_iter = 1:1000
        gamma = 1;
        W = subjectMatrix;  % set subject matrix
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
    communityAssignment = communityAssignment_Iter(:,maxModularity_Index); % save community assignment vector    
    numCommunities(i_subject) = max(communityAssignment);   % save number of communities
    modularity(i_subject) = modularity_Iter(maxModularity_Index); % save modularity value
    % calculate participation & diversity coefficients
    [Ppos Pneg] = participation_coef_sign(subjectMatrix,communityAssignment); % participation coefficient
    participation(i_subject,1) = mean(Ppos,1);
    participation(i_subject,2) = mean(Pneg,1);
    participationSD(i_subject,1) = std(Ppos,1);
    participationSD(i_subject,2) = std(Pneg,1);
    [Hpos Hneg] = diversity_coef_sign(subjectMatrix,communityAssignment); % diversity coefficient
    diversity(i_subject,1) = mean(Hpos,1);
    diversity(i_subject,2) = mean(Hneg,1);
    diversitySD(i_subject,1) = std(Hpos,1);
    diversitySD(i_subject,2) = std(Hneg,1);
    % calculate community radius
    locations = table2array(words(:,[3 4]));
    locations(missingColVarAll,:) = [];
    partitions = communityAssignment;
    comRadiusArray = comm_radius(partitions,locations);
    comRadius(:,i_subject) = comRadiusArray(end,end);
    % run community assignment
    subjectM0 = assignedM; % copy overall community affiliation vector
    subjectM0(missingColVarAll) = []; % delete missing nodes from community affiliation vector
    M0 = subjectM0; % set initial community affiliations vector
    [~, Q_VA] = community_louvain(W,gamma,M0,B);
    modularity_VA(i_subject) = Q_VA; % save modularity value
    fprintf('participant %d\n',subjectID)
    clearvars subjectMatrix_Iter subjectMatrix communityAssignment_Iter modularity_Iter communityAssignment
end

%% create summary table and write to file
networkMeasures = horzcat(clusterCoefGlobal', clusterCoefGlobalSD', kdenW', numCommunities', modularity', participation, participationSD, diversity, diversitySD, comRadius', modularity_VA');
variableNames = {'cluster','clusterSD','density','numCom','mod','pPos','pNeg','pPosSD','pNegSD','dPos','dNeg','dPosSD','dNegSD','comRadius','modVA'};
save(['dataSet' num2str(dataSet) '_networkMeasures_polychoric.mat'],'networkMeasures');
networkMeasures_Table = array2table(networkMeasures,'VariableNames',variableNames);
if print == 1
    writetable(networkMeasures_Table,['dataSet' num2str(dataSet) '_networkMeasures_polychoric.xlsx']);
end

%% run Kolmogorov-Smirnov test to assess normality; normalize distributions as needed
[~,numVariables] = size(networkMeasures);
for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(networkMeasures(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(networkMeasures(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm(:,i_variable) = norminv(fracRank,mean(networkMeasures(:,i_variable),'omitnan'),std(networkMeasures(:,i_variable),'omitnan'));
    elseif h(i_variable) == 0
        fracRankNorm(:,i_variable) = networkMeasures(:,i_variable);
    end
end

%% compile all (normalized) data and re-export
save(['dataSet' num2str(dataSet) '_networkMeasures_polychoric_normalized.mat'],'fracRankNorm');
fracRankNorm_Table = array2table(fracRankNorm,'VariableNames',variableNames);
if print == 1
    writetable(fracRankNorm_Table,['dataSet' num2str(dataSet) '_networkMeasures_polychoric_normalized.xlsx']);
end