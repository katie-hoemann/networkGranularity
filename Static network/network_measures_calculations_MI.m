% NOTE: Script is not configured to use gamma values other than the default of 1 for community detection

clear;
clc;

%% specify dataset and parameters
dataSet = 18; % set to 18, 39, or 88
print = 0; 
test = 0; % set to 1 to correlate network measures with outcome variables

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
    % calculate mutual information matrix
    rows = 1:1:size(subjectData,2);
    for i_var = 1:length(rows)
        for j_var = rows(rows~=i_var)
            subjectMatrix(i_var,j_var) = nmi(subjectData(:,i_var),subjectData(:,j_var));
        end
    end
    % delete zero-variance rows/columns from data
    colVar = var(subjectMatrix,0,1); % find columns (nodes) with zero variance
    rowVar = var(subjectMatrix,0,2); % find rows (instances) with zero variance
    missingColVar = find(colVar==0); % index zero-variance nodes
    missingRowVar = find(rowVar==0); % index zero-variance instances
    subjectMatrix(:,missingColVar) = []; % remove zero-variance nodes
    subjectMatrix(missingRowVar,:) = []; % remove zero-variance instances
    % track missing/removed nodes
    deletedNodes{i_subject} = wordList(missingColVar); % note which nodes are not in the network
    remainingNodes = wordList;
    remainingNodes(missingColVar) = [];
    numNodes(i_subject) = size(remainingNodes,1); % record number of remaining nodes
    % calculate clustering coefficient (original connection weights)
    [C,C_null,Ctot,Ctot_null] = clustering_coef_wu_sign(subjectMatrix,3);
    clusterCoefGlobal(i_subject) = Ctot; % save off global value
    clusterCoefGlobalSD(i_subject) = std(C); % save off SD of node-level measure
    % calculate weighted network density (original connection weigthts)
    kden(i_subject) = density_und(subjectMatrix);
    halfSubjectMatrix = triu(subjectMatrix);
    halfSubjectMatrix = reshape(halfSubjectMatrix,[numel(halfSubjectMatrix),1]);
    halfSubjectMatrix = halfSubjectMatrix(halfSubjectMatrix~=0);
    kdenW(i_subject) = mean(halfSubjectMatrix); 
    % run community detection (original connection weights)
    gamma = 1;
    W = subjectMatrix;  % set subject matrix
    B = 'modularity';
    n = size(W,1);      % number of nodes
    M = 1:n;            % initial community affiliations
    Q0 = -1; Q = 0;    % initialize modularity values
        while Q - Q0 > 1e-5    % while modularity increases
            Q0 = Q;            % perform community detection
            [M, Q] = community_louvain(W,gamma,[],B); 
        end
    communityAssignment = M; % save community assignment vector    
    numCommunities(i_subject) = max(M);   % save number of communities
    modularity(i_subject) = Q; % save modularity value
    % calculate participation & diversity coefficients (original connection weights)
    [Ppos Pneg] = participation_coef_sign(W,M); % participation coefficient
    participation(i_subject,1) = mean(Ppos,1);
    participation(i_subject,2) = mean(Pneg,1);
    participationSD(i_subject,1) = std(Ppos,1);
    participationSD(i_subject,2) = std(Pneg,1);
    [Hpos Hneg] = diversity_coef_sign(W,M); % diversity coefficient
    diversity(i_subject,1) = mean(Hpos,1);
    diversity(i_subject,2) = mean(Hneg,1);
    diversitySD(i_subject,1) = std(Hpos,1);
    diversitySD(i_subject,2) = std(Hneg,1);
    % calculate community radius (original connection weights)
    locations = table2array(words(:,[3 4]));
    locations(missingColVar,:) = [];
    partitions = communityAssignment;
    commRadiusArray = comm_radius(partitions,locations);
    commRadius(:,i_subject) = commRadiusArray(end,end);
    % run community assignment (original connection weights)
    subjectM0 = assignedM; % copy overall community affiliation vector
    subjectM0(missingColVar) = []; % delete missing nodes from community affiliation vector
    M0 = subjectM0; % set initial community affiliations vector
    [~, Q_VA] = community_louvain(W,gamma,M0,B);
    modularity_VA(i_subject) = Q_VA; % save modularity value
    clear subjectMatrix missingColVar
end

%% create summary table and write to file
networkMeasures = horzcat(clusterCoefGlobal', clusterCoefGlobalSD', kdenW', numCommunities', modularity', participation, participationSD, diversity, diversitySD, commRadius', modularity_VA');
variableNames = {'cluster','clusterSD','density','numCom','mod','pPos','pNeg','pPosSD','pNegSD','dPos','dNeg','dPosSD','dNegSD','comRadius','modVA'};
save(['dataSet' num2str(dataSet) '_networkMeasures.mat'],'networkMeasures');
networkMeasures_Table = array2table(networkMeasures,'VariableNames',variableNames);
if print == 1
    writetable(networkMeasures_Table,['dataSet' num2str(dataSet) '_networkMeasures.xlsx']);
end

%% run Kolmogorov-Smirnov test to assess normality; normalize distributions as needed
[~,numVariables] = size(networkMeasures);
for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(networkMeasures(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(networkMeasures(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm(:,i_variable) = norminv(fracRank,mean(networkMeasures(:,i_variable)),std(networkMeasures(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm(:,i_variable) = networkMeasures(:,i_variable);
    end
end

%% compile all (normalized) data and re-export
save(['dataSet' num2str(dataSet) '_networkMeasures_normalized.mat'],'fracRankNorm');
fracRankNorm_Table = array2table(fracRankNorm,'VariableNames',variableNames);
if print == 1
    writetable(fracRankNorm_Table,['dataSet' num2str(dataSet) '_networkMeasures_normalized.xlsx']);
end

%% test against zICC, anxiety, depression
if test == 1
    load(['dataSet' num2str(dataSet) '_zICC.mat']);
    zICC = rtoZ(:,3);
    load(['dataSet' num2str(dataSet) '_affect.mat']);
    neg = affect(:,2); % mean negative affect
    negSD = affect(:,4); % std negative affect
    pos = affect(:,1); % mean positive affect
    posSD = affect(:,3); % std positive affect
    if dataSet == 18
        questionnaires = importdata('18ARIVariations_QuestionnaireData.xlsx');
        GADtotal = questionnaires.data(:,12); 
        PHQ8total = questionnaires.data(:,23); %23 for session 1; 48 for session 2
        [r_GAD,p_GAD] = corr(networkMeasures,GADtotal,'rows','complete');
        [r_PHQ8,p_PHQ8] = corr(networkMeasures,PHQ8total,'rows','complete');
    elseif dataSet == 39
        questionnaires = importdata('39ESdata_QuestionnaireData.xlsx');
        BAItotal = questionnaires.data(:,9); 
        BDItotal = questionnaires.data(:,12);
        [r_BAI,p_BAI] = corr(networkMeasures,BAItotal,'rows','complete');
        [r_BDI,p_BDI] = corr(networkMeasures,BDItotal,'rows','complete');
    end
    [r_zICC,p_zICC] = corr(networkMeasures,zICC,'rows','complete');
    [r_neg,p_neg] = corr(networkMeasures,neg,'rows','complete');
    [r_negSD,p_negSD] = corr(networkMeasures,negSD,'rows','complete');
end
