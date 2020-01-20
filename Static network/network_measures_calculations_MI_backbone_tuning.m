% NOTE: Script is not configured to use gamma values other than the default of 1 for community detection

clear;
clc;

%% specify dataset and parameters
dataSet = 18; % set to 18, 39, or 88
print = 0;
test = 1; % set to 1 to correlate network measures with outcome variables

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
allData = rawData.data;
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
    % prune edges using backbone
    [nRow,nCol] = size(subjectMatrix);
%     testValues = 2:1:(nCol-1);
    testValues = [round(nCol-(3*nCol/4)) round(nCol-(2*nCol/3)) round(nCol-(nCol/2)) round(nCol-(nCol/3)) round(nCol-(nCol/4))];
    for i_value = 1:length(testValues)
        subjectMatrixTemp = round(subjectMatrix,4);
%         [~,subjectMatrixTemp] = backbone_wu(subjectMatrixTemp,testValues(i_value));
        % implement backbone function directly
        CIJ = subjectMatrixTemp;
        avgdeg = testValues(i_value);
        N = size(CIJ,1);
        CIJtree = zeros(N);
        % find strongest edge (note if multiple edges are tied, only use first one)
        [i,j,s] = find(max(max(CIJ))==CIJ);                     
        im = [i(1) i(2)];
        jm = [j(1) j(2)];
        % copy into tree graph
        CIJtree(im,jm) = CIJ(im,jm);
        in = im;
        out = setdiff(1:N,in);
        % repeat N-2 times
        for n=1:N-2
            % find strongest link between 'in' and 'out',ignore tied ranks
            [i,j,s] = find(max(max(CIJ(in,out)))==CIJ(in,out)); 
            im = in(i(1));
            jm = out(j(1));
            % copy into tree graph
            CIJtree(im,jm) = CIJ(im,jm); CIJtree(jm,im) = CIJ(jm,im);
            in = [in jm];                                    
            out = setdiff(1:N,in);
        end
        % now add connections back, with the total number of added connections determined by the desired 'avgdeg'
        CIJnotintree = CIJ.*~CIJtree;
        [a,b] = sort(nonzeros(CIJnotintree),'descend');
        cutoff = avgdeg*N - 2*(N-1);
        if cutoff>numel(a)
            continue
        end
        thr = a(cutoff);
        CIJclus = CIJtree + CIJnotintree.*(CIJnotintree>=thr);        
        subjectMatrixTemp = CIJclus;
        subjectMatrixNew = reshape(subjectMatrixTemp,nRow,nCol);
        subjectMatrixPruned(:,:,i_value) = subjectMatrixNew;
        clear subjectMatrixTemp subjectMatrixNew;
        % calculate clustering coefficient (original connection weights)
        [C,C_null,Ctot,Ctot_null] = clustering_coef_wu_sign(subjectMatrixPruned(:,:,i_value),3);
        clusterCoefGlobal(i_subject,i_value) = Ctot; % save off global value
        % calculate weighted network density (original connection weigthts)
        halfSubjectMatrix = triu(subjectMatrixPruned(:,:,i_value));
        halfSubjectMatrix = reshape(halfSubjectMatrix,[numel(halfSubjectMatrix),1]);
        halfSubjectMatrix = halfSubjectMatrix(halfSubjectMatrix~=0);
        kdenW(i_subject,i_value) = mean(halfSubjectMatrix); 
        % run community detection (original connection weights)
        gamma = 1;
        W = subjectMatrixPruned(:,:,i_value);  % set subject matrix
        B = 'modularity';
        n = size(W,1);      % number of nodes
        M = 1:n;            % initial community affiliations
        Q0 = -1; Q = 0;    % initialize modularity values
            while Q - Q0 > 1e-5    % while modularity increases
                Q0 = Q;            % perform community detection
                [M, Q] = community_louvain(W,gamma,[],B); 
            end
        communityAssignment = M; % save community assignment vector    
        numCommunities(i_subject,i_value) = max(M);   % save number of communities
        % calculate diversity coefficients (original connection weights)
        [Hpos Hneg] = diversity_coef_sign(W,M); % diversity coefficient
        dPos(i_subject,i_value) = mean(Hpos,1);
        % calculate community radius (original connection weights)
        locations = table2array(words(:,[3 4]));
        locations(missingColVar,:) = [];
        partitions = communityAssignment;
        commRadiusArray = comm_radius(partitions,locations);
        commRadius(i_subject,i_value) = commRadiusArray(end,end);
        % run community assignment (original connection weights)
        subjectM0 = assignedM; % copy overall community affiliation vector
        subjectM0(missingColVar) = []; % delete missing node from community affiliation vector
        M0 = subjectM0; % set initial community affiliations vector
        [~, Q_VA] = community_louvain(W,gamma,M0,B);
        modularity_VA(i_subject,i_value) = Q_VA; % save modularity value
    end
    clear subjectMatrixPruned
end

%% create summary table and write to file
clusterCoefGlobal(clusterCoefGlobal==0) = NaN;
clusterCoefGlobal = fillmissing(clusterCoefGlobal,'previous',2);

kdenW(kdenW==0) = NaN;
kdenW = fillmissing(kdenW,'previous',2);

numCommunities(numCommunities==0) = NaN;
numCommunities = fillmissing(numCommunities,'previous',2);

dPos(dPos==0) = NaN;
dPos = fillmissing(dPos,'previous',2);

commRadius(commRadius==0) = NaN;
commRadius = fillmissing(commRadius,'previous',2);

modularity_VA(modularity_VA==0) = NaN;
modularity_VA = fillmissing(modularity_VA,'previous',2);

for i_var = 1:size(numCommunities,2)
    networkMeasures(:,:,i_var) = horzcat(clusterCoefGlobal(:,i_var), kdenW(:,i_var), numCommunities(:,i_var), dPos(:,i_var), commRadius(:,i_var), modularity_VA(:,i_var));
end
% variableNames = {'cluster','density','numCom','dPos','comRadius','modVA'};
% save(['dataSet' num2str(dataSet) '_networkMeasures.mat'],'networkMeasures');
% networkMeasures_Table = array2table(networkMeasures,'VariableNames',variableNames);
% if print == 1
%     writetable(networkMeasures_Table,['dataSet' num2str(dataSet) '_networkMeasures.xlsx']);
% end

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
        for i_col = 1:size(networkMeasures,3)
            [r_GAD(:,i_col),p_GAD(:,i_col)] = corr(networkMeasures(:,:,i_col),GADtotal,'rows','complete');
            [r_PHQ8(:,i_col),p_PHQ8(:,i_col)] = corr(networkMeasures(:,:,i_col),PHQ8total,'rows','complete');
        end
        r_GAD_mean = mean(abs(r_GAD),1);
        r_PHQ8_mean = mean(abs(r_PHQ8),1);
        p_GAD_mean = mean(p_GAD,1);
        p_PHQ8_mean = mean(p_PHQ8,1);
        optValue_GAD = testValues(p_GAD_mean==min(p_GAD_mean));
        optValue_PHQ8 = testValues(p_PHQ8_mean==min(p_PHQ8_mean));
    elseif dataSet == 39
        questionnaires = importdata('39ESdata_QuestionnaireData.xlsx');
        BAItotal = questionnaires.data(:,9); 
        BDItotal = questionnaires.data(:,12);
        for i_col = 1:size(networkMeasures,3)
            [r_BAI(:,i_col),p_BAI(:,i_col)] = corr(networkMeasures(:,:,i_col),BAItotal,'rows','complete');
            [r_BDI(:,i_col),p_BDI(:,i_col)] = corr(networkMeasures(:,:,i_col),BDItotal,'rows','complete');
        end
        r_BAI_mean = mean(abs(r_BAI),1);
        r_BDI_mean = mean(abs(r_BDI),1);
        p_BAI_mean = mean(p_BAI,1);
        p_BDI_mean = mean(p_BDI,1);
        optValue_BAI = testValues(p_BAI_mean==min(p_BAI_mean));
        optValue_BDI = testValues(p_BDI_mean==min(p_BDI_mean));
    end
    for i_col = 1:size(networkMeasures,3)
        [r_zICC(:,i_col),p_zICC(:,i_col)] = corr(networkMeasures(:,:,i_col),zICC,'rows','complete');
        [r_neg(:,i_col),p_neg(:,i_col)] = corr(networkMeasures(:,:,i_col),neg,'rows','complete');
        [r_negSD(:,i_col),p_negSD(:,i_col)] = corr(networkMeasures(:,:,i_col),negSD,'rows','complete');
    end
    r_zICC_mean = mean(abs(r_zICC),1);
    r_neg_mean = mean(abs(r_neg),1);
    p_zICC_mean = mean(p_zICC,1);
    p_neg_mean = mean(p_neg,1);
    optValue_zICC = testValues(p_zICC_mean==min(p_zICC_mean));
    optValue_neg = testValues(p_neg_mean==min(p_neg_mean));
end