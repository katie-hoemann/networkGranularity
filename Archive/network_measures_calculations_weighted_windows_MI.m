% NOTES: Currently built out for 3-day weighted sliding window; modification to
% window size or weighting will require structural edits to script
% Script is not configured to identify deleted/remaining nodes in each window, or to
% use gamma values other than the default of 1 for community detection

clear;
clc;

%% specify dataset and parameters
dataSet = 39; % only 18 or 39 for now
printSummaries = 1; % set to 1 to save summary tables to file
printTables = 0; % set to 1 to write measure-specific tables to file
saveData = 0; % set to 1 to save measure-specific data matrices
test = 0; % set to 1 to correlate network measures with outcome variables

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
% else
%     dataFile = '88PANASdata.xlsx';
%     wordFile = 'words88.csv';
end

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
    dayIDlist = unique(subjectData(:,1));
    for i_day = 1:length(dayIDlist)
        dayIDlistUpdate = dayIDlist;
        dayData = [];
        dayID = dayIDlist(i_day);
        index2 = find(subjectData(:,1)==dayID);
        dayData = subjectData(index2,4:end);
        % remove invalid values and missing data
        dayData(dayData > maximumRating) = NaN;
        missingData = isnan(dayData); 
        missingData2 = any(missingData,2); 
        dayData = dayData(~missingData2,:); 
        if size(dayData,1) < 2
            dayIDlistUpdate(i_day) = [];
            continue
        end
        % if necessary, rescale data to start at 0
        if startRatingsat1 == 1
            dayData = dayData-1;
        end
        % calculate mutual information matrix
        rows = 1:1:size(dayData,2);
        for i_var = 1:length(rows)
            for j_var = rows(rows~=i_var)
                dayMatrix(i_var,j_var) = nmi(dayData(:,i_var),dayData(:,j_var));
            end
        end
        % delete zero-variance rows/columns from data
        colVar = var(dayMatrix,0,1); % find columns (nodes) with zero variance
        rowVar = var(dayMatrix,0,2); % find rows (instances) with zero variance
        missingColVar = find(colVar==0); % index zero-variance nodes
        missingRowVar = find(rowVar==0); % index zero-variance instances
        dayMatrix(:,missingColVar) = []; % remove zero-variance nodes
        dayMatrix(missingRowVar,:) = []; % remove zero-variance instances
        if numel(dayMatrix) == 0
            dayIDlistUpdate(i_day) = [];
            continue
        end
        % calculate clustering coefficient (original connection weights)
        [C,C_null,Ctot,Ctot_null] = clustering_coef_wu_sign(dayMatrix,3);
        clusterCoefGlobal(i_day) = Ctot; % save off global value
        clusterCoefGlobalSD(i_day) = std(C); % save off SD of node-level measure
        % calculate weighted network density (original connection weigthts)
        kden(i_day) = density_und(dayMatrix);
        halfdayMatrix = triu(dayMatrix);
        halfdayMatrix = reshape(halfdayMatrix,[numel(halfdayMatrix),1]);
        halfdayMatrix = halfdayMatrix(halfdayMatrix~=0);
        kdenW(i_day) = mean(halfdayMatrix); 
        % run community detection (original connection weights)
        gamma = 1;
        W = dayMatrix;  % set subject matrix
        B = 'modularity';
        n = size(W,1);      % number of nodes
        M = 1:n;            % initial community affiliations
        Q0 = -1; Q = 0;    % initialize modularity values
            while Q - Q0 > 1e-5    % while modularity increases
                Q0 = Q;            % perform community detection
                [M, Q] = community_louvain(W,gamma,[],B); 
            end
        communityAssignment = M; % save community assignment vector    
        numCommunities(i_day) = max(M);   % save number of communities
        modularity(i_day) = Q; % save modularity value
        % calculate participation & diversity coefficients (original connection weights)
        [Ppos Pneg] = participation_coef_sign(W,M); % participation coefficient
        participation(i_day,1) = mean(Ppos,1);
        participation(i_day,2) = mean(Pneg,1);
        participationSD(i_day,1) = std(Ppos,1);
        participationSD(i_day,2) = std(Pneg,1);
        [Hpos Hneg] = diversity_coef_sign(W,M); % diversity coefficient
        diversity(i_day,1) = mean(Hpos,1);
        diversity(i_day,2) = mean(Hneg,1);
        diversitySD(i_day,1) = std(Hpos,1);
        diversitySD(i_day,2) = std(Hneg,1);
        % calculate community radius (original connection weights)
        locations = table2array(words(:,[3 4]));
        locations(missingColVar,:) = [];
        partitions = communityAssignment;
        commRadiusArray = comm_radius(partitions,locations);
        commRadius(:,i_day) = commRadiusArray(end,end);
        % run community assignment (original connection weights)
        dayM0 = assignedM; % copy overall community affiliation vector
        dayM0(missingColVar) = []; % delete missing node from community affiliation vector
        M0 = dayM0; % set initial community affiliations vector
        [~, Q_VA] = community_louvain(W,gamma,M0,B);
        modularity_VA(i_day) = Q_VA; % save modularity value
        fprintf('participant %d, day %d\n',subjectID,i_day)
        clearvars dayMatrix communityAssignment missingColVar
    end
    %% combine measures
    networkMeasures(:,1) = clusterCoefGlobal;
    networkMeasures(:,2) = clusterCoefGlobalSD;
    networkMeasures(:,3) = kdenW;
    networkMeasures(:,4) = numCommunities;
    networkMeasures(:,5) = modularity;
    networkMeasures(:,6) = participation(:,1);
    networkMeasures(:,7) = participation(:,2);
    networkMeasures(:,8) = participationSD(:,1);
    networkMeasures(:,9) = participationSD(:,2);
    networkMeasures(:,10) = diversity(:,1);
    networkMeasures(:,11) = diversity(:,2);
    networkMeasures(:,12) = diversitySD(:,1);
    networkMeasures(:,13) = diversitySD(:,2);
    networkMeasures(:,14) = commRadius;
    networkMeasures(:,15) = modularity_VA;
    %% create table of values per day
    ppID = ['PP' num2str(subjectID)];
    day = (1:1:max(allData(:,2)))';
    day = array2table(day,'VariableNames',{'Day'});
    networkMeasures = horzcat(dayIDlistUpdate,networkMeasures);
    networkMeasures = array2table(networkMeasures,'VariableNames',{'ppDay','cluster','clusterSD','density','numCom','mod','pPos','pNeg','pPosSD','pNegSD','dPos','dNeg','dPosSD','dNegSD','comRadius','modVA'});
    networkMeasures = outerjoin(day,networkMeasures,'LeftKeys','Day','RightKeys','ppDay');
    networkMeasures = removevars(networkMeasures,'ppDay');
    %% impute any missing data
    networkMeasures = fillmissing(networkMeasures,'movmean',3); 
    %% create windowed values
    networkMeasures = table2array(networkMeasures(:,2:end));
    day = table2array(day);
    for i_window = 1:length(day)-2
        focusDay = i_window+1;
        windowData = [];
        windowData = [networkMeasures(focusDay-1,:); networkMeasures(focusDay,:); networkMeasures(focusDay+1,:)];
        windowSimpleMeans(i_window,:) = mean(windowData,1); % simple average across days in window
        windowWeightedMeans(i_window,:) = networkMeasures(focusDay-1,:)*1/6 + networkMeasures(focusDay,:)*4/6 + networkMeasures(focusDay+1,:)*1/6; % weighted average across days in window
    end 
    %% calculate subject-level mean and SD for measures
    networkMeasureSimpleMeans(i_subject,:) = mean(windowSimpleMeans,'omitnan');
    networkMeasureSimpleSDs(i_subject,:) = std(windowSimpleMeans,'omitnan');
    networkMeasureWeightedMeans(i_subject,:) = mean(windowWeightedMeans,'omitnan');
    networkMeasureWeightedSDs(i_subject,:) = std(windowWeightedMeans,'omitnan');
    %% add subject results to summary tables
    windowSimpleMeans_Table = array2table(windowSimpleMeans);
    windowWeightedMeans_Table = array2table(windowWeightedMeans);
    if i_subject == 1
        window = (1:1:max(allData(:,2))-2)';
        window = array2table(window);
     
        clustering_Table_S = horzcat(window,windowSimpleMeans_Table(:,1));
        clustering_Table_S.Properties.VariableNames = {'Window',ppID};
        clusteringSD_Table_S = horzcat(window,windowSimpleMeans_Table(:,2));
        clusteringSD_Table_S.Properties.VariableNames = {'Window',ppID};
        density_Table_S = horzcat(window,windowSimpleMeans_Table(:,3));
        density_Table_S.Properties.VariableNames = {'Window',ppID};
        numCom_Table_S = horzcat(window,windowSimpleMeans_Table(:,4));
        numCom_Table_S.Properties.VariableNames = {'Window',ppID};
        modularity_Table_S = horzcat(window,windowSimpleMeans_Table(:,5));
        modularity_Table_S.Properties.VariableNames = {'Window',ppID};
        participationPos_Table_S = horzcat(window,windowSimpleMeans_Table(:,6));
        participationPos_Table_S.Properties.VariableNames = {'Window',ppID};
        participationNeg_Table_S = horzcat(window,windowSimpleMeans_Table(:,7));
        participationNeg_Table_S.Properties.VariableNames = {'Window',ppID};
        participationPosSD_Table_S = horzcat(window,windowSimpleMeans_Table(:,8));
        participationPosSD_Table_S.Properties.VariableNames = {'Window',ppID};
        participationNegSD_Table_S = horzcat(window,windowSimpleMeans_Table(:,9));
        participationNegSD_Table_S.Properties.VariableNames = {'Window',ppID};
        diversityPos_Table_S = horzcat(window,windowSimpleMeans_Table(:,10));
        diversityPos_Table_S.Properties.VariableNames = {'Window',ppID};
        diversityNeg_Table_S = horzcat(window,windowSimpleMeans_Table(:,11));
        diversityNeg_Table_S.Properties.VariableNames = {'Window',ppID};
        diversityPosSD_Table_S = horzcat(window,windowSimpleMeans_Table(:,12));
        diversityPosSD_Table_S.Properties.VariableNames = {'Window',ppID};
        diversityNegSD_Table_S = horzcat(window,windowSimpleMeans_Table(:,13));
        diversityNegSD_Table_S.Properties.VariableNames = {'Window',ppID};
        comRadius_Table_S = horzcat(window,windowSimpleMeans_Table(:,14));
        comRadius_Table_S.Properties.VariableNames = {'Window',ppID};
        modVA_Table_S = horzcat(window,windowSimpleMeans_Table(:,15));
        modVA_Table_S.Properties.VariableNames = {'Window',ppID};

        clustering_Table_W = horzcat(window,windowWeightedMeans_Table(:,1));
        clustering_Table_W.Properties.VariableNames = {'Window',ppID};
        clusteringSD_Table_W = horzcat(window,windowWeightedMeans_Table(:,2));
        clusteringSD_Table_W.Properties.VariableNames = {'Window',ppID};
        density_Table_W = horzcat(window,windowWeightedMeans_Table(:,3));
        density_Table_W.Properties.VariableNames = {'Window',ppID};
        numCom_Table_W = horzcat(window,windowWeightedMeans_Table(:,4));
        numCom_Table_W.Properties.VariableNames = {'Window',ppID};
        modularity_Table_W = horzcat(window,windowWeightedMeans_Table(:,5));
        modularity_Table_W.Properties.VariableNames = {'Window',ppID};
        participationPos_Table_W = horzcat(window,windowWeightedMeans_Table(:,6));
        participationPos_Table_W.Properties.VariableNames = {'Window',ppID};
        participationNeg_Table_W = horzcat(window,windowWeightedMeans_Table(:,7));
        participationNeg_Table_W.Properties.VariableNames = {'Window',ppID};
        participationPosSD_Table_W = horzcat(window,windowWeightedMeans_Table(:,8));
        participationPosSD_Table_W.Properties.VariableNames = {'Window',ppID};
        participationNegSD_Table_W = horzcat(window,windowWeightedMeans_Table(:,9));
        participationNegSD_Table_W.Properties.VariableNames = {'Window',ppID};
        diversityPos_Table_W = horzcat(window,windowWeightedMeans_Table(:,10));
        diversityPos_Table_W.Properties.VariableNames = {'Window',ppID};
        diversityNeg_Table_W = horzcat(window,windowWeightedMeans_Table(:,11));
        diversityNeg_Table_W.Properties.VariableNames = {'Window',ppID};
        diversityPosSD_Table_W = horzcat(window,windowWeightedMeans_Table(:,12));
        diversityPosSD_Table_W.Properties.VariableNames = {'Window',ppID};
        diversityNegSD_Table_W = horzcat(window,windowWeightedMeans_Table(:,13));
        diversityNegSD_Table_W.Properties.VariableNames = {'Window',ppID};
        comRadius_Table_W = horzcat(window,windowWeightedMeans_Table(:,14));
        comRadius_Table_W.Properties.VariableNames = {'Window',ppID};
        modVA_Table_W = horzcat(window,windowWeightedMeans_Table(:,15));
        modVA_Table_W.Properties.VariableNames = {'Window',ppID};
    else
        clustering_Table_S = horzcat(clustering_Table_S,windowSimpleMeans_Table(:,1));
        clustering_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        clusteringSD_Table_S = horzcat(clusteringSD_Table_S,windowSimpleMeans_Table(:,2));
        clusteringSD_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        density_Table_S = horzcat(density_Table_S,windowSimpleMeans_Table(:,3));
        density_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        numCom_Table_S = horzcat(numCom_Table_S,windowSimpleMeans_Table(:,4));
        numCom_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        modularity_Table_S = horzcat(modularity_Table_S,windowSimpleMeans_Table(:,5));
        modularity_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        participationPos_Table_S = horzcat(participationPos_Table_S,windowSimpleMeans_Table(:,6));
        participationPos_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        participationNeg_Table_S = horzcat(participationNeg_Table_S,windowSimpleMeans_Table(:,7));
        participationNeg_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        participationPosSD_Table_S = horzcat(participationPosSD_Table_S,windowSimpleMeans_Table(:,8));
        participationPosSD_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        participationNegSD_Table_S = horzcat(participationNegSD_Table_S,windowSimpleMeans_Table(:,9));
        participationNegSD_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        diversityPos_Table_S = horzcat(diversityPos_Table_S,windowSimpleMeans_Table(:,10));
        diversityPos_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        diversityNeg_Table_S = horzcat(diversityNeg_Table_S,windowSimpleMeans_Table(:,11));
        diversityNeg_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        diversityPosSD_Table_S = horzcat(diversityPosSD_Table_S,windowSimpleMeans_Table(:,12));
        diversityPosSD_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        diversityNegSD_Table_S = horzcat(diversityNegSD_Table_S,windowSimpleMeans_Table(:,13));
        diversityNegSD_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        comRadius_Table_S = horzcat(comRadius_Table_S,windowSimpleMeans_Table(:,14));
        comRadius_Table_S.Properties.VariableNames{i_subject+1} = ppID;
        modVA_Table_S = horzcat(modVA_Table_S,windowSimpleMeans_Table(:,15));
        modVA_Table_S.Properties.VariableNames{i_subject+1} = ppID;

        clustering_Table_W = horzcat(clustering_Table_W,windowWeightedMeans_Table(:,1));
        clustering_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        clusteringSD_Table_W = horzcat(clusteringSD_Table_W,windowWeightedMeans_Table(:,2));
        clusteringSD_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        density_Table_W = horzcat(density_Table_W,windowWeightedMeans_Table(:,3));
        density_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        numCom_Table_W = horzcat(numCom_Table_W,windowWeightedMeans_Table(:,4));
        numCom_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        modularity_Table_W = horzcat(modularity_Table_W,windowWeightedMeans_Table(:,5));
        modularity_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        participationPos_Table_W = horzcat(participationPos_Table_W,windowWeightedMeans_Table(:,6));
        participationPos_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        participationNeg_Table_W = horzcat(participationNeg_Table_W,windowWeightedMeans_Table(:,7));
        participationNeg_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        participationPosSD_Table_W = horzcat(participationPosSD_Table_W,windowWeightedMeans_Table(:,8));
        participationPosSD_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        participationNegSD_Table_W = horzcat(participationNegSD_Table_W,windowWeightedMeans_Table(:,9));
        participationNegSD_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        diversityPos_Table_W = horzcat(diversityPos_Table_W,windowWeightedMeans_Table(:,10));
        diversityPos_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        diversityNeg_Table_W = horzcat(diversityNeg_Table_W,windowWeightedMeans_Table(:,11));
        diversityNeg_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        diversityPosSD_Table_W = horzcat(diversityPosSD_Table_W,windowWeightedMeans_Table(:,12));
        diversityPosSD_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        diversityNegSD_Table_W = horzcat(diversityNegSD_Table_W,windowWeightedMeans_Table(:,13));
        diversityNegSD_Table_W.Properties.VariableNames{i_subject+1} = ppID;
        comRadius_Table_W = horzcat(comRadius_Table_W,windowWeightedMeans_Table(:,14));
        comRadius_Table_W.Properties.VariableNames{i_subject+1} = ppID;;
        modVA_Table_W = horzcat(modVA_Table_W,windowWeightedMeans_Table(:,15));
        modVA_Table_W.Properties.VariableNames{i_subject+1} = ppID;
    end
    clearvars clusterCoefGlobal clusterCoefGlobalSD kdenW numCommunities modularity participation participationSD diversity diversitySD commRadius modularity_VA diameter numEdges
    clearvars networkMeasures windowSimpleMeans windowWeightedMeans dayIDlist
end

%% create summary tables
save(['dataSet' num2str(dataSet) '_networkMeasureSimpleMeans_windows.mat'],'networkMeasureSimpleMeans');
save(['dataSet' num2str(dataSet) '_networkMeasureSimpleSDs_windows.mat'],'networkMeasureSimpleSDs');
save(['dataSet' num2str(dataSet) '_networkMeasureWeightedMeans_windows.mat'],'networkMeasureWeightedMeans');
save(['dataSet' num2str(dataSet) '_networkMeasureWeightedSDs_windows.mat'],'networkMeasureWeightedSDs');
variableNames = {'cluster','clusterSD','density','numCom','mod','pPos','pNeg','pPosSD','pNegSD','dPos','dNeg','dPosSD','dNegSD','comRadius','modVA'};
networkMeasureSimpleMeans_Table = array2table(networkMeasureSimpleMeans,'VariableNames',variableNames);
networkMeasureSimpleSDs_Table = array2table(networkMeasureSimpleSDs,'VariableNames',variableNames);
networkMeasureWeightedMeans_Table = array2table(networkMeasureWeightedMeans,'VariableNames',variableNames);
networkMeasureWeightedSDs_Table = array2table(networkMeasureWeightedSDs,'VariableNames',variableNames);
if printSummaries == 1
    writetable(networkMeasureSimpleMeans_Table,['dataSet' num2str(dataSet) '_networkMeasureSimpleMeans_windows.xlsx']);
    writetable(networkMeasureSimpleSDs_Table,['dataSet' num2str(dataSet) '_networkMeasureSimpleSDs_windows.xlsx']);
    writetable(networkMeasureWeightedMeans_Table,['dataSet' num2str(dataSet) '_networkMeasureWeightedMeans_windows.xlsx']);
    writetable(networkMeasureWeightedSDs_Table,['dataSet' num2str(dataSet) '_networkMeasureWeightedSDs_windows.xlsx']);
end

%% write tables to file
if printTables == 1
    writetable(clustering_Table_S,['dataSet' num2str(dataSet) '_clustering_windows_simple.xlsx']);
    writetable(clusteringSD_Table_S,['dataSet' num2str(dataSet) '_clustering_SD_windows_simple.xlsx']);
    writetable(density_Table_S,['dataSet' num2str(dataSet) '_density_windows_simple.xlsx']);
    writetable(numCom_Table_S,['dataSet' num2str(dataSet) '_number_communities_windows_simple.xlsx']);
    writetable(modularity_Table_S,['dataSet' num2str(dataSet) '_modularity_windows_simple.xlsx']);
    writetable(participationPos_Table_S,['dataSet' num2str(dataSet) '_participation_pos_windows_simple.xlsx']);
    writetable(participationNeg_Table_S,['dataSet' num2str(dataSet) '_participation_neg_windows_simple.xlsx']);
    writetable(participationPosSD_Table_S,['dataSet' num2str(dataSet) '_participation_SD_pos_windows_simple.xlsx']);
    writetable(participationNegSD_Table_S,['dataSet' num2str(dataSet) '_participation_SD_neg_windows_simple.xlsx']);
    writetable(diversityPos_Table_S,['dataSet' num2str(dataSet) '_diversity_pos_windows_simple.xlsx']);
    writetable(diversityNeg_Table_S,['dataSet' num2str(dataSet) '_diversity_neg_windows_simple.xlsx']);
    writetable(diversityPosSD_Table_S,['dataSet' num2str(dataSet) '_diversity_SD_pos_windows_simple.xlsx']);
    writetable(diversityNegSD_Table_S,['dataSet' num2str(dataSet) '_diversity_SD_neg_windows_simple.xlsx']);
    writetable(comRadius_Table_S,['dataSet' num2str(dataSet) '_community_radius_windows_simple.xlsx']);
    writetable(modVA_Table_S,['dataSet' num2str(dataSet) '_modularity_VA_windows_simple.xlsx']);
   
    writetable(clustering_Table_W,['dataSet' num2str(dataSet) '_clustering_windows_weighted.xlsx']);
    writetable(clusteringSD_Table_W,['dataSet' num2str(dataSet) '_clustering_SD_windows_weighted.xlsx']);
    writetable(density_Table_W,['dataSet' num2str(dataSet) '_density_windows_weighted.xlsx']);
    writetable(numCom_Table_W,['dataSet' num2str(dataSet) '_number_communities_windows_weighted.xlsx']);
    writetable(modularity_Table_W,['dataSet' num2str(dataSet) '_modularity_windows_weighted.xlsx']);
    writetable(participationPos_Table_W,['dataSet' num2str(dataSet) '_participation_pos_windows_weighted.xlsx']);
    writetable(participationNeg_Table_W,['dataSet' num2str(dataSet) '_participation_neg_windows_weighted.xlsx']);
    writetable(participationPosSD_Table_W,['dataSet' num2str(dataSet) '_participation_SD_pos_windows_weighted.xlsx']);
    writetable(participationNegSD_Table_W,['dataSet' num2str(dataSet) '_participation_SD_neg_windows_weighted.xlsx']);
    writetable(diversityPos_Table_W,['dataSet' num2str(dataSet) '_diversity_pos_windows_weighted.xlsx']);
    writetable(diversityNeg_Table_W,['dataSet' num2str(dataSet) '_diversity_neg_windows_weighted.xlsx']);
    writetable(diversityPosSD_Table_W,['dataSet' num2str(dataSet) '_diversity_SD_pos_windows_weighted.xlsx']);
    writetable(diversityNegSD_Table_W,['dataSet' num2str(dataSet) '_diversity_SD_neg_windows_weighted.xlsx']);
    writetable(comRadius_Table_W,['dataSet' num2str(dataSet) '_community_radius_windows_weighted.xlsx']);
    writetable(modVA_Table_W,['dataSet' num2str(dataSet) '_modularity_VA_windows_weighted.xlsx']);
end

%% run Kolmogorov-Smirnov test to assess normality; normalize distributions as needed
[~,numVariables] = size(networkMeasureSimpleMeans);
for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(networkMeasureSimpleMeans(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(networkMeasureSimpleMeans(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm_SimpleMeans(:,i_variable) = norminv(fracRank,mean(networkMeasureSimpleMeans(:,i_variable)),std(networkMeasureSimpleMeans(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm_SimpleMeans(:,i_variable) = networkMeasureSimpleMeans(:,i_variable);
    end
end

for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(networkMeasureSimpleSDs(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(networkMeasureSimpleSDs(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm_SimpleSDs(:,i_variable) = norminv(fracRank,mean(networkMeasureSimpleSDs(:,i_variable)),std(networkMeasureSimpleSDs(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm_SimpleSDs(:,i_variable) = networkMeasureSimpleSDs(:,i_variable);
    end
end

for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(networkMeasureWeightedMeans(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(networkMeasureWeightedMeans(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm_WeightedMeans(:,i_variable) = norminv(fracRank,mean(networkMeasureWeightedMeans(:,i_variable)),std(networkMeasureWeightedMeans(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm_WeightedMeans(:,i_variable) = networkMeasureWeightedMeans(:,i_variable);
    end
end

for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(networkMeasureWeightedSDs(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(networkMeasureWeightedSDs(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm_WeightedSDs(:,i_variable) = norminv(fracRank,mean(networkMeasureWeightedSDs(:,i_variable)),std(networkMeasureWeightedSDs(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm_WeightedSDs(:,i_variable) = networkMeasureWeightedSDs(:,i_variable);
    end
end

%% compile all (normalized) data and re-export
save(['dataSet' num2str(dataSet) '_networkMeasureSimpleMeans_normalized_windows.mat'],'fracRankNorm_SimpleMeans');
save(['dataSet' num2str(dataSet) '_networkMeasureSimpleSDs_normalized_windows.mat'],'fracRankNorm_SimpleSDs');
save(['dataSet' num2str(dataSet) '_networkMeasureWeightedMeans_normalized_windows.mat'],'fracRankNorm_WeightedMeans');
save(['dataSet' num2str(dataSet) '_networkMeasureWeightedSDs_normalized_windows.mat'],'fracRankNorm_WeightedSDs');
fracRankNorm_SimpleMeans_Table = array2table(fracRankNorm_SimpleMeans,'VariableNames',variableNames);
fracRankNorm_SimpleSDs_Table = array2table(fracRankNorm_SimpleSDs,'VariableNames',variableNames);
fracRankNorm_WeightedMeans_Table = array2table(fracRankNorm_WeightedMeans,'VariableNames',variableNames);
fracRankNorm_WeightedSDs_Table = array2table(fracRankNorm_WeightedSDs,'VariableNames',variableNames);
if printSummaries == 1
    writetable(fracRankNorm_SimpleMeans_Table,['dataSet' num2str(dataSet) '_networkMeasureSimpleMeans_normalized_windows.xlsx']);
    writetable(fracRankNorm_SimpleSDs_Table,['dataSet' num2str(dataSet) '_networkMeasureSimpleSDs_normalized_windows.xlsx']);
    writetable(fracRankNorm_WeightedMeans_Table,['dataSet' num2str(dataSet) '_networkMeasureWeightedMeans_normalized_windows.xlsx']);
    writetable(fracRankNorm_WeightedSDs_Table,['dataSet' num2str(dataSet) '_networkMeasureWeightedSDs_normalized_windows.xlsx']);
end

%% save tables as matrices
if saveData == 1
    subjectCol = [0; subjectIDlist];

    clustering_Array_S = table2array(clustering_Table_S)';
    clustering_Array_S = horzcat(subjectCol,clustering_Array_S);
    save(['dataSet' num2str(dataSet) '_clustering_windows_simple.mat'],'clustering_Array_S');

    clusteringSD_Array_S = table2array(clusteringSD_Table_S)';
    clusteringSD_Array_S = horzcat(subjectCol,clusteringSD_Array_S);
    save(['dataSet' num2str(dataSet) '_clusteringSD_windows_simple.mat'],'clusteringSD_Array_S');

    density_Array_S = table2array(density_Table_S)';
    density_Array_S = horzcat(subjectCol,density_Array_S);
    save(['dataSet' num2str(dataSet) '_density_windows_simple.mat'],'density_Array_S');

    numCom_Array_S = table2array(numCom_Table_S)';
    numCom_Array_S = horzcat(subjectCol,numCom_Array_S);
    save(['dataSet' num2str(dataSet) '_number_communities_windows_simple.mat'],'numCom_Array_S');

    modularity_Array_S = table2array(modularity_Table_S)';
    modularity_Array_S = horzcat(subjectCol,modularity_Array_S);
    save(['dataSet' num2str(dataSet) '_modularity_windows_simple.mat'],'modularity_Array_S');

    participationPos_Array_S = table2array(participationPos_Table_S)';
    participationPos_Array_S = horzcat(subjectCol,participationPos_Array_S);
    save(['dataSet' num2str(dataSet) '_participation_pos_windows_simple.mat'],'participationPos_Array_S');

    participationNeg_Array_S = table2array(participationNeg_Table_S)';
    participationNeg_Array_S = horzcat(subjectCol,participationNeg_Array_S);
    save(['dataSet' num2str(dataSet) '_participation_neg_windows_simple.mat'],'participationNeg_Array_S');

    participationPosSD_Array_S = table2array(participationPosSD_Table_S)';
    participationPosSD_Array_S = horzcat(subjectCol,participationPosSD_Array_S);
    save(['dataSet' num2str(dataSet) '_participation_pos_SD_windows_simple.mat'],'participationPosSD_Array_S');

    participationNegSD_Array_S = table2array(participationNegSD_Table_S)';
    participationNegSD_Array_S = horzcat(subjectCol,participationNegSD_Array_S);
    save(['dataSet' num2str(dataSet) '_participation_neg_SD_windows_simple.mat'],'participationNegSD_Array_S');

    diversityPos_Array_S = table2array(diversityPos_Table_S)';
    diversityPos_Array_S = horzcat(subjectCol,diversityPos_Array_S);
    save(['dataSet' num2str(dataSet) '_diversity_pos_windows_simple.mat'],'diversityPos_Array_S');

    diversityNeg_Array_S = table2array(diversityNeg_Table_S)';
    diversityNeg_Array_S = horzcat(subjectCol,diversityNeg_Array_S);
    save(['dataSet' num2str(dataSet) '_diversity_neg_windows_simple.mat'],'diversityNeg_Array_S');

    diversityPosSD_Array_S = table2array(diversityPosSD_Table_S)';
    diversityPosSD_Array_S = horzcat(subjectCol,diversityPosSD_Array_S);
    save(['dataSet' num2str(dataSet) '_diversity_pos_SD_windows_simple.mat'],'diversityPosSD_Array_S');

    diversityNegSD_Array_S = table2array(diversityNegSD_Table_S)';
    diversityNegSD_Array_S = horzcat(subjectCol,diversityNegSD_Array_S);
    save(['dataSet' num2str(dataSet) '_diversity_neg_SD_windows_simple.mat'],'diversityNegSD_Array_S');

    comRadius_Array_S = table2array(comRadius_Table_S)';
    comRadius_Array_S = horzcat(subjectCol,comRadius_Array_S);
    save(['dataSet' num2str(dataSet) '_community_radius_windows_simple.mat'],'comRadius_Array_S');

    modVA_Array_S = table2array(modVA_Table_S)';
    modVA_Array_S = horzcat(subjectCol,modVA_Array_S);
    save(['dataSet' num2str(dataSet) '_modularity_VA_windows_simple.mat'],'modVA_Array_S');
    
    clustering_Array_W = table2array(clustering_Table_W)';
    clustering_Array_W = horzcat(subjectCol,clustering_Array_W);
    save(['dataSet' num2str(dataSet) '_clustering_windows_weighted.mat'],'clustering_Array_W');

    clusteringSD_Array_W = table2array(clusteringSD_Table_W)';
    clusteringSD_Array_W = horzcat(subjectCol,clusteringSD_Array_W);
    save(['dataSet' num2str(dataSet) '_clusteringSD_windows_weighted.mat'],'clusteringSD_Array_W');

    density_Array_W = table2array(density_Table_W)';
    density_Array_W = horzcat(subjectCol,density_Array_W);
    save(['dataSet' num2str(dataSet) '_density_windows_weighted.mat'],'density_Array_W');

    numCom_Array_W = table2array(numCom_Table_W)';
    numCom_Array_W = horzcat(subjectCol,numCom_Array_W);
    save(['dataSet' num2str(dataSet) '_number_communities_windows_weighted.mat'],'numCom_Array_W');

    modularity_Array_W = table2array(modularity_Table_W)';
    modularity_Array_W = horzcat(subjectCol,modularity_Array_W);
    save(['dataSet' num2str(dataSet) '_modularity_windows_weighted.mat'],'modularity_Array_W');

    participationPos_Array_W = table2array(participationPos_Table_W)';
    participationPos_Array_W = horzcat(subjectCol,participationPos_Array_W);
    save(['dataSet' num2str(dataSet) '_participation_pos_windows_weighted.mat'],'participationPos_Array_W');

    participationNeg_Array_W = table2array(participationNeg_Table_W)';
    participationNeg_Array_W = horzcat(subjectCol,participationNeg_Array_W);
    save(['dataSet' num2str(dataSet) '_participation_neg_windows_weighted.mat'],'participationNeg_Array_W');

    participationPosSD_Array_W = table2array(participationPosSD_Table_W)';
    participationPosSD_Array_W = horzcat(subjectCol,participationPosSD_Array_W);
    save(['dataSet' num2str(dataSet) '_participation_pos_SD_windows_weighted.mat'],'participationPosSD_Array_W');

    participationNegSD_Array_W = table2array(participationNegSD_Table_W)';
    participationNegSD_Array_W = horzcat(subjectCol,participationNegSD_Array_W);
    save(['dataSet' num2str(dataSet) '_participation_neg_SD_windows_weighted.mat'],'participationNegSD_Array_W');

    diversityPos_Array_W = table2array(diversityPos_Table_W)';
    diversityPos_Array_W = horzcat(subjectCol,diversityPos_Array_W);
    save(['dataSet' num2str(dataSet) '_diversity_pos_windows_weighted.mat'],'diversityPos_Array_W');

    diversityNeg_Array_W = table2array(diversityNeg_Table_W)';
    diversityNeg_Array_W = horzcat(subjectCol,diversityNeg_Array_W);
    save(['dataSet' num2str(dataSet) '_diversity_neg_windows_weighted.mat'],'diversityNeg_Array_W');

    diversityPosSD_Array_W = table2array(diversityPosSD_Table_W)';
    diversityPosSD_Array_W = horzcat(subjectCol,diversityPosSD_Array_W);
    save(['dataSet' num2str(dataSet) '_diversity_pos_SD_windows_weighted.mat'],'diversityPosSD_Array_W');

    diversityNegSD_Array_W = table2array(diversityNegSD_Table_W)';
    diversityNegSD_Array_W = horzcat(subjectCol,diversityNegSD_Array_W);
    save(['dataSet' num2str(dataSet) '_diversity_neg_SD_windows_weighted.mat'],'diversityNegSD_Array_W');

    comRadius_Array_W = table2array(comRadius_Table_W)';
    comRadius_Array_W = horzcat(subjectCol,comRadius_Array_W);
    save(['dataSet' num2str(dataSet) '_community_radius_windows_weighted.mat'],'comRadius_Array_W');

    modVA_Array_W = table2array(modVA_Table_W)';
    modVA_Array_W = horzcat(subjectCol,modVA_Array_W);
    save(['dataSet' num2str(dataSet) '_modularity_VA_windows_weighted.mat'],'modVA_Array_W'); 
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
        [r_GADm_S,p_GADm_S] = corr(networkMeasureSimpleMeans,GADtotal,'rows','complete');
        [r_PHQ8m_S,p_PHQ8m_S] = corr(networkMeasureSimpleMeans,PHQ8total,'rows','complete');
        [r_GADsd_S,p_GADsd_S] = corr(networkMeasureSimpleSDs,GADtotal,'rows','complete');
        [r_PHQ8sd_S,p_PHQ8sd_S] = corr(networkMeasureSimpleSDs,PHQ8total,'rows','complete');
        
        [r_GADm_W,p_GADm_W] = corr(networkMeasureWeightedMeans,GADtotal,'rows','complete');
        [r_PHQ8m_W,p_PHQ8m_W] = corr(networkMeasureWeightedMeans,PHQ8total,'rows','complete');
        [r_GADsd_W,p_GADsd_W] = corr(networkMeasureWeightedSDs,GADtotal,'rows','complete');
        [r_PHQ8sd_W,p_PHQ8sd_W] = corr(networkMeasureWeightedSDs,PHQ8total,'rows','complete');
    elseif dataSet == 39
        questionnaires = importdata('39ESdata_QuestionnaireData.xlsx');
        BAItotal = questionnaires.data(:,9); 
        BDItotal = questionnaires.data(:,12);
        [r_BAIm_S,p_BAIm_S] = corr(networkMeasureSimpleMeans,BAItotal,'rows','complete');
        [r_BDIm_S,p_BDIm_S] = corr(networkMeasureSimpleMeans,BDItotal,'rows','complete');
        [r_BAIsd_S,p_BAIsd_S] = corr(networkMeasureSimpleSDs,BAItotal,'rows','complete');
        [r_BDIsd_S,p_BDIsd_S] = corr(networkMeasureSimpleSDs,BDItotal,'rows','complete');
        
        [r_BAIm_W,p_BAIm_W] = corr(networkMeasureWeightedMeans,BAItotal,'rows','complete');
        [r_BDIm_W,p_BDIm_W] = corr(networkMeasureWeightedMeans,BDItotal,'rows','complete');
        [r_BAIsd_W,p_BAIsd_W] = corr(networkMeasureWeightedSDs,BAItotal,'rows','complete');
        [r_BDIsd_W,p_BDIsd_W] = corr(networkMeasureWeightedSDs,BDItotal,'rows','complete');
    end
    [r_zICCm_S,p_zICCm_S] = corr(networkMeasureSimpleMeans,zICC,'rows','complete');
    [r_negm_S,p_negm_S] = corr(networkMeasureSimpleMeans,neg,'rows','complete');
    [r_negSDm_S,p_negSDm_S] = corr(networkMeasureSimpleMeans,negSD,'rows','complete');
    [r_zICCsd_S,p_zICCsd_S] = corr(networkMeasureSimpleSDs,zICC,'rows','complete');
    [r_negsd_S,p_negsd_S] = corr(networkMeasureSimpleSDs,neg,'rows','complete');
    [r_negSDsd_S,p_negSDsd_S] = corr(networkMeasureSimpleSDs,negSD,'rows','complete');
    
    [r_zICCm_W,p_zICCm_W] = corr(networkMeasureWeightedMeans,zICC,'rows','complete');
    [r_negm_W,p_negm_W] = corr(networkMeasureWeightedMeans,neg,'rows','complete');
    [r_negSDm_W,p_negSDm_W] = corr(networkMeasureWeightedMeans,negSD,'rows','complete');
    [r_zICCsd_W,p_zICCsd_W] = corr(networkMeasureWeightedSDs,zICC,'rows','complete');
    [r_negsd_W,p_negsd_W] = corr(networkMeasureWeightedSDs,neg,'rows','complete');
    [r_negSDsd_W,p_negSDsd_W] = corr(networkMeasureWeightedSDs,negSD,'rows','complete');
end