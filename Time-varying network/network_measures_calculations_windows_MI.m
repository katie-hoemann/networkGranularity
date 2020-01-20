% NOTES: Currently built out for 3-day sliding window; modification to
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
        % calculate mutual information matrix
        rows = 1:1:size(windowData,2);
        for i_var = 1:length(rows)
            for j_var = rows(rows~=i_var)
                windowMatrix(i_var,j_var) = nmi(windowData(:,i_var),windowData(:,j_var));
            end
        end
        % delete zero-variance rows/columns from data
        colVar = var(windowMatrix,0,1); % find columns (nodes) with zero variance
        rowVar = var(windowMatrix,0,2); % find rows (instances) with zero variance
        missingColVar = find(colVar==0); % index zero-variance nodes
        missingRowVar = find(rowVar==0); % index zero-variance instances
        windowMatrix(:,missingColVar) = []; % remove zero-variance nodes
        windowMatrix(missingRowVar,:) = []; % remove zero-variance instances
        if numel(windowMatrix) == 0
            continue
        end
        % calculate clustering coefficient (original connection weights)
        [C,C_null,Ctot,Ctot_null] = clustering_coef_wu_sign(windowMatrix,3);
        clusterCoefGlobal(i_window) = Ctot; % save off global value
        clusterCoefGlobalSD(i_window) = std(C); % save off SD of node-level measure
        % calculate weighted network density (original connection weigthts)
        kden(i_window) = density_und(windowMatrix);
        halfWindowMatrix = triu(windowMatrix);
        halfWindowMatrix = reshape(halfWindowMatrix,[numel(halfWindowMatrix),1]);
        halfWindowMatrix = halfWindowMatrix(halfWindowMatrix~=0);
        kdenW(i_window) = mean(halfWindowMatrix); 
        % run community detection (original connection weights)
        gamma = 1;
        W = windowMatrix;  % set subject matrix
        B = 'modularity';
        n = size(W,1);      % number of nodes
        M = 1:n;            % initial community affiliations
        Q0 = -1; Q = 0;    % initialize modularity values
            while Q - Q0 > 1e-5    % while modularity increases
                Q0 = Q;            % perform community detection
                [M, Q] = community_louvain(W,gamma,[],B); 
            end
        communityAssignment = M; % save community assignment vector    
        numCommunities(i_window) = max(M);   % save number of communities
        modularity(i_window) = Q; % save modularity value
        % calculate participation & diversity coefficients (original connection weights)
        [Ppos Pneg] = participation_coef_sign(W,M); % participation coefficient
        participation(i_window,1) = mean(Ppos,1);
        participation(i_window,2) = mean(Pneg,1);
        participationSD(i_window,1) = std(Ppos,1);
        participationSD(i_window,2) = std(Pneg,1);
        [Hpos Hneg] = diversity_coef_sign(W,M); % diversity coefficient
        diversity(i_window,1) = mean(Hpos,1);
        diversity(i_window,2) = mean(Hneg,1);
        diversitySD(i_window,1) = std(Hpos,1);
        diversitySD(i_window,2) = std(Hneg,1);
        % calculate community radius (original connection weights)
        locations = table2array(words(:,[3 4]));
        locations(missingColVar,:) = [];
        partitions = communityAssignment;
        commRadiusArray = comm_radius(partitions,locations);
        commRadius(:,i_window) = commRadiusArray(end,end);
        % run community assignment (original connection weights)
        windowM0 = assignedM; % copy overall community affiliation vector
        windowM0(missingColVar) = []; % delete missing nodes from community affiliation vector
        M0 = windowM0; % set initial community affiliations vector
        [~, Q_VA] = community_louvain(W,gamma,M0,B);
        modularity_VA(i_window) = Q_VA; % save modularity value
        fprintf('participant %d, window %d\n',subjectID,i_window)
        clearvars windowMatrix communityAssignment
    end
    %% calculate subject-level mean and SD for measures
    networkMeasureMeans(i_subject,1) = mean(clusterCoefGlobal,'omitnan');
    networkMeasureMeans(i_subject,2) = mean(clusterCoefGlobalSD,'omitnan');
    networkMeasureMeans(i_subject,3) = mean(kdenW,'omitnan');
    networkMeasureMeans(i_subject,4) = mean(numCommunities,'omitnan');
    networkMeasureMeans(i_subject,5) = mean(modularity,'omitnan');
    networkMeasureMeans(i_subject,6) = mean(participation(:,1),'omitnan');
    networkMeasureMeans(i_subject,7) = mean(participation(:,2),'omitnan');
    networkMeasureMeans(i_subject,8) = mean(participationSD(:,1),'omitnan');
    networkMeasureMeans(i_subject,9) = mean(participationSD(:,2),'omitnan');
    networkMeasureMeans(i_subject,10) = mean(diversity(:,1),'omitnan');
    networkMeasureMeans(i_subject,11) = mean(diversity(:,2),'omitnan');
    networkMeasureMeans(i_subject,12) = mean(diversitySD(:,1),'omitnan');
    networkMeasureMeans(i_subject,13) = mean(diversitySD(:,2),'omitnan');
    networkMeasureMeans(i_subject,14) = mean(commRadius,'omitnan');
    networkMeasureMeans(i_subject,15) = mean(modularity_VA,'omitnan');
    
    networkMeasureSDs(i_subject,1) = std(clusterCoefGlobal,'omitnan');
    networkMeasureSDs(i_subject,2) = std(clusterCoefGlobalSD,'omitnan');
    networkMeasureSDs(i_subject,3) = std(kdenW,'omitnan');
    networkMeasureSDs(i_subject,4) = std(numCommunities,'omitnan');
    networkMeasureSDs(i_subject,5) = std(modularity,'omitnan');
    networkMeasureSDs(i_subject,6) = std(participation(:,1),'omitnan');
    networkMeasureSDs(i_subject,7) = std(participation(:,2),'omitnan');
    networkMeasureSDs(i_subject,8) = std(participationSD(:,1),'omitnan');
    networkMeasureSDs(i_subject,9) = std(participationSD(:,2),'omitnan');
    networkMeasureSDs(i_subject,10) = std(diversity(:,1),'omitnan');
    networkMeasureSDs(i_subject,11) = std(diversity(:,2),'omitnan');
    networkMeasureSDs(i_subject,12) = std(diversitySD(:,1),'omitnan');   
    networkMeasureSDs(i_subject,13) = std(diversitySD(:,2),'omitnan'); 
    networkMeasureSDs(i_subject,14) = std(commRadius,'omitnan'); 
    networkMeasureSDs(i_subject,15) = std(modularity_VA,'omitnan'); 
    %% add subject results to summary table
    ppID = ['PP' num2str(subjectID)];
    windowNumbers = dayIDlist(1:end-2);
    if i_subject == 1
        window = (1:1:max(allData(:,2))-2)';
        window = array2table(window,'VariableNames',{'Window'});
     
        clustering_Table = tablecompile_start(window,windowNumbers,clusterCoefGlobal',ppID,'Window','ppWindow','ppWindow');
        clusteringSD_Table = tablecompile_start(window,windowNumbers,clusterCoefGlobalSD',ppID,'Window','ppWindow','ppWindow');
        density_Table = tablecompile_start(window,windowNumbers,kdenW',ppID,'Window','ppWindow','ppWindow');
        numCom_Table = tablecompile_start(window,windowNumbers,numCommunities',ppID,'Window','ppWindow','ppWindow');
        modularity_Table = tablecompile_start(window,windowNumbers,modularity',ppID,'Window','ppWindow','ppWindow');
        participationPos_Table = tablecompile_start(window,windowNumbers,participation(:,1),ppID,'Window','ppWindow','ppWindow');
        participationNeg_Table = tablecompile_start(window,windowNumbers,participation(:,2),ppID,'Window','ppWindow','ppWindow');
        participationPosSD_Table = tablecompile_start(window,windowNumbers,participationSD(:,1),ppID,'Window','ppWindow','ppWindow');
        participationNegSD_Table = tablecompile_start(window,windowNumbers,participationSD(:,2),ppID,'Window','ppWindow','ppWindow');
        diversityPos_Table = tablecompile_start(window,windowNumbers,diversity(:,1),ppID,'Window','ppWindow','ppWindow');
        diversityNeg_Table = tablecompile_start(window,windowNumbers,diversity(:,2),ppID,'Window','ppWindow','ppWindow');
        diversityPosSD_Table = tablecompile_start(window,windowNumbers,diversitySD(:,1),ppID,'Window','ppWindow','ppWindow');
        diversityNegSD_Table = tablecompile_start(window,windowNumbers,diversitySD(:,2),ppID,'Window','ppWindow','ppWindow');
        comRadius_Table = tablecompile_start(window,windowNumbers,commRadius',ppID,'Window','ppWindow','ppWindow');
        modVA_Table = tablecompile_start(window,windowNumbers,modularity_VA',ppID,'Window','ppWindow','ppWindow');
    else
        clustering_Table = tablecompile_iter(clustering_Table,windowNumbers,clusterCoefGlobal',ppID,'Window','ppWindow','ppWindow');
        clusteringSD_Table = tablecompile_iter(clusteringSD_Table,windowNumbers,clusterCoefGlobalSD',ppID,'Window','ppWindow','ppWindow');
        density_Table = tablecompile_iter(density_Table,windowNumbers,kdenW',ppID,'Window','ppWindow','ppWindow');
        numCom_Table = tablecompile_iter(numCom_Table,windowNumbers,numCommunities',ppID,'Window','ppWindow','ppWindow');
        modularity_Table = tablecompile_iter(modularity_Table,windowNumbers,modularity',ppID,'Window','ppWindow','ppWindow');
        participationPos_Table = tablecompile_iter(participationPos_Table,windowNumbers,participation(:,1),ppID,'Window','ppWindow','ppWindow');
        participationNeg_Table = tablecompile_iter(participationNeg_Table,windowNumbers,participation(:,2),ppID,'Window','ppWindow','ppWindow');
        participationPosSD_Table = tablecompile_iter(participationPosSD_Table,windowNumbers,participationSD(:,1),ppID,'Window','ppWindow','ppWindow');
        participationNegSD_Table = tablecompile_iter(participationNegSD_Table,windowNumbers,participationSD(:,2),ppID,'Window','ppWindow','ppWindow');
        diversityPos_Table = tablecompile_iter(diversityPos_Table,windowNumbers,diversity(:,1),ppID,'Window','ppWindow','ppWindow');
        diversityNeg_Table = tablecompile_iter(diversityNeg_Table,windowNumbers,diversity(:,2),ppID,'Window','ppWindow','ppWindow');
        diversityPosSD_Table = tablecompile_iter(diversityPosSD_Table,windowNumbers,diversitySD(:,1),ppID,'Window','ppWindow','ppWindow');
        diversityNegSD_Table = tablecompile_iter(diversityNegSD_Table,windowNumbers,diversitySD(:,2),ppID,'Window','ppWindow','ppWindow');
        comRadius_Table = tablecompile_iter(comRadius_Table,windowNumbers,commRadius',ppID,'Window','ppWindow','ppWindow');
        modVA_Table = tablecompile_iter(modVA_Table,windowNumbers,modularity_VA',ppID,'Window','ppWindow','ppWindow');
    end
    clear clusterCoefGlobal clusterCoefGlobalSD kdenW numCommunities modularity participation participationSD diversity diversitySD commRadius modularity_VA diameter numEdges
end

%% create summary tables and write to file
save(['dataSet' num2str(dataSet) '_networkMeasureMeans_windows.mat'],'networkMeasureMeans');
save(['dataSet' num2str(dataSet) '_networkMeasureSDs_windows.mat'],'networkMeasureSDs');
variableNames = {'cluster','clusterSD','density','numCom','mod','pPos','pNeg','pPosSD','pNegSD','dPos','dNeg','dPosSD','dNegSD','comRadius','modVA'};
networkMeasureMeans_Table = array2table(networkMeasureMeans,'VariableNames',variableNames);
networkMeasureSDs_Table = array2table(networkMeasureSDs,'VariableNames',variableNames);
if printSummaries == 1
    writetable(networkMeasureMeans_Table,['dataSet' num2str(dataSet) '_networkMeasureMeans_windows.xlsx']);
    writetable(networkMeasureSDs_Table,['dataSet' num2str(dataSet) '_networkMeasureSDs_windows.xlsx']);
end

%% write measure-specific tables to file
if printTables == 1    
    writetable(clustering_Table,['dataSet' num2str(dataSet) '_clustering_windows.xlsx']);
    writetable(clusteringSD_Table,['dataSet' num2str(dataSet) '_clustering_SD_windows.xlsx']);
    writetable(density_Table,['dataSet' num2str(dataSet) '_density_windows.xlsx']);
    writetable(numCom_Table,['dataSet' num2str(dataSet) '_number_communities_windows.xlsx']);
    writetable(modularity_Table,['dataSet' num2str(dataSet) '_modularity_windows.xlsx']);
    writetable(participationPos_Table,['dataSet' num2str(dataSet) '_participation_pos_windows.xlsx']);
    writetable(participationNeg_Table,['dataSet' num2str(dataSet) '_participation_neg_windows.xlsx']);
    writetable(participationPosSD_Table,['dataSet' num2str(dataSet) '_participation_SD_pos_windows.xlsx']);
    writetable(participationNegSD_Table,['dataSet' num2str(dataSet) '_participation_SD_neg_windows.xlsx']);
    writetable(diversityPos_Table,['dataSet' num2str(dataSet) '_diversity_pos_windows.xlsx']);
    writetable(diversityNeg_Table,['dataSet' num2str(dataSet) '_diversity_neg_windows.xlsx']);
    writetable(diversityPosSD_Table,['dataSet' num2str(dataSet) '_diversity_SD_pos_windows.xlsx']);
    writetable(diversityNegSD_Table,['dataSet' num2str(dataSet) '_diversity_SD_neg_windows.xlsx']);
    writetable(comRadius_Table,['dataSet' num2str(dataSet) '_community_radius_windows.xlsx']);
    writetable(modVA_Table,['dataSet' num2str(dataSet) '_modularity_VA_windows.xlsx']);
end

%% run Kolmogorov-Smirnov test to assess normality; normalize distributions as needed
[~,numVariables] = size(networkMeasureMeans);
for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(networkMeasureMeans(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(networkMeasureMeans(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm_Means(:,i_variable) = norminv(fracRank,mean(networkMeasureMeans(:,i_variable)),std(networkMeasureMeans(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm_Means(:,i_variable) = networkMeasureMeans(:,i_variable);
    end
end

for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(networkMeasureSDs(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(networkMeasureSDs(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm_SDs(:,i_variable) = norminv(fracRank,mean(networkMeasureSDs(:,i_variable)),std(networkMeasureSDs(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm_SDs(:,i_variable) = networkMeasureSDs(:,i_variable);
    end
end

%% compile all (normalized) data and re-export
save(['dataSet' num2str(dataSet) '_networkMeasureMeans_normalized_windows.mat'],'fracRankNorm_Means');
save(['dataSet' num2str(dataSet) '_networkMeasureSDs_normalized_windows.mat'],'fracRankNorm_SDs');
fracRankNorm_Means_Table = array2table(fracRankNorm_Means,'VariableNames',variableNames);
fracRankNorm_SDs_Table = array2table(fracRankNorm_SDs,'VariableNames',variableNames);
if printSummaries == 1
    writetable(fracRankNorm_Means_Table,['dataSet' num2str(dataSet) '_networkMeasureMeans_normalized_windows.xlsx']);
    writetable(fracRankNorm_SDs_Table,['dataSet' num2str(dataSet) '_networkMeasureSDs_normalized_windows.xlsx']);
end

%% save tables as matrices
if saveData == 1
    subjectCol = [0; subjectIDlist];

    clustering_Array = table2array(clustering_Table)';
    clustering_Array = horzcat(subjectCol,clustering_Array);
    save(['dataSet' num2str(dataSet) '_clustering_windows.mat'],'clustering_Array');

    clusteringSD_Array = table2array(clusteringSD_Table)';
    clusteringSD_Array = horzcat(subjectCol,clusteringSD_Array);
    save(['dataSet' num2str(dataSet) '_clusteringSD_windows.mat'],'clusteringSD_Array');

    density_Array = table2array(density_Table)';
    density_Array = horzcat(subjectCol,density_Array);
    save(['dataSet' num2str(dataSet) '_density_windows.mat'],'density_Array');

    numCom_Array = table2array(numCom_Table)';
    numCom_Array = horzcat(subjectCol,numCom_Array);
    save(['dataSet' num2str(dataSet) '_number_communities_windows.mat'],'numCom_Array');

    modularity_Array = table2array(modularity_Table)';
    modularity_Array = horzcat(subjectCol,modularity_Array);
    save(['dataSet' num2str(dataSet) '_modularity_windows.mat'],'modularity_Array');

    participationPos_Array = table2array(participationPos_Table)';
    participationPos_Array = horzcat(subjectCol,participationPos_Array);
    save(['dataSet' num2str(dataSet) '_participation_pos_windows.mat'],'participationPos_Array');

    participationNeg_Array = table2array(participationNeg_Table)';
    participationNeg_Array = horzcat(subjectCol,participationNeg_Array);
    save(['dataSet' num2str(dataSet) '_participation_neg_windows.mat'],'participationNeg_Array');

    participationPosSD_Array = table2array(participationPosSD_Table)';
    participationPosSD_Array = horzcat(subjectCol,participationPosSD_Array);
    save(['dataSet' num2str(dataSet) '_participation_pos_SD_windows.mat'],'participationPosSD_Array');

    participationNegSD_Array = table2array(participationNegSD_Table)';
    participationNegSD_Array = horzcat(subjectCol,participationNegSD_Array);
    save(['dataSet' num2str(dataSet) '_participation_neg_SD_windows.mat'],'participationNegSD_Array');

    diversityPos_Array = table2array(diversityPos_Table)';
    diversityPos_Array = horzcat(subjectCol,diversityPos_Array);
    save(['dataSet' num2str(dataSet) '_diversity_pos_windows.mat'],'diversityPos_Array');

    diversityNeg_Array = table2array(diversityNeg_Table)';
    diversityNeg_Array = horzcat(subjectCol,diversityNeg_Array);
    save(['dataSet' num2str(dataSet) '_diversity_neg_windows.mat'],'diversityNeg_Array');

    diversityPosSD_Array = table2array(diversityPosSD_Table)';
    diversityPosSD_Array = horzcat(subjectCol,diversityPosSD_Array);
    save(['dataSet' num2str(dataSet) '_diversity_pos_SD_windows.mat'],'diversityPosSD_Array');

    diversityNegSD_Array = table2array(diversityNegSD_Table)';
    diversityNegSD_Array = horzcat(subjectCol,diversityNegSD_Array);
    save(['dataSet' num2str(dataSet) '_diversity_neg_SD_windows.mat'],'diversityNegSD_Array');

    comRadius_Array = table2array(comRadius_Table)';
    comRadius_Array = horzcat(subjectCol,comRadius_Array);
    save(['dataSet' num2str(dataSet) '_community_radius_windows.mat'],'comRadius_Array');

    modVA_Array = table2array(modVA_Table)';
    modVA_Array = horzcat(subjectCol,modVA_Array);
    save(['dataSet' num2str(dataSet) '_modularity_VA_windows.mat'],'modVA_Array');
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
        [r_GADm,p_GADm] = corr(networkMeasureMeans,GADtotal,'rows','complete');
        [r_PHQ8m,p_PHQ8m] = corr(networkMeasureMeans,PHQ8total,'rows','complete');
        [r_GADsd,p_GADsd] = corr(networkMeasureSDs,GADtotal,'rows','complete');
        [r_PHQ8sd,p_PHQ8sd] = corr(networkMeasureSDs,PHQ8total,'rows','complete');
    elseif dataSet == 39
        questionnaires = importdata('39ESdata_QuestionnaireData.xlsx');
        BAItotal = questionnaires.data(:,9); 
        BDItotal = questionnaires.data(:,12);
        [r_BAIm,p_BAIm] = corr(networkMeasureMeans,BAItotal,'rows','complete');
        [r_BDIm,p_BDIm] = corr(networkMeasureMeans,BDItotal,'rows','complete');
        [r_BAIsd,p_BAIsd] = corr(networkMeasureSDs,BAItotal,'rows','complete');
        [r_BDIsd,p_BDIsd] = corr(networkMeasureSDs,BDItotal,'rows','complete');
    end
    [r_zICCm,p_zICCm] = corr(networkMeasureMeans,zICC,'rows','complete');
    [r_negm,p_negm] = corr(networkMeasureMeans,neg,'rows','complete');
    [r_negSDm,p_negSDm] = corr(networkMeasureMeans,negSD,'rows','complete');
    [r_zICCsd,p_zICCsd] = corr(networkMeasureSDs,zICC,'rows','complete');
    [r_negsd,p_negsd] = corr(networkMeasureSDs,neg,'rows','complete');
    [r_negSDsd,p_negSDsd] = corr(networkMeasureSDs,negSD,'rows','complete');
end