% This script generates, for every subject in the specified data set, a
% spreadsheet for 3 participation coefficients per node, each at 3 
% pre-specified connection weight thresholds
% Across all subjects in the data set, this script generates 3 spreadsheets
% of the participation coefficients at each threshold,
% as well as histograms of the distribution for each coefficient at each threshold
% This script also generates a spreadsheet of the standard deviations of
% the coefficients at each threshold (no histograms)

clear;
clc;

%% specify dataset
dataSet = 88; % set to 18, 39, or 88
polychoric = 1; % set to 0 to use standard Pearson correlation matrix
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

    %% participation coefficients on original connection weights
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
    B = 'negative_asym';
    n = size(W,1);      % number of nodes
    M = 1:n;            % initial community affiliations
    Q0 = -1; Q = 0;    % initialize modularity values
        while Q - Q0 > 1e-5    % while modularity increases
            Q0 = Q;            % perform community detection
            [M, Q] = community_louvain(W,gamma,[],B); 
        end
    [Ppos(:,1) Pneg(:,1)] = participation_coef_sign(W,M);
    [Gpos1(:,1) Gneg1(:,1)] = gateway_coef_sign(W,M,1); % gateway coefficient with node strength centrality
    [Gpos2(:,1) Gneg2(:,1)] = gateway_coef_sign(W,M,2); % gateway coefficient with betweenness centrality
    [Hpos(:,1) Hneg(:,1)] = diversity_coef_sign(W,M); % diversity coefficient
    
    %% participation coefficients on thresholded connection weights, first test value
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
    B = 'negative_asym';
    n = size(W,1);      % number of nodes
    M = 1:n;            % initial community affiliations
    Q0 = -1; Q = 0;    % initialize modularity values
        while Q - Q0 > 1e-5    % while modularity increases
            Q0 = Q;            % perform community detection
            [M, Q] = community_louvain(W,gamma,[],B); 
        end
    [Ppos(:,2) Pneg(:,2)] = participation_coef_sign(W,M);
    [Gpos1(:,2) Gneg1(:,2)] = gateway_coef_sign(W,M,1); % gateway coefficient with node strength centrality
    [Gpos2(:,2) Gneg2(:,2)] = gateway_coef_sign(W,M,2); % gateway coefficient with betweenness centrality
    [Hpos(:,2) Hneg(:,2)] = diversity_coef_sign(W,M); % diversity coefficient

    %% participation coefficients on thresholded connection weights, second test value
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
    B = 'negative_asym';
    n = size(W,1);      % number of nodes
    M = 1:n;            % initial community affiliations
    Q0 = -1; Q = 0;    % initialize modularity values
        while Q - Q0 > 1e-5    % while modularity increases
            Q0 = Q;            % perform community detection
            [M, Q] = community_louvain(W,gamma,[],B); 
        end
    [Ppos(:,3) Pneg(:,3)] = participation_coef_sign(W,M);
    [Gpos1(:,3) Gneg1(:,3)] = gateway_coef_sign(W,M,1); % gateway coefficient with node strength centrality
    [Gpos2(:,3) Gneg2(:,3)] = gateway_coef_sign(W,M,2); % gateway coefficient with betweenness centrality
    [Hpos(:,3) Hneg(:,3)] = diversity_coef_sign(W,M); % diversity coefficient
    
    %% calculate summary statistics for participation coefficients per subject
    participation(i_subject,1:3) = mean(Ppos,1);
    participation(i_subject,4:6) = mean(Pneg,1);
    gateway(i_subject,1:3) = mean(Gpos1,1);
    gateway(i_subject,4:6) = mean(Gneg1,1);
    gateway(i_subject,7:9) = mean(Gpos2,1);
    gateway(i_subject,10:12) = mean(Gneg2,1);
    diversity(i_subject,1:3) = mean(Hpos,1);
    diversity(i_subject,4:6) = mean(Hneg,1);
    
    participationSD(i_subject,1:3) = std(Ppos,1);
    participationSD(i_subject,4:6) = std(Pneg,1);
    gatewaySD(i_subject,1:3) = std(Gpos1,1);
    gatewaySD(i_subject,4:6) = std(Gneg1,1);
    gatewaySD(i_subject,7:9) = std(Gpos2,1);
    gatewaySD(i_subject,10:12) = std(Gneg2,1);
    diversitySD(i_subject,1:3) = std(Hpos,1);
    diversitySD(i_subject,4:6) = std(Hneg,1);
    
    %% write results for participation coefficients per subject
    if subjectFiles == 1
        rowNames = cell2table(remainingNodes);

        filenameP = ['subject' num2str(subjectID) '_participation_coefficients.xlsx'];
        pParticipationTable = array2table(Ppos);
        nParticipationTable = array2table(Pneg);
        participationTable = horzcat(rowNames,pParticipationTable,nParticipationTable);
        participationTable.Properties.VariableNames = {'Emotion','PosOrig','PosT1','PosT2','NegOrig','NegT1','NegT2'};
        writetable(participationTable,filenameP);

        filenameG = ['subject' num2str(subjectID) '_gateway_coefficients.xlsx'];
        pGateway1Table = array2table(Gpos1);
        nGateway1Table = array2table(Gneg1);
        pGateway2Table = array2table(Gpos2);
        nGateway2Table = array2table(Gneg2);
        gatewayTable = horzcat(rowNames,pGateway1Table,nGateway1Table,pGateway2Table,nGateway2Table);
        gatewayTable.Properties.VariableNames = {'Emotion','Pos1Orig','Pos1T1','Pos1T2','Neg1Orig','Neg1T1','Neg1T2',...
            'Pos2Orig','Pos2T1','Pos2T2','Neg2Orig','Neg2T1','Neg2T2'};
        writetable(gatewayTable,filenameG);

        filenameH = ['subject' num2str(subjectID) '_diversity_coefficients.xlsx'];
        pDiversityTable = array2table(Hpos);
        nDiversityTable = array2table(Hneg);
        diversityTable = horzcat(rowNames,pDiversityTable,nDiversityTable);
        diversityTable.Properties.VariableNames = {'Emotion','PosOrig','PosT1','PosT2','NegOrig','NegT1','NegT2'};
        writetable(diversityTable,filenameH);
    end
    
    %% clear variables
    clear remainingNodes Ppos Pneg Gpos1 Gneg1 Gpos2 Gneg2 Hpos Hneg;
end

%% write results for participation coefficients across subjects
participation = horzcat(subjectIDlist,participation);
colNames = {'SubjectID','PosOrig','PosT1','PosT2','NegOrig','NegT1','NegT2'};
participationTableAll = array2table(participation,'VariableNames',colNames);
writetable(participationTableAll,'participation_coefficients_distribution.xlsx');
filename = ['dataSet' num2str(dataSet) '_participation.mat'];
save(filename,'participation');

gateway = horzcat(subjectIDlist,gateway);
colNames = {'SubjectID','Pos1Orig','Pos1T1','Pos1T2','Neg1Orig','Neg1T1','Neg1T2',...
        'Pos2Orig','Pos2T1','Pos2T2','Neg2Orig','Neg2T1','Neg2T2'};
gatewayTableAll = array2table(gateway,'VariableNames',colNames);
writetable(gatewayTableAll,'gateway_coefficients_distribution.xlsx');
filename = ['dataSet' num2str(dataSet) '_gateway.mat'];
save(filename','gateway');

diversity = horzcat(subjectIDlist,diversity);
colNames = {'SubjectID','PosOrig','PosT1','PosT2','NegOrig','NegT1','NegT2'};
diversityTableAll = array2table(diversity,'VariableNames',colNames);
writetable(diversityTableAll,'diversity_coefficients_distribution.xlsx');
filename = ['dataSet' num2str(dataSet) '_diversity.mat'];
save(filename,'diversity');

participationSD = horzcat(subjectIDlist,participationSD);
colNames = {'SubjectID','PosOrig','PosT1','PosT2','NegOrig','NegT1','NegT2'};
participationSDTableAll = array2table(participationSD,'VariableNames',colNames);
writetable(participationSDTableAll,'participation_coefficients_SD_distribution.xlsx');
filename = ['dataSet' num2str(dataSet) '_participation_SD.mat'];
save(filename,'participationSD');

gatewaySD = horzcat(subjectIDlist,gatewaySD);
colNames = {'SubjectID','Pos1Orig','Pos1T1','Pos1T2','Neg1Orig','Neg1T1','Neg1T2',...
        'Pos2Orig','Pos2T1','Pos2T2','Neg2Orig','Neg2T1','Neg2T2'};
gatewaySDTableAll = array2table(gatewaySD,'VariableNames',colNames);
writetable(gatewaySDTableAll,'gateway_coefficients_SD_distribution.xlsx');
filename = ['dataSet' num2str(dataSet) '_gateway_SD.mat'];
save(filename','gatewaySD');

diversitySD = horzcat(subjectIDlist,diversitySD);
colNames = {'SubjectID','PosOrig','PosT1','PosT2','NegOrig','NegT1','NegT2'};
diversitySDTableAll = array2table(diversitySD,'VariableNames',colNames);
writetable(diversitySDTableAll,'diversity_coefficients_SD_distribution.xlsx');
filename = ['dataSet' num2str(dataSet) '_diversity_SD.mat'];
save(filename,'diversitySD');

%% plot figures for participation distributions across subjects
figure;
histogram1 = histogram(participation(:,2));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Positive participation coefficient in original data');
saveas(histogram1,'Positive_participation_distribution_original','tiff');

figure;
histogram2 = histogram(participation(:,3));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Positive participation coefficient at first threshold');
saveas(histogram2,'Positive_participation_distribution_threshold1','tiff');

figure;
histogram3 = histogram(participation(:,4));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Positive participation coefficient at second threshold');
saveas(histogram3,'Positive_participation_distribution_threshold2','tiff');

figure;
histogram4 = histogram(participation(:,5));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Negative participation coefficient in original data');
saveas(histogram4,'Negative_participation_distribution_original','tiff');

figure;
histogram5 = histogram(participation(:,6));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Negative participation coefficient at first threshold');
saveas(histogram5,'Negative_participation_distribution_threshold1','tiff');

figure;
histogram6 = histogram(participation(:,7));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Negative participation coefficient at second threshold');
saveas(histogram6,'Negative_participation_distribution_threshold2','tiff');

close all;

%% plot figures for diversity distributions across subjects
figure;
histogram7 = histogram(diversity(:,2));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Positive diversity coefficient in original data');
saveas(histogram7,'Positive_diversity_distribution_original','tiff');

figure;
histogram8 = histogram(diversity(:,3));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Positive diversity coefficient at first threshold');
saveas(histogram8,'Positive_diversity_distribution_threshold1','tiff');

figure;
histogram9 = histogram(diversity(:,4));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Positive diversity coefficient at second threshold');
saveas(histogram9,'Positive_diversity_distribution_threshold2','tiff');

figure;
histogram10 = histogram(diversity(:,5));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Negative diversity coefficient in original data');
saveas(histogram10,'Negative_diversity_distribution_original','tiff');

figure;
histogram11 = histogram(diversity(:,6));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Negative diversity coefficient at first threshold');
saveas(histogram11,'Negative_diversity_distribution_threshold1','tiff');

figure;
histogram12 = histogram(diversity(:,7));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Negative diversity coefficient at second threshold');
saveas(histogram12,'Negative_diversity_distribution_threshold2','tiff');

close all;

%% plot figures for gateway distributions across subjects
figure;
histogram13 = histogram(gateway(:,2));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Positive gateway coefficient (node strength) in original data');
saveas(histogram13,'Positive_gateway_node_distribution_original','tiff');

figure;
histogram14 = histogram(gateway(:,3));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Positive gateway coefficient (node strength) at first threshold');
saveas(histogram14,'Positive_gateway_node_distribution_threshold1','tiff');

figure;
histogram15 = histogram(gateway(:,4));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Positive gateway coefficient (node strength) at second threshold');
saveas(histogram15,'Positive_gateway_node_distribution_threshold2','tiff');

figure;
histogram16 = histogram(gateway(:,5));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Negative gateway coefficient (node strength) in original data');
saveas(histogram16,'Negative_gateway_node_distribution_original','tiff');

figure;
histogram17 = histogram(gateway(:,6));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Negative gateway coefficient (node strength) at first threshold');
saveas(histogram17,'Negative_gateway_node_distribution_threshold1','tiff');

figure;
histogram18 = histogram(gateway(:,7));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Negative gateway coefficient (node strength) at second threshold');
saveas(histogram18,'Negative_gateway_node_distribution_threshold2','tiff');

figure;
histogram19 = histogram(gateway(:,8));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Positive gateway coefficient (betweennness) in original data');
saveas(histogram19,'Positive_gateway_between_distribution_original','tiff');

figure;
histogram20 = histogram(gateway(:,9));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Positive gateway coefficient (betweenness) at first threshold');
saveas(histogram20,'Positive_gateway_between_distribution_threshold1','tiff');

figure;
histogram21 = histogram(gateway(:,10));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Positive gateway coefficient (betweenness) at second threshold');
saveas(histogram21,'Positive_gateway_between_distribution_threshold2','tiff');

figure;
histogram22 = histogram(gateway(:,11));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Negative gateway coefficient (betweenness) in original data');
saveas(histogram22,'Negative_gateway_between_distribution_original','tiff');

figure;
histogram23 = histogram(gateway(:,12));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Negative gateway coefficient (betweenness) at first threshold');
saveas(histogram23,'Negative_gateway_between_distribution_threshold1','tiff');

figure;
histogram24 = histogram(gateway(:,13));
xlim([0 1]);
ylim([1 length(subjectIDlist)]);
title('Negative gateway coefficient (betweenness) at second threshold');
saveas(histogram24,'Negative_gateway_between_distribution_threshold2','tiff');

close all;
