clear;
clc;

%% load data file, set tuning parameter, initialize variables
%datafile = 'ARIEOD_n17_040718.xlsx';
datafile = 'ES39datav2.xlsx';
rawdata = importdata(datafile);
subjectIDlist = unique(rawdata.data(:,1)); % grab subject IDs from first column
wordlist = rawdata.colheaders(2:end)';  % grab sampled words from top row
rawICC_mat = zeros(length(subjectIDlist),3);
rtoZ_mat = zeros(length(subjectIDlist),3);
gran_mat = zeros(length(subjectIDlist),3);
granFZ_mat = zeros(length(subjectIDlist),3);
corrMatrix = zeros(length(wordlist),length(wordlist),length(subjectIDlist));
deletedCorr = nan(length(subjectIDlist),length(wordlist));
%deletedCorr = zeros(1,length(subjectIDlist));
dropSubject = nan(1,length(subjectIDlist));
gMin = 0; % set starting value for tuning parameter gamma
gMax = 3; % set ending value for tuning parameter gamma
gStepSize = .05;  % set step size for tuning parameter gamma
gArray = [gMin:gStepSize:gMax];
nGamma = length(gArray);
nodes = nan(length(subjectIDlist),1);
modularity = nan(length(subjectIDlist),nGamma);
community = nan(length(subjectIDlist),nGamma);

%% set valence and arousal categories for sampled words
%wordfile = 'words18.csv'; % load word list that includes raw norms
wordfile = 'words39v2.csv';
words = readtable(wordfile); 

for i_word=1:height(words) % define valence categories
if words.Valence(i_word)>5
    valCat(i_word)={'Positive'};
    positive(i_word)=1;
else
    valCat(i_word)={'Negative'};
    positive(i_word)=0;
end
end 

for i_word=1:height(words) % define arousal categories
if words.Arousal(i_word)>4.6
    aroCat(i_word)={'High'};
    high(i_word)=1;
else
    aroCat(i_word)={'Low'};
    high(i_word)=0;
end
end 

words=[words valCat' aroCat']; % append table with category assignments
words.Properties.VariableNames(5:end)={'ValCat' 'AroCat'}; % label new variables
labels=[positive' high']; % create matrix for logical indexing in ICC commands

%% grab data for each subject and run through calculations
for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index = find(rawdata.data(:,1)==subjectID);
    subjectData = rawdata.data(index,2:end);
    %score = data(:,:); % removed; unnecessary given updated source data file
    %% remove missing data
    missing_data = isnan(subjectData); % edited 'score' to 'data'
    missing_data2 = any(missing_data,2); % edited 'score' to 'data'
    subjectData = subjectData(~missing_data2,:); % edited 'score' to 'data'
    %% Rescale to start at 0
    if min(subjectData) > 0;
        subjectData = subjectData-1;
    end
    %% compute ICCs - positive, negative, and valence average
    rawICC_mat(i_subject,1) = ICC(subjectData(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
    rawICC_mat(i_subject,2) = ICC(subjectData(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
    rawICC_mat(i_subject,3) = (rawICC_mat(i_subject,1)+rawICC_mat(i_subject,2))/2; % valence average ICC
    %% Fisher transform ICCs
    rtoZ_mat(i_subject,1) = 0.5*log((1+rawICC_mat(i_subject,1))/(1-rawICC_mat(i_subject,1)));
    rtoZ_mat(i_subject,2) = 0.5*log((1+rawICC_mat(i_subject,2))/(1-rawICC_mat(i_subject,2)));
    rtoZ_mat(i_subject,3) = 0.5*log((1+rawICC_mat(i_subject,3))/(1-rawICC_mat(i_subject,3)));
    %% compute granularity for both raw and Fisher transformed ICCs
    gran_mat(i_subject,1) = 1-rawICC_mat(i_subject,1);
    gran_mat(i_subject,2) = 1-rawICC_mat(i_subject,2);
    gran_mat(i_subject,3) = 1-rawICC_mat(i_subject,3);
    granFZ_mat(i_subject,1) = -1*rtoZ_mat(i_subject,1);
    granFZ_mat(i_subject,2) = -1*rtoZ_mat(i_subject,2);
    granFZ_mat(i_subject,3) = -1*rtoZ_mat(i_subject,3);
    %% create correlation matrix
    subjectMatrix = corrcoef(subjectData); % calculate n-by-n correlation matrix
    subjectMatrix(logical(eye(size(subjectMatrix)))) = 0; % set diagonal to 0 for BCT functions
% Now normalize lengths to weights -- should remove overall space usage as
% factor. This is also necessary to compute clustering coefficient measures
% on the matrix (see BCT documentation).
    %subjectMatrix = weight_conversion(threshold_proportional(subjectMatrix, 1),'normalize');
    %% save off correlation matrix before continuing
    %subjectIDstr = int2str(subjectID);
    %corrName = strcat('corrMatrix',subjectIDstr);
    %save(corrName,'subjectMatrix');
    %load(corrName);
    %ppIDstr=['pp' subjectIDstr];
    %structMatrix = struct(ppIDstr,subjectMatrix);
    corrMatrix(:,:,i_subject) = subjectMatrix;
    %% delete NaN values from correlation matrix
    missing_corr = find(isnan(subjectMatrix(1,:)));
    if numel(missing_corr) >= .2*length(wordlist); % add subject to the cut list
        dropSubject(i_subject) = numel(missing_corr);
    end
    if missing_corr >= 1        % save off which row/col were removed
        %deletedCorr(i_subject) = missing_corr;
        deletedCorr(i_subject,1:numel(missing_corr)) = missing_corr;
        for i_value = 1:length(missing_corr);
            subjectMatrix(missing_corr(i_value),:) = []; % delete missing row
            subjectMatrix(:,missing_corr(i_value)) = []; % delete missing column
            if i_value < length(missing_corr)
                missing_corr(i_value+1) = missing_corr(i_value+1)-i_value; % decrement subsequent value by 1 to account for offset
            end
        end
    else
        deletedCorr(i_subject) = 0;
    end
    %% community detection finetuning loop
    for i_gamma = 1:nGamma
            W = subjectMatrix;  % set subject matrix
            n = size(W,1);      % number of nodes
            M = 1:n;            % initial community affiliations
            Q0 = -1; Q1 = 0;    % initialize modularity values
                while Q1 - Q0 > 1e-5    % while modularity increases
                    Q0 = Q1;            % perform community detection
                    [M, Q1] = community_louvain(W,gArray(i_gamma),[],'negative_sym');
                end
            nodes(i_subject) = n;       % save number of nodes
            %if any(isnan(M))
            %    missing_comm(i_subject) = subjectID;    % drop participants who produce NaN values
            %else
                community(i_subject,i_gamma) = max(M);   % save number of communities
            %end
            %if any(isnan(Q1))
            %    missing_mod(i_subject) = subjectID;         % drop participants who produce NaN values
            %else
                modularity(i_subject,i_gamma) = Q1;         % save modularity (q value)
            %end
            if max(M) == n              % stop if number of communities equals number of nodes
                break;
            end
    end
    %% calculate clustering coefficient
    [C,C_null,Ctot,Ctot_null] = clustering_coef_wu_sign(W,3);
    clusterCoefGlobal(i_subject) = Ctot; % save off global value
    clusterCoefNode{i_subject} = C; % save off node level vector
end

%% calculate distributional statistics on the modularity values
meanModularity = mean(modularity,1,'omitnan');
sdModularity = std(modularity,1,'omitnan');
varModularity = var(modularity,1,'omitnan');
varModMin = min(varModularity);                     % find minimum variance of q
gammaOptimalbyModVariance = gArray(varModularity==varModMin);    % find gamma at which variance of q is minimum

%% calculate distributional statistics on the number of communities
meanCommunity = mean(community,1,'omitnan');
sdCommunity = std(community,1,'omitnan');
varCommunity = var(community,1,'omitnan');
%medCommunity = median(community);
minCommunity = min(community);
maxCommunity = max(community);
rangeCommunity = maxCommunity - minCommunity;

%% run normality tests on distributions of number of communities
%communitySWtest = nan(1,nGamma); % Shapiro-Wilk test of normality
%communityKStest = nan(1,nGamma); % Kolmogorov-Smirnov test of normality
%for g = 1:nGamma
    %communitySWtest(g) = swtest(community(:,g));
    %communityKStest(g) = kstest(community(:,g));
%end

%% create data tables of results
%gammaTuning = table(gArray', meanCommunity', sdCommunity', rangeCommunity', communitySWtest', meanModularity', sdModularity', varModularity');
%gammaTuning.Properties.VariableNames = {'Gamma' 'Mean_Num_Com' 'SD_Num_Com' 'Range_Num_Com' 'SW_Num_Com' 'Mean_Mod' 'SD_Mod' 'Var_Mod'};

rawICCs = array2table(rawICC_mat,'RowNames',cellstr(num2str(subjectIDlist)),'VariableNames',{'Neg_V','Pos_V','Average_V'});

cutList = subjectIDlist(~isnan(dropSubject)); % save off list of subjects who had > 20% missing correlations
clusterCoefGlobal = clusterCoefGlobal';

%% compare modularity and number of communities against ICC

ICCtoModularity = corr(rtoZ_mat(:,3),modularity,'rows','complete');
ICCtoCommunity = corr(rtoZ_mat(:,3),community,'rows','complete');

corrModMax = max(ICCtoModularity);                  % find maximum positive correlation between ICC and modularity (expected direction)
corrModMin = min(ICCtoModularity);                  % find maximum negative correlation between ICC and modularity
gammaOptimalbyModCorrPos = gArray(ICCtoModularity==corrModMax);    % find gamma at which positive correlation is highest
gammaOptimalbyModCorrNeg = gArray(ICCtoModularity==corrModMin);    % find gamma at which negative correlation is highest
corrComMax = max(ICCtoCommunity);                   % find maximum positive correlation between ICC and number of communities
corrComMin = min(ICCtoCommunity);                   % find maximum negative correlation between ICC and number of communities (expected direction)
gammaOptimalbyComCorrPos = gArray(ICCtoCommunity==corrComMax);    % find gamma at which positive correlation is highest
gammaOptimalbyComCorrNeg = gArray(ICCtoCommunity==corrComMin);    % find gamma at which negative correlation is highest

%% generate plots for results
% plot for mean modularity with SD error bars
errorbar(gArray,meanModularity,sdModularity); 
title ('Mean Modularity (+/- SD) by Gamma Value')
xlabel ('Gamma Value')
ylabel ('Modularity')

% plot for mean number of communities with SD error bars
figure;
errorbar(gArray,meanCommunity,sdCommunity); 
title ('Mean Number of Communities (+/- SD) by Gamma Value')
xlabel ('Gamma Value')
ylabel ('Number of Communities')

% plot of correlation between modularity and ICC by gamma value
figure;
plot(gArray,ICCtoModularity);
title ('Correlation between Modularity and ICC by Gamma Value')
xlabel ('Gamma Value')
ylabel ('r Value')
line ([gammaOptimalbyModCorrNeg gammaOptimalbyModCorrNeg], get(gca, 'ylim'), 'Color','r');
line ([gammaOptimalbyModCorrPos gammaOptimalbyModCorrPos], get(gca, 'ylim'), 'Color','r');

% plot of correlation between community and ICC by gamma value
figure;
plot(gArray,ICCtoCommunity);
title ('Correlation between Number of Communities and ICC by Gamma Value')
xlabel ('Gamma Value')
ylabel ('r Value')
line ([gammaOptimalbyComCorrNeg gammaOptimalbyComCorrNeg], get(gca, 'ylim'), 'Color','r');
line ([gammaOptimalbyComCorrPos gammaOptimalbyComCorrPos], get(gca, 'ylim'), 'Color','r');

% plot of correlation between ICC and community, against range in community
figure;
scatter(ICCtoCommunity,rangeCommunity);
title ('Correlation between ICC and Number of Communities, Against Range in Number of Communities')
xlabel ('r Value')
ylabel ('Range in Number of Communities')


% To-do list EOD 22 April 2018:
% - Optimize how/if variables are initialized and values are saved off
% - Figure out a better way to drop/clean bad data along the way
% e.g. clustering coefficients that are 0 could be dropped immediately
% - Frequently getting errors in swtest; is this necessary to keep?

% Predicted relationships:
% More communities = higher granularity (lower ICC); negative correlation
% More communities = lower modularity = lower ICC; positive correlation
% So we want to optimize gamma by ComCorrNeg and/or ModCorrPos
