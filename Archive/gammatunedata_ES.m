%% This is a script for looking at the impact of gamma values on the number of communities detected 

% Data Needed: dataname.dat and wordslist input files must be in the same
% folder as this one when you run this script.

% Toolboxes: Toolbox for community_louvain algorithm should be included.

% Set values in following 'General setup' section before running script

%% General setup

clear

% If computing instance network, see instance variable to 1
compinstance = 1;

% for thresholding the matrix, set the proportion
proportion= 1;

%set datafile name
filename = '39ES_DATA_CLEAN.dat';

%set wordlist filename
wordlistname = 'wordslist_39.mat';

%gamma tune parameters
g_start = 1;
%g_end= 0; %variable is set in the script based on the max number of
%communities possible.
increment_g=.01;

%% Import raw data & wordlist, setup files for ICC

rawdata = importdata(filename);
subjectIDList = unique(rawdata.data(:,1));
num_subs= length(subjectIDList);

load(wordlistname, 'words','M')
ICC_mat = zeros(length(subjectIDList),4);

% [M, I1] = sort(words);
% [B, I2] = sort(M(:,1));
% [~, I3] = sort(I1);
% 
% strcmp(M, B);
% M = M(I2,:);
% M = M(I3,:);

labels = cellfun(@conv_label,M(:,2:end));

%% ICC Computation (Needed for tuning the network metrics)

for i_subject = 1:length(subjectIDList)
    data = [];
    subjectID = subjectIDList(i_subject);
    index = find(rawdata.data(:,1)==subjectID);
    data = rawdata.data(index,2:end);
    score = data(:,2:end);
    %% Removing missing data
    missing_data = isnan(score);
    missing_data2 = any(missing_data,2);
    score = score(~missing_data2,:);
    
    %% Computation of ICC
    ICC_mat(i_subject,1) = ICC(score(:,(labels(:,1)==0)),'A-k');
    ICC_mat(i_subject,2) = ICC(score(:,(labels(:,1)==1)),'A-k');
    ICC_mat(i_subject,3) = ICC(score(:,(labels(:,2)==0)),'A-k');
    ICC_mat(i_subject,4) = ICC(score(:,(labels(:,2)==1)),'A-k');
end

%% Fisher transform, Compute global zICC
cICC_mat = num2cell(ICC_mat);
zICC_mat = cellfun(@r2z,cICC_mat);

posnegICC=(zICC_mat(:,1)+zICC_mat(:,2))/2;

%% Save ICC values to a datatable
%ICCTable = array2table(ICC_mat,'RowNames',cellstr(num2str(subjectIDList)),'VariableNames',{'Neg_V',...
%     'Pos_V','Low_A','High_A'});
% writetable(ICCTable,'ICC.dat','WriteRowNames',true);

%% Initialization
id = [];
% modularity = zeros(num_subs,1);
% Ctot_posmat = zeros(num_subs,1);
% Sposmat = zeros(num_subs,1);
% Snegmat = zeros(num_subs,1);
% vposmat = zeros(num_subs,1);
% vnegmat = zeros(num_subs,1);
% pos_diversity = zeros(num_subs,1);
% neg_diversity = zeros(num_subs,1);
% pos_gateway = zeros(num_subs,1);
% neg_gateway = zeros(num_subs,1);
% pos_Participation = zeros(num_subs,1);
% neg_Participation = zeros(num_subs,1);

for subjectID = 1:84  % read data from one subject each time
    [X,words] = getDataFor1Subject(rawdata, subjectID,-1);
    % ignore IDs with no data points  
    if isempty(X)
        continue
    end
    
    % Invert matrix if computing instance networks- CHECK THIS CODE
    if compinstance == 1
        X=transpose(X);
        X = X(:,all(~isnan(X)));
    else
        X=X(all(~isnan(X),2),:);
    end
    
    
    % compute Pearson correlation
    W = corr(X,'rows','pairwise');
    % delete the NaN columns and rows in W
    S_first_row = W(1,:);
    deleteTime = find(isnan(S_first_row)); %find the column or row that is NaN
    W(deleteTime, :) = [];
    W(:, deleteTime) = [];
    if isempty(W)
        continue
    end
    
    id = [id subjectID];%record the remianing subject ids after omitting
    
    
%% Matrix adjustment
    
% Set diagnonals to 0 (Necessary for BCT functions to run properly)
W(logical(eye(size(W))))=0; 

% Now normalize lengths to weights -- should remove overall space usage as
% factor. This is also necessary to compute clustering coefficient measures
% on the matrix (see BCT documentation).
%W= weight_conversion(threshold_proportional(W,proportion), 'normalize');

xcount=0; %set x-count on each loop

while e_gamma=
for x=1:0.01:1.5 %set steps for gamma values (can be modified to look at different ranges)

    W= W;   % set subject matrix
    n  = size(W,1);          % number of nodes
    M  = 1:n;                % initial community affiliations
    Q0 = -1; Q1 = 0;         % initialize modularity values

    %Community detection fune-tuning loop
    while Q1-Q0>1e-5           % while modularity increases
        Q0 = Q1;
        [M, Q1] = community_louvain(W,x,[],'negative_sym');
    end

%Data storing
maxcom(x) = max(M);
modularity(x) = Q1';

xcount = xcount+1; %add .01 to x for next loop
allmaxcom(:,xcount) = maxcom; %save off maxcom data for gamma value
allmodularity(:,xcount) = modularity; %save off modularity data for gamma value
end

clear W n M Q0 Q1 maxcom modularity x xcount

end

%% Create Data Table of Results

meanmaxcom=mean(allmaxcom(:,:))';
sdmaxcom=std(allmaxcom(:,:))';
meanmodularity=mean(allmodularity(:,:))';
sdmodularity=std(allmodularity(:,:))';

gammatune=table(meanmaxcom, sdmaxcom, meanmodularity, sdmodularity, z');
gammatune.Properties.VariableNames = {'Mean_Num_Communities' 'SD_Num_Communities' 'Mean_Modularity' 'SD_Modularity' 'Gamma'};

clearvars -except gammatune allmaxcom allmodularity ICCtable

%% Join files and calculate correlations across levels of gamma

tune_assess_table=join(gammatune,ICCTable); %Does the ICC table have IDs in it? Should perhaps only make table with subset of ICC



%% PLOTTING

z=1:0.01:1.5;

% generate plot for maxcom with SD error bars
errorbar(mean(allmaxcom(:,:)),z,std(allmaxcom(:,:))); %plot with SD
title ('Mean # of Communities (+/- SD) by Gamma Value')


% generate plot for maxcom with SEM error bars
figure;
sem=std(allmaxcom(:,:))/sqrt(length(allmaxcom(:,:)));
errorbar(mean(allmaxcom(:,:)),z,sem);
title ('Mean # of Communities (+/- SEM) by Gamma Value')

% generate plot for maxcom with SD error bars
errorbar(mean(allmodularity(:,:)),z,std(allmodularity(:,:))); %plot with SD
title ('Mean Modularity (+/- SD) by Gamma Value')

% generate plot for maxcom with SEM error bars
figure;
sem=std(allmodularity(:,:))/sqrt(length(allmodularity(:,:)));
errorbar(mean(allmodularity(:,:)),z,sem);
title ('Mean Modularity (+/- SEM) by Gamma Value')


%% 

function [X,words] = getDataFor1Subject(rawdata,subjectID,numOfShownWord)
% Input:
%     1) subject ID: pick a num in (1~84), or -1 for all subjects
%     2) numOfShownWord : pick a num in (1~39), or -1 for all words
% Output: X, words

    %% choose the num of words you want to show
    if numOfShownWord == -1   % for all words
        numOfShownWord = 39;
    end
    numOfShownWord = numOfShownWord +2;
    words = rawdata.textdata(3:numOfShownWord);
    
    %% choose the subject ID and days
    if subjectID == -1  % for all subjects
        X = rawdata.data(:,3:numOfShownWord); % for all the data
    else     
        subjectIndex = subjectID == rawdata.data(:,1);
        X = rawdata.data(subjectIndex,3:numOfShownWord);
    end

end


%%
function y = conv_label(x)
switch x
    case {'Negative','Low'}
        y = 0;
    case {'Positive','High'}
        y = 1;
    otherwise
        y = 4;
end
end

%%

function z = r2z(r)
    z = 0.5*log((1+r)/(1-r));
end
