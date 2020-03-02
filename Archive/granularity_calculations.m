% This script calculates tradiational, ICC-based granularity measures for
% all subjects in the specified data set

clear;
clc;

%% specify dataset
dataSet = 18; % set to 18, 39, or 88 and the rest will take care of itself!
print = 1; % set to 0 to turn of spreadsheet generation

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
elseif dataSet == 88
    startRatingsat1 = 0;
    maximumRating = 6;
    skipSubjectRange = 500; % set to invalid subject ID range
    subjectIDlist(subjectIDlist >= skipSubjectRange,:) = [];
else
    startRatingsat1 = 0;
    maximumRating = 6;
end

%% initialize variables
rawICC = zeros(length(subjectIDlist),7);
gran = zeros(length(subjectIDlist),7);
rtoZ = zeros(length(subjectIDlist),7);
zInv = zeros(length(subjectIDlist),7);
corrMatrix = zeros(length(wordList),length(wordList),length(subjectIDlist));
deletedCorr = nan(length(subjectIDlist),length(wordList));
deletedNodes = cell(length(subjectIDlist),1);

%% set valence and arousal categories for sampled words
for i_word = 1:height(words) % define valence categories
if words.Valence(i_word) > 5 % derived based on the database mean for Warriner et al (2013)
    valCat(i_word) = {'Positive'};
    positive(i_word) = 1;
else
    valCat(i_word) = {'Negative'};
    positive(i_word) = 0;
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

%% grab data for each subject and run through calculations
for i_subject = 1:length(subjectIDlist)
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
    %% compute ICCs - positive, negative, valence average; combinations of valence x arousal
    rawICC(i_subject,1) = ICC(subjectData(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
    rawICC(i_subject,2) = ICC(subjectData(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
    rawICC(i_subject,3) = (rawICC(i_subject,1)+rawICC(i_subject,2))/2; % valence mean ICC
    rawICC(i_subject,4) = ICC(subjectData(:,(labels(:,1)==0 & labels(:,2)==0)),'A-k'); % negative valence, low arousal ICC
    rawICC(i_subject,5) = ICC(subjectData(:,(labels(:,1)==1 & labels(:,2)==1)),'A-k'); % positive valence, high arousal ICC
    rawICC(i_subject,6) = ICC(subjectData(:,(labels(:,1)==0 & labels(:,2)==1)),'A-k'); % negative valence, high arousal ICC
    rawICC(i_subject,7) = ICC(subjectData(:,(labels(:,1)==1 & labels(:,2)==0)),'A-k'); % positive valence, low arousal ICC
    %% calculate granularity from ICCs
    gran(i_subject,1) = 1-rawICC(i_subject,1);
    gran(i_subject,2) = 1-rawICC(i_subject,2);
    gran(i_subject,3) = 1-rawICC(i_subject,3);
    gran(i_subject,4) = 1-rawICC(i_subject,4);
    gran(i_subject,5) = 1-rawICC(i_subject,5);
    gran(i_subject,6) = 1-rawICC(i_subject,6);
    gran(i_subject,7) = 1-rawICC(i_subject,7);
    %% Fisher transform ICCs (noted as zICCs)
    rtoZ(i_subject,1) = 0.5*log((1+rawICC(i_subject,1))/(1-rawICC(i_subject,1)));
    rtoZ(i_subject,2) = 0.5*log((1+rawICC(i_subject,2))/(1-rawICC(i_subject,2)));
    rtoZ(i_subject,3) = 0.5*log((1+rawICC(i_subject,3))/(1-rawICC(i_subject,3)));
    rtoZ(i_subject,4) = 0.5*log((1+rawICC(i_subject,4))/(1-rawICC(i_subject,4)));
    rtoZ(i_subject,5) = 0.5*log((1+rawICC(i_subject,5))/(1-rawICC(i_subject,5)));
    rtoZ(i_subject,6) = 0.5*log((1+rawICC(i_subject,6))/(1-rawICC(i_subject,6)));
    rtoZ(i_subject,7) = 0.5*log((1+rawICC(i_subject,7))/(1-rawICC(i_subject,7)));
    %% invert zICCs for intuitive directionality
    zInv(i_subject,1) = rtoZ(i_subject,1)*-1;
    zInv(i_subject,2) = rtoZ(i_subject,2)*-1;
    zInv(i_subject,3) = rtoZ(i_subject,3)*-1;
    zInv(i_subject,4) = rtoZ(i_subject,4)*-1;
    zInv(i_subject,5) = rtoZ(i_subject,5)*-1;
    zInv(i_subject,6) = rtoZ(i_subject,6)*-1;
    zInv(i_subject,7) = rtoZ(i_subject,7)*-1;
end

%% create (and optionally write) tables of results
rawICCs = horzcat(subjectIDlist,rawICC);
rawICCs = array2table(rawICCs,'RowNames',cellstr(num2str(subjectIDlist)),'VariableNames',{'PPID','NegV','PosV','MeanV','NegV_LowA','PosV_HighA','NegV_HighA','PosV_LowA'});
granularity = horzcat(subjectIDlist,gran);
granularity = array2table(granularity,'RowNames',cellstr(num2str(subjectIDlist)),'VariableNames',{'PPID','NegV','PosV','MeanV','NegV_LowA','PosV_HighA','NegV_HighA','PosV_LowA'});
zICCs = horzcat(subjectIDlist,rtoZ);
zICCs = array2table(rtoZ,'RowNames',cellstr(num2str(subjectIDlist)),'VariableNames',{'PPID','NegV','PosV','MeanV','NegV_LowA','PosV_HighA','NegV_HighA','PosV_LowA'});
zInvs = horzcat(subjectIDlist,zInv);
zInvs = array2table(zInvs,'RowNames',cellstr(num2str(subjectIDlist)),'VariableNames',{'PPID','NegV','PosV','MeanV','NegV_LowA','PosV_HighA','NegV_HighA','PosV_LowA'});

if print == 1
    %writetable(rawICCs,['dataSet' num2str(dataSet) '_raw_ICC_values.xlsx']);
    writetable(granularity,['dataSet' num2str(dataSet) '_granularity_values.xlsx']);
    %writetable(zICCs,['dataSet' num2str(dataSet) '_Fisher-z_transformed_ICC_values.xlsx']);
    writetable(zInvs,['dataSet' num2str(dataSet) '_inverted_Fisher-z_granularity.xlsx']);
end

%% save off matrix for construct validity checks
filename = ['dataSet' num2str(dataSet) '_zICC.mat'];
save(filename,'rtoZ');