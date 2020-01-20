% This script calculates traditional, ICC-based granularity measures for
% all subjects and every day in the specified data set

clear;
clc;

%% specify dataset
dataSet = 18; % only 18 or 39 for now
imputeData = 1; % set to 1 to impute values for missing days
print = 0; % set to 1 to print data to file
saveData = 0; % set to 1 to save measure-specific data matrices

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
    index1 = find(allData(:,1)==subjectID);
    subjectData = allData(index1,2:end);
    dayIDlist = unique(subjectData(:,1));
    for i_day = 1:length(dayIDlist)
        dayData = [];
        dayID = dayIDlist(i_day);
        index2 = find(subjectData(:,1)==dayID);
        dayData = subjectData(index2,4:end);
        %% remove invalid values and missing data
        dayData(dayData > maximumRating) = NaN;
        missingData = isnan(dayData); 
        missingData2 = any(missingData,2); 
        dayData = dayData(~missingData2,:); 
        %% if necessary, rescale data to start at 0
        if startRatingsat1 == 1
            dayData = dayData-1;
        end
        %% compute ICCs - positive, negative, valence average; combinations of valence x arousal
        rawICC(i_day,1) = ICC(dayData(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
        rawICC(i_day,2) = ICC(dayData(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
        rawICC(i_day,3) = (rawICC(i_day,1)+rawICC(i_day,2))/2; % valence mean ICC
        rawICC(i_day,4) = ICC(dayData(:,(labels(:,1)==0 & labels(:,2)==0)),'A-k'); % negative valence, low arousal ICC
        rawICC(i_day,5) = ICC(dayData(:,(labels(:,1)==1 & labels(:,2)==1)),'A-k'); % positive valence, high arousal ICC
        rawICC(i_day,6) = ICC(dayData(:,(labels(:,1)==0 & labels(:,2)==1)),'A-k'); % negative valence, high arousal ICC
        rawICC(i_day,7) = ICC(dayData(:,(labels(:,1)==1 & labels(:,2)==0)),'A-k'); % positive valence, low arousal ICC
        %% set lower bound of ICCs to 0 (no negative values)
        rawICC(rawICC<0) = 0;
        %% calculate granularity from ICCs
        gran(i_day,1) = 1-rawICC(i_day,1);
        gran(i_day,2) = 1-rawICC(i_day,2);
        gran(i_day,3) = 1-rawICC(i_day,3);
        gran(i_day,4) = 1-rawICC(i_day,4);
        gran(i_day,5) = 1-rawICC(i_day,5);
        gran(i_day,6) = 1-rawICC(i_day,6);
        gran(i_day,7) = 1-rawICC(i_day,7);
        %% Fisher transform ICCs (noted as zICCs)
        rtoZ(i_day,1) = 0.5*log((1+rawICC(i_day,1))/(1-rawICC(i_day,1)));
        rtoZ(i_day,2) = 0.5*log((1+rawICC(i_day,2))/(1-rawICC(i_day,2)));
        rtoZ(i_day,3) = 0.5*log((1+rawICC(i_day,3))/(1-rawICC(i_day,3)));
        rtoZ(i_day,4) = 0.5*log((1+rawICC(i_day,4))/(1-rawICC(i_day,4)));
        rtoZ(i_day,5) = 0.5*log((1+rawICC(i_day,5))/(1-rawICC(i_day,5)));
        rtoZ(i_day,6) = 0.5*log((1+rawICC(i_day,6))/(1-rawICC(i_day,6)));
        rtoZ(i_day,7) = 0.5*log((1+rawICC(i_day,7))/(1-rawICC(i_day,7)));
        %% invert zICCs for intuitive directionality
        zInv(i_day,1) = rtoZ(i_day,1)*-1;
        zInv(i_day,2) = rtoZ(i_day,2)*-1;
        zInv(i_day,3) = rtoZ(i_day,3)*-1;
        zInv(i_day,4) = rtoZ(i_day,4)*-1;
        zInv(i_day,5) = rtoZ(i_day,5)*-1;
        zInv(i_day,6) = rtoZ(i_day,6)*-1;
        zInv(i_day,7) = rtoZ(i_day,7)*-1;
    end
    %% impute missing data, if desired
    if imputeData == 1
        gran = fillmissing(gran,'movmean',3,1);
        zInv = fillmissing(zInv,'movmean',3,1);
    end
    %% add subject results to summary table
    ppID = ['PP' num2str(subjectID)];
    if i_subject == 1
        day = (1:1:max(allData(:,2)))';
        day = array2table(day,'VariableNames',{'Day'});
                
        gran_Neg = tablecompile_start(day,dayIDlist,gran(:,1),ppID,'Day','ppDay','ppDay');
        gran_Pos = tablecompile_start(day,dayIDlist,gran(:,2),ppID,'Day','ppDay','ppDay');
        gran_M = tablecompile_start(day,dayIDlist,gran(:,3),ppID,'Day','ppDay','ppDay');
        
        zInv_Neg = tablecompile_start(day,dayIDlist,zInv(:,1),ppID,'Day','ppDay','ppDay');
        zInv_Pos = tablecompile_start(day,dayIDlist,zInv(:,2),ppID,'Day','ppDay','ppDay');
        zInv_M = tablecompile_start(day,dayIDlist,zInv(:,3),ppID,'Day','ppDay','ppDay');
    else
        gran_Neg = tablecompile_iter(gran_Neg,dayIDlist,gran(:,1),ppID,'Day','ppDay','ppDay');
        gran_Pos = tablecompile_iter(gran_Pos,dayIDlist,gran(:,2),ppID,'Day','ppDay','ppDay');
        gran_M = tablecompile_iter(gran_M,dayIDlist,gran(:,3),ppID,'Day','ppDay','ppDay');
        
        zInv_Neg = tablecompile_iter(zInv_Neg,dayIDlist,zInv(:,1),ppID,'Day','ppDay','ppDay');
        zInv_Pos = tablecompile_iter(zInv_Pos,dayIDlist,zInv(:,2),ppID,'Day','ppDay','ppDay');
        zInv_M = tablecompile_iter(zInv_M,dayIDlist,zInv(:,3),ppID,'Day','ppDay','ppDay');
    end
    clear rawICC gran rtoZ zInv
end

%% write tables to file
if print == 1
    writetable(gran_Neg,['dataSet' num2str(dataSet) '_gran_Neg_daily.xlsx']);
    writetable(gran_Pos,['dataSet' num2str(dataSet) '_gran_Pos_daily.xlsx']);
    writetable(gran_M,['dataSet' num2str(dataSet) '_gran_M_daily.xlsx']);

    writetable(zInv_Neg,['dataSet' num2str(dataSet) '_zInv_Neg_daily.xlsx']);
    writetable(zInv_Pos,['dataSet' num2str(dataSet) '_zInv_Pos_daily.xlsx']);
    writetable(zInv_M,['dataSet' num2str(dataSet) '_zInv_M_daily.xlsx']);
end

%% save tables as matrices
gran_Neg_array = table2array(gran_Neg)';
gran_Pos_array = table2array(gran_Pos)';
gran_M_array = table2array(gran_M)';
zInv_Neg_array = table2array(zInv_Neg)';
zInv_Pos_array = table2array(zInv_Pos)';
zInv_M_array = table2array(zInv_M)';

% impute missing data, if desired
if imputeData == 1
    gran_Neg_array(2:end,:) = fillmissing(gran_Neg_array(2:end,:),'movmean',3,2);
    gran_Pos_array(2:end,:) = fillmissing(gran_Pos_array(2:end,:),'movmean',3,2);
    gran_M_array(2:end,:) = fillmissing(gran_M_array(2:end,:),'movmean',3,2);
    zInv_Neg_array(2:end,:) = fillmissing(zInv_Neg_array(2:end,:),'movmean',3,2);
    zInv_Pos_array(2:end,:) = fillmissing(zInv_Pos_array(2:end,:),'movmean',3,2);
    zInv_M_array(2:end,:) = fillmissing(zInv_M_array(2:end,:),'movmean',3,2);
end

if saveData == 1
    subjectCol = [0; subjectIDlist];

    gran_Neg_array = horzcat(subjectCol,gran_Neg_array);
    save(['dataSet' num2str(dataSet) '_gran_Neg_daily.mat'],'gran_Neg_array');

    gran_Pos_array = horzcat(subjectCol,gran_Pos_array);
    save(['dataSet' num2str(dataSet) '_gran_Pos_daily.mat'],'gran_Pos_array');

    gran_M_array = horzcat(subjectCol,gran_M_array);
    save(['dataSet' num2str(dataSet) '_gran_M_daily.mat'],'gran_M_array');
    
    zInv_Neg_array = horzcat(subjectCol,zInv_Neg_array);
    save(['dataSet' num2str(dataSet) '_zInv_Neg_daily.mat'],'zInv_Neg_array');

    zInv_Pos_array = horzcat(subjectCol,zInv_Pos_array);
    save(['dataSet' num2str(dataSet) '_zInv_Pos_daily.mat'],'zInv_Pos_array');

    zInv_M_array = horzcat(subjectCol,zInv_M_array);
    save(['dataSet' num2str(dataSet) '_zInv_M_daily.mat'],'zInv_M_array');
end