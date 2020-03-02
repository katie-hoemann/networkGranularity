% NOTE: Currently built out for 3-day weighted sliding window; modification to
% window size or weighting will require structural edits to script

clear;
clc;

%% specify dataset and parameters
dataSet = 18; % only 18 or 39
printSummaries = 0; % set to 1 to save summary tables to file
printTables = 0; % set to 1 to write measure-specific tables to file
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

if startRatingsat1 == 1
    scaleMin = 0;
    scaleMax = maximumRating-1;
else
    scaleMin = 0;
    scaleMax = maximumRating;
end
scaleMidpoint = (scaleMax-scaleMin)/2;
scaleSteps = maximumRating+1;

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
        % calculate granularity (ICCs and zICCs)
        rawICC(i_day,1) = ICC(dayData(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
        rawICC(i_day,2) = ICC(dayData(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
        rawICC(i_day,3) = (rawICC(i_day,1)+rawICC(i_day,2))/2; % valence mean ICC
        rawICC(rawICC<0) = 0;
        rtoZ(i_day,1) = 0.5*log((1+rawICC(i_day,1))/(1-rawICC(i_day,1)));
        rtoZ(i_day,2) = 0.5*log((1+rawICC(i_day,2))/(1-rawICC(i_day,2)));
        rtoZ(i_day,3) = 0.5*log((1+rawICC(i_day,3))/(1-rawICC(i_day,3)));
        % calculate emodiversity
        totalCount = sum(sum(dayData>0));
        for i_emotion = 1:length(wordList)
            countEmotion = sum(dayData(:,i_emotion)>0);
            propEmotion(i_emotion) = (countEmotion./totalCount)*log(countEmotion./totalCount);
        end
        propEmotion = propEmotion(isfinite(propEmotion));
        emodiversity(i_day) = -1*sum(propEmotion);
        % calculate instability
        for i_emotion = 1:length(wordList)
            diffEmotion = diff(dayData(:,i_emotion));
            MSSD_emotion(i_emotion) = mean(diffEmotion.^2);
        end
        instability(i_day) = mean(MSSD_emotion);
        % calculate average emotional intensity
        emotionIntensity(i_day) = mean(mean(dayData));
        % calculate single emotion experience
        [prompts,~] = size(dayData);
        for i_prompt = 1:prompts
            highestEmotion = max(dayData(i_prompt,:));
            distances = highestEmotion-dayData(i_prompt,:);
            distancesClean = distances(distances ~= 0);
            meanDistance(i_prompt) = mean(distancesClean);
        end
        meanDistance(isnan(meanDistance)) = 0; % for instances where all emotions are rated as 0, meanDistance=0    
        singleEmotion(i_day) = mean(meanDistance)/scaleMax; % calculate as ratio to max distance
        % calculate extreme scale use
        midpointDistance = abs(dayData-scaleMidpoint);
        extremeScaleUse(i_day) = mean(mean(midpointDistance))/scaleMidpoint; % calculate as ratio to possible distance from scale midpoint
        % calculate mean and SD for positive and negative valence
        mPositive(i_day) = mean(mean(dayData(:,(valence==1))));
        mNegative(i_day) = mean(mean(dayData(:,(valence==2))));
        sdPositive(i_day) = mean(std(dayData(:,(valence==1))));
        sdNegative(i_day) = mean(std(dayData(:,(valence==2))));
    end
    %% combine measures
    measures(:,1:3) = rtoZ(:,1:3);
    measures(:,4) = emodiversity;
    measures(:,5) = instability;
    measures(:,6) = emotionIntensity;
    measures(:,7) = singleEmotion;
    measures(:,8) = extremeScaleUse;
    measures(:,9) = mPositive;
    measures(:,10) = mNegative;
    measures(:,11) = sdPositive;
    measures(:,12) = sdNegative;
    %% create table of values per day
    ppID = ['PP' num2str(subjectID)];
    day = (1:1:max(allData(:,2)))';
    day = array2table(day,'VariableNames',{'Day'});
    measures = horzcat(dayIDlistUpdate,measures);
    measures = array2table(measures,'VariableNames',{'ppDay','zICC_N','zICC_P','zICC_M','emodiv','instab','intens','singleEmo','extreme','mPos','mNeg','sdPos','sdNeg'});
    measures = outerjoin(day,measures,'LeftKeys','Day','RightKeys','ppDay');
    measures = removevars(measures,'ppDay');
    %% impute any missing data
    measures = fillmissing(measures,'movmean',3);
    %% create windowed values
    measures = table2array(measures(:,2:end));
    day = table2array(day);
    for i_window = 1:length(day)-2
        focusDay = i_window+1;
        windowData = [];
        windowData = [measures(focusDay-1,:); measures(focusDay,:); measures(focusDay+1,:)];
        windowSimpleMeans(i_window,:) = mean(windowData,1); % simple average across days in window
        windowWeightedMeans(i_window,:) = measures(focusDay-1,:)*1/6 + measures(focusDay,:)*4/6 + measures(focusDay+1,:)*1/6; % weighted average across days in window
    end
    %% calculate subject-level mean and SD for measures
    measureSimpleMeans(i_subject,:) = mean(windowSimpleMeans,'omitnan');
    measureSimpleSDs(i_subject,:) = std(windowSimpleMeans,'omitnan');
    measureWeightedMeans(i_subject,:) = mean(windowWeightedMeans,'omitnan');
    measureWeightedSDs(i_subject,:) = std(windowWeightedMeans,'omitnan'); 
    %% add subject results to summary table
    windowSimpleMeans_Table = array2table(windowSimpleMeans);
    windowWeightedMeans_Table = array2table(windowWeightedMeans);
    if i_subject == 1
        window = (1:1:max(allData(:,2))-2)';
        window = array2table(window);
        variableNames = {'Window',ppID};
        
        zICC_Neg_Table_S = tableconcat_start(window,windowSimpleMeans_Table(:,1),variableNames);
        zICC_Pos_Table_S = tableconcat_start(window,windowSimpleMeans_Table(:,2),variableNames);
        zICC_Overall_Table_S = tableconcat_start(window,windowSimpleMeans_Table(:,3),variableNames);
        emodiversity_Table_S = tableconcat_start(window,windowSimpleMeans_Table(:,4),variableNames);
        instability_Table_S = tableconcat_start(window,windowSimpleMeans_Table(:,5),variableNames);
        emotionIntensity_Table_S = tableconcat_start(window,windowSimpleMeans_Table(:,6),variableNames);
        singleEmotion_Table_S = tableconcat_start(window,windowSimpleMeans_Table(:,7),variableNames);
        extremeScaleUse_Table_S = tableconcat_start(window,windowSimpleMeans_Table(:,8),variableNames);
        mPositive_Table_S = tableconcat_start(window,windowSimpleMeans_Table(:,9),variableNames);
        mNegative_Table_S = tableconcat_start(window,windowSimpleMeans_Table(:,10),variableNames);
        sdPositive_Table_S = tableconcat_start(window,windowSimpleMeans_Table(:,11),variableNames);
        sdNegative_Table_S = tableconcat_start(window,windowSimpleMeans_Table(:,12),variableNames);
        
        zICC_Neg_Table_W = tableconcat_start(window,windowWeightedMeans_Table(:,1),variableNames);
        zICC_Pos_Table_W = tableconcat_start(window,windowWeightedMeans_Table(:,2),variableNames);
        zICC_Overall_Table_W = tableconcat_start(window,windowWeightedMeans_Table(:,3),variableNames);
        emodiversity_Table_W = tableconcat_start(window,windowWeightedMeans_Table(:,4),variableNames);
        instability_Table_W = tableconcat_start(window,windowWeightedMeans_Table(:,5),variableNames);
        emotionIntensity_Table_W = tableconcat_start(window,windowWeightedMeans_Table(:,6),variableNames);
        singleEmotion_Table_W = tableconcat_start(window,windowWeightedMeans_Table(:,7),variableNames);
        extremeScaleUse_Table_W = tableconcat_start(window,windowWeightedMeans_Table(:,8),variableNames);
        mPositive_Table_W = tableconcat_start(window,windowWeightedMeans_Table(:,9),variableNames);
        mNegative_Table_W = tableconcat_start(window,windowWeightedMeans_Table(:,10),variableNames);
        sdPositive_Table_W = tableconcat_start(window,windowWeightedMeans_Table(:,11),variableNames);
        sdNegative_Table_W = tableconcat_start(window,windowWeightedMeans_Table(:,12),variableNames);
    else
        zICC_Neg_Table_S = tableconcat_iter(zICC_Neg_Table_S,windowSimpleMeans_Table(:,1),i_subject,ppID);
        zICC_Pos_Table_S = tableconcat_iter(zICC_Pos_Table_S,windowSimpleMeans_Table(:,2),i_subject,ppID);
        zICC_Overall_Table_S = tableconcat_iter(zICC_Overall_Table_S,windowSimpleMeans_Table(:,3),i_subject,ppID);
        emodiversity_Table_S = tableconcat_iter(emodiversity_Table_S,windowSimpleMeans_Table(:,4),i_subject,ppID);
        instability_Table_S = tableconcat_iter(instability_Table_S,windowSimpleMeans_Table(:,5),i_subject,ppID);
        emotionIntensity_Table_S = tableconcat_iter(emotionIntensity_Table_S,windowSimpleMeans_Table(:,6),i_subject,ppID);
        singleEmotion_Table_S = tableconcat_iter(singleEmotion_Table_S,windowSimpleMeans_Table(:,7),i_subject,ppID);
        extremeScaleUse_Table_S = tableconcat_iter(extremeScaleUse_Table_S,windowSimpleMeans_Table(:,8),i_subject,ppID);
        mPositive_Table_S = tableconcat_iter(mPositive_Table_S,windowSimpleMeans_Table(:,9),i_subject,ppID);
        mNegative_Table_S = tableconcat_iter(mNegative_Table_S,windowSimpleMeans_Table(:,10),i_subject,ppID);
        sdPositive_Table_S = tableconcat_iter(sdPositive_Table_S,windowSimpleMeans_Table(:,11),i_subject,ppID);
        sdNegative_Table_S = tableconcat_iter(sdNegative_Table_S,windowSimpleMeans_Table(:,12),i_subject,ppID);
        
        zICC_Neg_Table_W = tableconcat_iter(zICC_Neg_Table_W,windowWeightedMeans_Table(:,1),i_subject,ppID);
        zICC_Pos_Table_W = tableconcat_iter(zICC_Pos_Table_W,windowWeightedMeans_Table(:,2),i_subject,ppID);
        zICC_Overall_Table_W = tableconcat_iter(zICC_Overall_Table_W,windowWeightedMeans_Table(:,3),i_subject,ppID);
        emodiversity_Table_W = tableconcat_iter(emodiversity_Table_W,windowWeightedMeans_Table(:,4),i_subject,ppID);
        instability_Table_W = tableconcat_iter(instability_Table_W,windowWeightedMeans_Table(:,5),i_subject,ppID);
        emotionIntensity_Table_W = tableconcat_iter(emotionIntensity_Table_W,windowWeightedMeans_Table(:,6),i_subject,ppID);
        singleEmotion_Table_W = tableconcat_iter(singleEmotion_Table_W,windowWeightedMeans_Table(:,7),i_subject,ppID);
        extremeScaleUse_Table_W = tableconcat_iter(extremeScaleUse_Table_W,windowWeightedMeans_Table(:,8),i_subject,ppID);
        mPositive_Table_W = tableconcat_iter(mPositive_Table_W,windowWeightedMeans_Table(:,9),i_subject,ppID);
        mNegative_Table_W = tableconcat_iter(mNegative_Table_W,windowWeightedMeans_Table(:,10),i_subject,ppID);
        sdPositive_Table_W = tableconcat_iter(sdPositive_Table_W,windowWeightedMeans_Table(:,11),i_subject,ppID);
        sdNegative_Table_W = tableconcat_iter(sdNegative_Table_W,windowWeightedMeans_Table(:,12),i_subject,ppID);
    end
    clearvars rawICC rtoZ emodiversity instability emotionIntensity singleEmotion extremeScaleUse mPositive mNegative sdPositive sdNegative
    clearvars measures windowSimpleMeans windowWeightedMeans dayIDlist
end

%% create summary tables
save(['dataSet' num2str(dataSet) '_measureSimpleMeans_windows.mat'],'measureSimpleMeans');
save(['dataSet' num2str(dataSet) '_measureSimpleSDs_windows.mat'],'measureSimpleSDs');
save(['dataSet' num2str(dataSet) '_measureWeightedMeans_windows.mat'],'measureWeightedMeans');
save(['dataSet' num2str(dataSet) '_measureWeightedSDs_windows.mat'],'measureWeightedSDs');
variableNames = {'zICC_N','zICC_P','zICC_M','emodiv','instab','intens','singleEmo','extreme','mPos','mNeg','sdPos','sdNeg'};
measureSimpleMeans_Table = array2table(measureSimpleMeans,'VariableNames',variableNames);
measureSimpleSDs_Table = array2table(measureSimpleSDs,'VariableNames',variableNames);
measureWeightedMeans_Table = array2table(measureWeightedMeans,'VariableNames',variableNames);
measureWeightedSDs_Table = array2table(measureWeightedSDs,'VariableNames',variableNames);
if printSummaries == 1
    writetable(measureSimpleMeans_Table,['dataSet' num2str(dataSet) '_measureSimpleMeans_windows.xlsx']);
    writetable(measureSimpleSDs_Table,['dataSet' num2str(dataSet) '_measureSimpleSDs_windows.xlsx']);
    writetable(measureWeightedMeans_Table,['dataSet' num2str(dataSet) '_measureWeightedMeans_windows.xlsx']);
    writetable(measureWeightedSDs_Table,['dataSet' num2str(dataSet) '_measureWeightedSDs_windows.xlsx']);
end

%% write tables to file
if printTables == 1 
    writetable(zICC_Neg_Table_S,['dataSet' num2str(dataSet) '_zICC_Neg_windows_simple.xlsx']);
    writetable(zICC_Pos_Table_S,['dataSet' num2str(dataSet) '_zICC_Pos_windows_simple.xlsx']);
    writetable(zICC_Overall_Table_S,['dataSet' num2str(dataSet) '_zICC_Overall_windows_simple.xlsx']);
    writetable(emodiversity_Table_S,['dataSet' num2str(dataSet) '_emodiversity_windows_simple.xlsx']);
    writetable(instability_Table_S,['dataSet' num2str(dataSet) '_instability_windows_simple.xlsx']);
    writetable(emotionIntensity_Table_S,['dataSet' num2str(dataSet) '_emotionIntensity_windows_simple.xlsx']);
    writetable(singleEmotion_Table_S,['dataSet' num2str(dataSet) '_singleEmotion_windows_simple.xlsx']);
    writetable(extremeScaleUse_Table_S,['dataSet' num2str(dataSet) '_extremeScaleUse_windows_simple.xlsx']);
    writetable(mPositive_Table_S,['dataSet' num2str(dataSet) '_mPositive_windows_simple.xlsx']);
    writetable(mNegative_Table_S,['dataSet' num2str(dataSet) '_mNegative_windows_simple.xlsx']);
    writetable(sdPositive_Table_S,['dataSet' num2str(dataSet) '_sdPositive_windows_simple.xlsx']);
    writetable(sdNegative_Table_S,['dataSet' num2str(dataSet) '_sdNegative_windows_simple.xlsx']);  
    
    writetable(zICC_Neg_Table_W,['dataSet' num2str(dataSet) '_zICC_Neg_windows_weighted.xlsx']);
    writetable(zICC_Pos_Table_W,['dataSet' num2str(dataSet) '_zICC_Pos_windows_weighted.xlsx']);
    writetable(zICC_Overall_Table_W,['dataSet' num2str(dataSet) '_zICC_Overall_windows_weighted.xlsx']);
    writetable(emodiversity_Table_W,['dataSet' num2str(dataSet) '_emodiversity_windows_weighted.xlsx']);
    writetable(instability_Table_W,['dataSet' num2str(dataSet) '_instability_windows_weighted.xlsx']);
    writetable(emotionIntensity_Table_W,['dataSet' num2str(dataSet) '_emotionIntensity_windows_weighted.xlsx']);
    writetable(singleEmotion_Table_W,['dataSet' num2str(dataSet) '_singleEmotion_windows_weighted.xlsx']);
    writetable(extremeScaleUse_Table_W,['dataSet' num2str(dataSet) '_extremeScaleUse_windows_weighted.xlsx']);
    writetable(mPositive_Table_W,['dataSet' num2str(dataSet) '_mPositive_windows_weighted.xlsx']);
    writetable(mNegative_Table_W,['dataSet' num2str(dataSet) '_mNegative_windows_weighted.xlsx']);
    writetable(sdPositive_Table_W,['dataSet' num2str(dataSet) '_sdPositive_windows_weighted.xlsx']);
    writetable(sdNegative_Table_W,['dataSet' num2str(dataSet) '_sdNegative_windows_weighted.xlsx']);  
end

%% run Kolmogorov-Smirnov test to assess normality; normalize distributions as needed
[~,numVariables] = size(measureSimpleMeans);
for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(measureSimpleMeans(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(measureSimpleMeans(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm_SimpleMeans(:,i_variable) = norminv(fracRank,mean(measureSimpleMeans(:,i_variable)),std(measureSimpleMeans(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm_SimpleMeans(:,i_variable) = measureSimpleMeans(:,i_variable);
    end
end

for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(measureSimpleSDs(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(measureSimpleSDs(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm_SimpleSDs(:,i_variable) = norminv(fracRank,mean(measureSimpleSDs(:,i_variable)),std(measureSimpleSDs(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm_SimpleSDs(:,i_variable) = measureSimpleSDs(:,i_variable);
    end
end

for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(measureWeightedMeans(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(measureWeightedMeans(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm_WeightedMeans(:,i_variable) = norminv(fracRank,mean(measureWeightedMeans(:,i_variable)),std(measureWeightedMeans(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm_WeightedMeans(:,i_variable) = measureWeightedMeans(:,i_variable);
    end
end

for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(measureWeightedSDs(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(measureWeightedSDs(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm_WeightedSDs(:,i_variable) = norminv(fracRank,mean(measureWeightedSDs(:,i_variable)),std(measureWeightedSDs(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm_WeightedSDs(:,i_variable) = measureWeightedSDs(:,i_variable);
    end
end

%% compile all (normalized) data and re-export
save(['dataSet' num2str(dataSet) '_measureSimpleMeans_normalized_windows.mat'],'fracRankNorm_SimpleMeans');
save(['dataSet' num2str(dataSet) '_measureSimpleSDs_normalized_windows.mat'],'fracRankNorm_SimpleSDs');
save(['dataSet' num2str(dataSet) '_measureWeightedMeans_normalized_windows.mat'],'fracRankNorm_WeightedMeans');
save(['dataSet' num2str(dataSet) '_measureWeightedSDs_normalized_windows.mat'],'fracRankNorm_WeightedSDs');
fracRankNorm_SimpleMeans_Table = array2table(fracRankNorm_SimpleMeans,'VariableNames',variableNames);
fracRankNorm_SimpleSDs_Table = array2table(fracRankNorm_SimpleSDs,'VariableNames',variableNames);
fracRankNorm_WeightedMeans_Table = array2table(fracRankNorm_WeightedMeans,'VariableNames',variableNames);
fracRankNorm_WeightedSDs_Table = array2table(fracRankNorm_WeightedSDs,'VariableNames',variableNames);
if printSummaries == 1
    writetable(fracRankNorm_SimpleMeans_Table,['dataSet' num2str(dataSet) '_measureSimpleMeans_normalized_windows.xlsx']);
    writetable(fracRankNorm_SimpleSDs_Table,['dataSet' num2str(dataSet) '_measureSimpleSDs_normalized_windows.xlsx']);
    writetable(fracRankNorm_WeightedMeans_Table,['dataSet' num2str(dataSet) '_measureWeightedMeans_normalized_windows.xlsx']);
    writetable(fracRankNorm_WeightedSDs_Table,['dataSet' num2str(dataSet) '_measureWeightedSDs_normalized_windows.xlsx']);
end

%% save tables as matrices
if saveData == 1
    subjectCol = [0; subjectIDlist];

    zICC_Neg_Array_S = savedata(zICC_Neg_Table_S,subjectCol,'zICC_Neg_windows_simple',dataSet);
    zICC_Pos_Array_S = savedata(zICC_Pos_Table_S,subjectCol,'zICC_Pos_windows_simple',dataSet);
    zICC_Overall_Array_S = savedata(zICC_Overall_Table_S,subjectCol,'zICC_Overall_windows_simple',dataSet);
    emodiversity_Array_S = savedata(emodiversity_Table_S,subjectCol,'emodiversity_windows_simple',dataSet);
    instability_Array_S = savedata(instability_Table_S,subjectCol,'instability_windows_simple',dataSet);
    emotionIntensity_Array_S = savedata(emotionIntensity_Table_S,subjectCol,'emotionIntensity_windows_simple',dataSet);
    singleEmotion_Array_S = savedata(singleEmotion_Table_S,subjectCol,'singleEmotion_windows_simple',dataSet);
    extremeScaleUse_Array_S = savedata(extremeScaleUse_Table_S,subjectCol,'extremeScaleUse_windows_simple',dataSet);
    mPositive_Array_S = savedata(mPositive_Table_S,subjectCol,'mPositive_windows_simple',dataSet);
    mNegative_Array_S = savedata(mNegative_Table_S,subjectCol,'mNegative_windows_simple',dataSet);
    sdPositive_Array_S = savedata(sdPositive_Table_S,subjectCol,'sdPositive_windows_simple',dataSet);
    sdNegative_Array_S = savedata(sdNegative_Table_S,subjectCol,'sdNegative_windows_simple',dataSet);

    zICC_Neg_Array_W = savedata(zICC_Neg_Table_W,subjectCol,'zICC_Neg_windows_weighted',dataSet);
    zICC_Pos_Array_W = savedata(zICC_Pos_Table_W,subjectCol,'zICC_Pos_windows_weighted',dataSet);
    zICC_Overall_Array_W = savedata(zICC_Overall_Table_W,subjectCol,'zICC_Overall_windows_weighted',dataSet);
    emodiversity_Array_W = savedata(emodiversity_Table_W,subjectCol,'emodiversity_windows_weighted',dataSet);
    instability_Array_W = savedata(instability_Table_W,subjectCol,'instability_windows_weighted',dataSet);
    emotionIntensity_Array_W = savedata(emotionIntensity_Table_W,subjectCol,'emotionIntensity_windows_weighted',dataSet);
    singleEmotion_Array_W = savedata(singleEmotion_Table_W,subjectCol,'singleEmotion_windows_weighted',dataSet);
    extremeScaleUse_Array_W = savedata(extremeScaleUse_Table_W,subjectCol,'extremeScaleUse_windows_weighted',dataSet);
    mPositive_Array_W = savedata(mPositive_Table_W,subjectCol,'mPositive_windows_weighted',dataSet);
    mNegative_Array_W = savedata(mNegative_Table_W,subjectCol,'mNegative_windows_weighted',dataSet);
    sdPositive_Array_W = savedata(sdPositive_Table_W,subjectCol,'sdPositive_windows_weighted',dataSet);
    sdNegative_Array_W = savedata(sdNegative_Table_W,subjectCol,'sdNegative_windows_weighted',dataSet);
end