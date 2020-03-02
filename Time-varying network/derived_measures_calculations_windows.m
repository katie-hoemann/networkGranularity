% NOTE: Currently built out for 3-day sliding window; modification to
% window size or weighting will require structural edits to script

clear;
clc;

%% specify dataset and parameters
dataSet = 18; % only 18 or 39
%windowSize = 3; % size in days
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
        % if necessary, rescale data to start at 0
        if startRatingsat1 == 1
            windowData = windowData-1;
        end
        % calculate granularity (ICCs and zICCs)
        rawICC(i_window,1) = ICC(windowData(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
        rawICC(i_window,2) = ICC(windowData(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
        rawICC(i_window,3) = (rawICC(i_window,1)+rawICC(i_window,2))/2; % valence mean ICC
        rawICC(rawICC<0) = 0;
        rtoZ(i_window,1) = 0.5*log((1+rawICC(i_window,1))/(1-rawICC(i_window,1)));
        rtoZ(i_window,2) = 0.5*log((1+rawICC(i_window,2))/(1-rawICC(i_window,2)));
        rtoZ(i_window,3) = 0.5*log((1+rawICC(i_window,3))/(1-rawICC(i_window,3)));
        % calculate emodiversity
        totalCount = sum(sum(windowData>0));
        for i_emotion = 1:length(wordList)
            countEmotion = sum(windowData(:,i_emotion)>0);
            propEmotion(i_emotion) = (countEmotion./totalCount)*log(countEmotion./totalCount);
        end
        propEmotion = propEmotion(isfinite(propEmotion));
        emodiversity(i_window) = -1*sum(propEmotion);
        % calculate instability
        for i_emotion = 1:length(wordList)
            diffEmotion = diff(windowData(:,i_emotion));
            MSSD_emotion(i_emotion) = mean(diffEmotion.^2);
        end
        instability(i_window) = mean(MSSD_emotion);
        % calculate average emotional intensity
        emotionIntensity(i_window) = mean(mean(windowData));
        % calculate single emotion experience
        [prompts,~] = size(windowData);
        for i_prompt = 1:prompts
            highestEmotion = max(windowData(i_prompt,:));
            distances = highestEmotion-windowData(i_prompt,:);
            distancesClean = distances(distances ~= 0);
            meanDistance(i_prompt) = mean(distancesClean);
        end
        meanDistance(isnan(meanDistance)) = 0; % for instances where all emotions are rated as 0, meanDistance=0    
        singleEmotion(i_window) = mean(meanDistance)/scaleMax; % calculate as ratio to max distance
        % calculate extreme scale use
        midpointDistance = abs(windowData-scaleMidpoint);
        extremeScaleUse(i_window) = mean(mean(midpointDistance))/scaleMidpoint; % calculate as ratio to possible distance from scale midpoint
        % calculate mean and SD for positive and negative valence
        mPositive(i_window) = mean(mean(windowData(:,(valence==1))));
        mNegative(i_window) = mean(mean(windowData(:,(valence==2))));
        sdPositive(i_window) = mean(std(windowData(:,(valence==1))));
        sdNegative(i_window) = mean(std(windowData(:,(valence==2))));
    end
    %% calculate subject-level mean and SD for measures
    measureMeans(i_subject,1) = mean(rtoZ(:,1),'omitnan');
    measureMeans(i_subject,2) = mean(rtoZ(:,2),'omitnan');
    measureMeans(i_subject,3) = mean(rtoZ(:,3),'omitnan');
    measureMeans(i_subject,4) = mean(emodiversity,'omitnan');
    measureMeans(i_subject,5) = mean(instability,'omitnan');
    measureMeans(i_subject,6) = mean(emotionIntensity,'omitnan');
    measureMeans(i_subject,7) = mean(singleEmotion,'omitnan');
    measureMeans(i_subject,8) = mean(extremeScaleUse,'omitnan');
    measureMeans(i_subject,9) = mean(mPositive,'omitnan');
    measureMeans(i_subject,10) = mean(mNegative,'omitnan');
    measureMeans(i_subject,11) = mean(sdPositive,'omitnan');
    measureMeans(i_subject,12) = mean(sdNegative,'omitnan');
    
    measureSDs(i_subject,1) = std(rtoZ(:,1),'omitnan');
    measureSDs(i_subject,2) = std(rtoZ(:,2),'omitnan');
    measureSDs(i_subject,3) = std(rtoZ(:,3),'omitnan');
    measureSDs(i_subject,4) = std(emodiversity,'omitnan');
    measureSDs(i_subject,5) = std(instability,'omitnan');
    measureSDs(i_subject,6) = std(emotionIntensity,'omitnan');
    measureSDs(i_subject,7) = std(singleEmotion,'omitnan');
    measureSDs(i_subject,8) = std(extremeScaleUse,'omitnan');
    measureSDs(i_subject,9) = std(mPositive,'omitnan');
    measureSDs(i_subject,10) = std(mNegative,'omitnan');
    measureSDs(i_subject,11) = std(sdPositive,'omitnan');
    measureSDs(i_subject,12) = std(sdNegative,'omitnan');    
    %% add subject results to summary table
    ppID = ['PP' num2str(subjectID)];
    windowNumbers = dayIDlist(1:end-2);
    if i_subject == 1
        window = (1:1:max(allData(:,2))-2)';
        window = array2table(window,'VariableNames',{'Window'});
     
        zICC_Neg_Table = tablecompile_start(window,windowNumbers,rtoZ(:,1),ppID,'Window','ppWindow','ppWindow');
        zICC_Pos_Table = tablecompile_start(window,windowNumbers,rtoZ(:,2),ppID,'Window','ppWindow','ppWindow');
        zICC_Overall_Table = tablecompile_start(window,windowNumbers,rtoZ(:,3),ppID,'Window','ppWindow','ppWindow');
        emodiversity_Table = tablecompile_start(window,windowNumbers,emodiversity(:),ppID,'Window','ppWindow','ppWindow');
        instability_Table = tablecompile_start(window,windowNumbers,instability(:),ppID,'Window','ppWindow','ppWindow');
        emotionIntensity_Table = tablecompile_start(window,windowNumbers,emotionIntensity(:),ppID,'Window','ppWindow','ppWindow');
        singleEmotion_Table = tablecompile_start(window,windowNumbers,singleEmotion(:),ppID,'Window','ppWindow','ppWindow');
        extremeScaleUse_Table = tablecompile_start(window,windowNumbers,extremeScaleUse(:),ppID,'Window','ppWindow','ppWindow');
        mPositive_Table = tablecompile_start(window,windowNumbers,mPositive(:),ppID,'Window','ppWindow','ppWindow');
        mNegative_Table = tablecompile_start(window,windowNumbers,mNegative(:),ppID,'Window','ppWindow','ppWindow');
        sdPositive_Table = tablecompile_start(window,windowNumbers,sdPositive(:),ppID,'Window','ppWindow','ppWindow');
        sdNegative_Table = tablecompile_start(window,windowNumbers,sdNegative(:),ppID,'Window','ppWindow','ppWindow');
    else
        zICC_Neg_Table = tablecompile_iter(zICC_Neg_Table,windowNumbers,rtoZ(:,1),ppID,'Window','ppWindow','ppWindow');
        zICC_Pos_Table = tablecompile_iter(zICC_Pos_Table,windowNumbers,rtoZ(:,2),ppID,'Window','ppWindow','ppWindow');
        zICC_Overall_Table = tablecompile_iter(zICC_Overall_Table,windowNumbers,rtoZ(:,3),ppID,'Window','ppWindow','ppWindow');
        emodiversity_Table = tablecompile_iter(emodiversity_Table,windowNumbers,emodiversity(:),ppID,'Window','ppWindow','ppWindow');
        instability_Table = tablecompile_iter(instability_Table,windowNumbers,instability(:),ppID,'Window','ppWindow','ppWindow');
        emotionIntensity_Table = tablecompile_iter(emotionIntensity_Table,windowNumbers,emotionIntensity(:),ppID,'Window','ppWindow','ppWindow');
        singleEmotion_Table = tablecompile_iter(singleEmotion_Table,windowNumbers,singleEmotion(:),ppID,'Window','ppWindow','ppWindow');
        extremeScaleUse_Table = tablecompile_iter(extremeScaleUse_Table,windowNumbers,extremeScaleUse(:),ppID,'Window','ppWindow','ppWindow');
        mPositive_Table = tablecompile_iter(mPositive_Table,windowNumbers,mPositive(:),ppID,'Window','ppWindow','ppWindow');
        mNegative_Table = tablecompile_iter(mNegative_Table,windowNumbers,mNegative(:),ppID,'Window','ppWindow','ppWindow');
        sdPositive_Table = tablecompile_iter(sdPositive_Table,windowNumbers,sdPositive(:),ppID,'Window','ppWindow','ppWindow');
        sdNegative_Table = tablecompile_iter(sdNegative_Table,windowNumbers,sdNegative(:),ppID,'Window','ppWindow','ppWindow');
    end
    clear rawICC rtoZ emodiversity instability emotionIntensity singleEmotion extremeScaleUse mPositive mNegative sdPositive sdNegative
end

%% create summary tables
save(['dataSet' num2str(dataSet) '_measureMeans_windows.mat'],'measureMeans');
save(['dataSet' num2str(dataSet) '_measureSDs_windows.mat'],'measureSDs');
variableNames = {'zICC_N','zICC_P','zICC_M','emodiv','instab','intens','singleEmo','extreme','mPos','mNeg','sdPos','sdNeg'};
measureMeans_Table = array2table(measureMeans,'VariableNames',variableNames);
measureSDs_Table = array2table(measureSDs,'VariableNames',variableNames);
if printSummaries == 1
    writetable(measureMeans_Table,['dataSet' num2str(dataSet) '_measureMeans_windows.xlsx']);
    writetable(measureSDs_Table,['dataSet' num2str(dataSet) '_measureSDs_windows.xlsx']);
end

%% write tables to file
if printTables == 1 
    writetable(zICC_Neg_Table,['dataSet' num2str(dataSet) '_zICC_Neg_windows.xlsx']);
    writetable(zICC_Pos_Table,['dataSet' num2str(dataSet) '_zICC_Pos_windows.xlsx']);
    writetable(zICC_Overall_Table,['dataSet' num2str(dataSet) '_zICC_Overall_windows.xlsx']);
    writetable(emodiversity_Table,['dataSet' num2str(dataSet) '_emodiversity_windows.xlsx']);
    writetable(instability_Table,['dataSet' num2str(dataSet) '_instability_windows.xlsx']);
    writetable(emotionIntensity_Table,['dataSet' num2str(dataSet) '_emotionIntensity_windows.xlsx']);
    writetable(singleEmotion_Table,['dataSet' num2str(dataSet) '_singleEmotion_windows.xlsx']);
    writetable(extremeScaleUse_Table,['dataSet' num2str(dataSet) '_extremeScaleUse_windows.xlsx']);
    writetable(mPositive_Table,['dataSet' num2str(dataSet) '_mPositive_windows.xlsx']);
    writetable(mNegative_Table,['dataSet' num2str(dataSet) '_mNegative_windows.xlsx']);
    writetable(sdPositive_Table,['dataSet' num2str(dataSet) '_sdPositive_windows.xlsx']);
    writetable(sdNegative_Table,['dataSet' num2str(dataSet) '_sdNegative_windows.xlsx']);  
end

%% run Kolmogorov-Smirnov test to assess normality; normalize distributions as needed
[~,numVariables] = size(measureMeans);
for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(measureMeans(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(measureMeans(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm_Means(:,i_variable) = norminv(fracRank,mean(measureMeans(:,i_variable)),std(measureMeans(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm_Means(:,i_variable) = measureMeans(:,i_variable);
    end
end

for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(measureSDs(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(measureSDs(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm_SDs(:,i_variable) = norminv(fracRank,mean(measureSDs(:,i_variable)),std(measureSDs(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm_SDs(:,i_variable) = measureSDs(:,i_variable);
    end
end

%% compile all (normalized) data and re-export
save(['dataSet' num2str(dataSet) '_measureMeans_normalized_windows.mat'],'fracRankNorm_Means');
save(['dataSet' num2str(dataSet) '_measureSDs_normalized_windows.mat'],'fracRankNorm_SDs');
fracRankNorm_Means_Table = array2table(fracRankNorm_Means,'VariableNames',variableNames);
fracRankNorm_SDs_Table = array2table(fracRankNorm_SDs,'VariableNames',variableNames);
if printSummaries == 1
    writetable(fracRankNorm_Means_Table,['dataSet' num2str(dataSet) '_measureMeans_normalized_windows.xlsx']);
    writetable(fracRankNorm_SDs_Table,['dataSet' num2str(dataSet) '_measuresSDs_normalized_windows.xlsx']);
end

%% save tables as matrices
if saveData == 1
    subjectCol = [0; subjectIDlist];

    zICC_Neg_Array = table2array(zICC_Neg_Table)';
    zICC_Neg_Array = horzcat(subjectCol,zICC_Neg_Array);
    save(['dataSet' num2str(dataSet) '_zICC_Neg_windows.mat'],'zICC_Neg_Array');

    zICC_Pos_Array = table2array(zICC_Pos_Table)';
    zICC_Pos_Array = horzcat(subjectCol,zICC_Pos_Array);
    save(['dataSet' num2str(dataSet) '_zICC_Pos_windows.mat'],'zICC_Pos_Array');

    zICC_Overall_Array = table2array(zICC_Overall_Table)';
    zICC_Overall_Array = horzcat(subjectCol,zICC_Overall_Array);
    save(['dataSet' num2str(dataSet) '_zICC_Overall_windows.mat'],'zICC_Overall_Array');

    emodiversity_Array = table2array(emodiversity_Table)';
    emodiversity_Array = horzcat(subjectCol,emodiversity_Array);
    save(['dataSet' num2str(dataSet) '_emodiversity_windows.mat'],'emodiversity_Array');

    instability_Array = table2array(instability_Table)';
    instability_Array = horzcat(subjectCol,instability_Array);
    save(['dataSet' num2str(dataSet) '_instability_windows.mat'],'instability_Array');

    emotionIntensity_Array = table2array(emotionIntensity_Table)';
    emotionIntensity_Array = horzcat(subjectCol,emotionIntensity_Array);
    save(['dataSet' num2str(dataSet) '_emotionIntensity_windows.mat'],'emotionIntensity_Array');

    singleEmotion_Array = table2array(singleEmotion_Table)';
    singleEmotion_Array = horzcat(subjectCol,singleEmotion_Array);
    save(['dataSet' num2str(dataSet) '_singleEmotion_windows.mat'],'singleEmotion_Array');

    extremeScaleUse_Array = table2array(extremeScaleUse_Table)';
    extremeScaleUse_Array = horzcat(subjectCol,extremeScaleUse_Array);
    save(['dataSet' num2str(dataSet) '_extremeScaleUse_windows.mat'],'extremeScaleUse_Array');

    mPositive_Array = table2array(mPositive_Table)';
    mPositive_Array = horzcat(subjectCol,mPositive_Array);
    save(['dataSet' num2str(dataSet) '_mPositive_windows.mat'],'mPositive_Array');

    mNegative_Array = table2array(mNegative_Table)';
    mNegative_Array = horzcat(subjectCol,mNegative_Array);
    save(['dataSet' num2str(dataSet) '_mNegative_windows.mat'],'mNegative_Array');

    sdPositive_Array = table2array(sdPositive_Table)';
    sdPositive_Array = horzcat(subjectCol,sdPositive_Array);
    save(['dataSet' num2str(dataSet) '_sdPositive_windows.mat'],'sdPositive_Array');

    sdNegative_Array = table2array(sdNegative_Table)';
    sdNegative_Array = horzcat(subjectCol,sdNegative_Array);
    save(['dataSet' num2str(dataSet) '_sdNegative_windows.mat'],'sdNegative_Array');
end