clear;
clc;

%% specify dataset and parameters
dataSet = 18; % set to 18, 39, or 88
print = 1;

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
    % remove invalid values and missing data
    subjectData(subjectData > maximumRating) = NaN;
    missingData = isnan(subjectData); 
    missingData2 = any(missingData,2); 
    subjectData = subjectData(~missingData2,:); 
    % if necessary, rescale data to start at 0
    if startRatingsat1 == 1
        subjectData = subjectData-1;
    end
    % calculate granularity (ICCs and zICCs)
    rawICC(i_subject,1) = ICC(subjectData(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
    rawICC(i_subject,2) = ICC(subjectData(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
    rawICC(i_subject,3) = (rawICC(i_subject,1)+rawICC(i_subject,2))/2; % valence mean ICC
    rawICC(rawICC<0) = 0;
    rtoZ(i_subject,1) = 0.5*log((1+rawICC(i_subject,1))/(1-rawICC(i_subject,1)));
    rtoZ(i_subject,2) = 0.5*log((1+rawICC(i_subject,2))/(1-rawICC(i_subject,2)));
    rtoZ(i_subject,3) = 0.5*log((1+rawICC(i_subject,3))/(1-rawICC(i_subject,3)));
    % calculate emodiversity
    totalCount = sum(sum(subjectData>0));
    for i_emotion = 1:length(wordList)
        countEmotion = sum(subjectData(:,i_emotion)>0);
        propEmotion(i_emotion) = (countEmotion./totalCount)*log(countEmotion./totalCount);
    end
    propEmotion = propEmotion(isfinite(propEmotion));
    emodiversity(i_subject) = -1*sum(propEmotion);
    % calculate instability
    for i_emotion = 1:length(wordList)
        diffEmotion = diff(subjectData(:,i_emotion));
        MSSD_emotion(i_emotion) = mean(diffEmotion.^2);
    end
    instability(i_subject) = mean(MSSD_emotion);
    % calculate average emotional intensity
    emotionIntensity(i_subject) = mean(mean(subjectData));
    % calculate single emotion experience
    [prompts,~] = size(subjectData);
    for i_prompt = 1:prompts
        highestEmotion = max(subjectData(i_prompt,:));
        distances = highestEmotion-subjectData(i_prompt,:);
        distancesClean = distances(distances ~= 0);
        meanDistance(i_prompt) = mean(distancesClean);
    end
    meanDistance(isnan(meanDistance)) = 0; % for instances where all emotions are rated as 0, meanDistance=0    
    singleEmotion(i_subject) = mean(meanDistance)/scaleMax; % calculate as ratio to max distance
    % calculate extreme scale use
    midpointDistance = abs(subjectData-scaleMidpoint);
    extremeScaleUse(i_subject) = mean(mean(midpointDistance))/scaleMidpoint; % calculate as ratio to possible distance from scale midpoint
    % calculate mean and SD for positive and negative valence
    mPositive(i_subject) = mean(mean(subjectData(:,(valence==1))));
    mNegative(i_subject) = mean(mean(subjectData(:,(valence==2))));
    sdPositive(i_subject) = mean(std(subjectData(:,(valence==1))));
    sdNegative(i_subject) = mean(std(subjectData(:,(valence==2)))); 
end

%% create summary table
measures = horzcat(rtoZ, emodiversity', instability', emotionIntensity', singleEmotion', extremeScaleUse', mPositive', mNegative', sdPositive', sdNegative');
save(['dataSet' num2str(dataSet) '_measures.mat'],'measures');
variableNames = {'zICC_N','zICC_P','zICC_M','emodiv','instab','intens','singleEmo','extreme','mPos','mNeg','sdPos','sdNeg'};
measures_Table = array2table(measures,'VariableNames',variableNames);
if print == 1
    writetable(measures_Table,['dataSet' num2str(dataSet) '_measures.xlsx']);
end

%% run Kolmogorov-Smirnov test to assess normality; normalize distributions as needed
[~,numVariables] = size(measures);
for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(measures(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(measures(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm(:,i_variable) = norminv(fracRank,mean(measures(:,i_variable)),std(measures(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm(:,i_variable) = measures(:,i_variable);
    end
end

%% compile all (normalized) data and re-export
save(['dataSet' num2str(dataSet) '_measures_normalized.mat'],'fracRankNorm');
fracRankNorm_Table = array2table(fracRankNorm,'VariableNames',variableNames);
if print == 1
    writetable(fracRankNorm_Table,['dataSet' num2str(dataSet) '_measures_normalized.xlsx']);
end