% This script calculates emodiversity (based on both Shannon's formula and the Gini coefficient) for all subjects in the specified data set

clear;
clc;

%% specify dataset
dataSet = 18; % set to 18, 39, or 88 and the rest will take care of itself!

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

%% grab data for each subject and run through calculation
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
    %% calculate emodiversity using Shannon formula (Quoidbach et al, 2014)
    totalCount = sum(sum(subjectData>0));
    for i_emotion = 1:length(wordList)
        countEmotion = sum(subjectData(:,i_emotion)>0);
        propEmotion(i_emotion) = (countEmotion./totalCount)*log(countEmotion./totalCount);
    end
    propEmotion = propEmotion(isfinite(propEmotion)); % remove values that are not finite (i.e., emotions never experienced)
    emodiversity(i_subject) = -1*sum(propEmotion);
    %% calculate emodiversity using Gini coefficient (Benson et al, 2018)
    numEmotions = length(wordList);
    for i_emotion = 1:numEmotions
        countEmotionRanked(i_emotion) = sum(subjectData(:,i_emotion)>0);
    end
    countEmotionRanked = sort(countEmotionRanked);
    index = 1:1:numEmotions;
    for i_emotion = 1:numEmotions
        weightedCount(i_emotion) = countEmotionRanked(i_emotion)*index(i_emotion);
    end
    giniCoef(i_subject) = 1-(((2*sum(weightedCount))/(numEmotions*sum(countEmotionRanked)))-((numEmotions+1)/numEmotions));
end

%% save off matrices for construct validity checks
emodiversity = emodiversity';
filename = ['dataSet' num2str(dataSet) '_emodiversity.mat'];
save(filename,'emodiversity');

giniCoef = giniCoef';
filename = ['dataSet' num2str(dataSet) '_giniCoef.mat'];
save(filename,'giniCoef');