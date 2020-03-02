% This script calculates single emotion experience, average emotion intensity, and extreme scale use and for all subjects in the specified data set

clear;
clc;

%% specify dataset
dataSet = 88; % set to 18, 39, or 88 and the rest will take care of itself!

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
    %% calculate single emotion experience
    [prompts,~] = size(subjectData);
    for i_prompt = 1:prompts
        highestEmotion = max(subjectData(i_prompt,:));
        distances = highestEmotion-subjectData(i_prompt,:);
        distancesClean = distances(distances ~= 0);
        meanDistance(i_prompt) = mean(distancesClean);
    end
%     meanDistanceValues = ~isnan(meanDistance);
%     meanDistanceClean = meanDistance(meanDistanceValues);
    meanDistance(isnan(meanDistance)) = 0; % for instances where all emotions are rated as 0, meanDistance=0    
    singleEmotion(i_subject) = mean(meanDistance)/scaleMax; % calculate as ratio to max distance
    %% calculate average emotion intensity
    emotionIntensity(i_subject) = mean(mean(subjectData));
    %% calculate extreme scale use
    midpointDistance = abs(subjectData-scaleMidpoint);
    extremeScaleUse(i_subject) = mean(mean(midpointDistance))/scaleMidpoint; % calculate as ratio to possible distance from scale midpoint
end

%% compile and save off matrix for construct validity checks
singleEmotion = singleEmotion';
emotionIntensity = emotionIntensity';
extremeScaleUse = extremeScaleUse';
NookMeasures = [singleEmotion emotionIntensity extremeScaleUse];
filename = ['dataSet' num2str(dataSet) '_NookMeasures.mat'];
save(filename,'NookMeasures');