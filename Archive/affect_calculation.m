clear;
clc;

%% specify dataset
dataSet = 88; % set to 18, 39, or 88

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

%% determine valence of sampled words
for i_word = 1:height(words) % define valence categories
if words.Valence(i_word) > 5 % derived based on the database mean for Warriner et al (2013)
    valence(i_word) = 1;
else
    valence(i_word) = 2;
end
end

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
    %% get mean and std for valence
    mPositive(i_subject) = mean(mean(subjectData(:,(valence==1))));
    mNegative(i_subject) = mean(mean(subjectData(:,(valence==2))));
    sdPositive(i_subject) = mean(std(subjectData(:,(valence==1))));
    sdNegative(i_subject) = mean(std(subjectData(:,(valence==2))));
end

%% save off matrix for construct validity checks
affect = [mPositive' mNegative' sdPositive' sdNegative'];
filename = ['dataSet' num2str(dataSet) '_affect.mat'];
save(filename,'affect');