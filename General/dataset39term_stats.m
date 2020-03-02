clear;
clc;

rawData = importdata('39ESdata_daily_filtered.xlsx');
allData = rawData.data.Sheet1;
subjectIDlist = unique(allData(:,1));

for i_subject = 1:length(subjectIDlist)
    subjectID = subjectIDlist(i_subject);
    index1 = find(allData(:,1)==subjectID);
    subjectData = allData(index1,2:end);
    dayIDlist = unique(subjectData(:,1));
    numDays(i_subject) = max(dayIDlist);
    for i_day = 1:length(dayIDlist)
        dayID = dayIDlist(i_day);
        index2 = find(subjectData(:,1)==dayID);
        numPrompts(i_day) = numel(index2);
    end
    meanPrompts(i_subject) = mean(numPrompts);
    sdPrompts(i_subject) = std(numPrompts);
    if i_subject == 1
        numPrompts_Total = numPrompts';
    else
        numPrompts_Total = [numPrompts_Total; numPrompts'];
    end
    clear numPrompts
end

numDays(29) = [];
meanPrompts(29) = [];
sdPrompts(29) = [];

minDays = min(numDays);
maxDays = max(numDays);
meanDays = mean(numDays);
sdDays = std(numDays);
gmeanPrompts = mean(meanPrompts);
gsdPrompts = mean(sdPrompts);
meanPrompts_Total = mean(numPrompts_Total);
sdPrompts_Total = std(numPrompts_Total);
