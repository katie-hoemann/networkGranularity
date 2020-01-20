clear;
clc;

rawData = importdata('39ESdata_daily_v2.xlsx');
subjectIDlist = unique(rawData.data(:,1)); 
allData = [rawData.data zeros(length(rawData.data),1)];

for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index1 = find(allData(:,1)==subjectID);
    subjectData = allData(index1,:);
    dayIDlist = unique(subjectData(:,3));
    for i_day = 1:length(dayIDlist)
        dayID = dayIDlist(i_day);
        index2 = find(subjectData(:,3)==dayID);
        subjectData(index2,end) = i_day;
    end
    if i_subject == 1
        allData_recompile = subjectData;
    else
        allData_recompile = vertcat(allData_recompile,subjectData);
    end
end

columnNames = {rawData.colheaders{:} 'Day'}; 
dataSet = array2table(allData_recompile,'VariableNames',columnNames); 
writetable(dataSet,'39ESdata_daily_relabeled.xlsx');