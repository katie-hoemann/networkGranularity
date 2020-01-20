function results = savedata(dataTable,subjectCol,dataName,dataSet)
%saves array of data with first column identifying subject
filename = ['_' num2str(dataName) '.mat'];
results = table2array(dataTable)';
results = horzcat(subjectCol,results);
resultsData = results;
results = 1;
getname = @(x) inputname(1);
resultsString = getname(results);
results = resultsData;
save(['dataSet' num2str(dataSet) filename],resultsString);
end