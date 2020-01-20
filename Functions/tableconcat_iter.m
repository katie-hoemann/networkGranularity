function results = tableconcat_iter(existingData,newData,index,colName)
%joins iterative data to overall key variable
results = horzcat(existingData,newData);
results.Properties.VariableNames{index+1} = colName;
end