function results = tableconcat_start(dataKey,data,variableNames)
%joins iterative data to overall key variable
results = horzcat(dataKey,data);
results.Properties.VariableNames = variableNames;
end