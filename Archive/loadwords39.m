words39=readtable('words39.csv');

for i=1:height(words39)
if words39.Valence(i)>5
    ValCat(i)={'Positive'};
    Positive(i)=1;
else
    ValCat(i)={'Negative'};
    Positive(i)=0;
end
end

for i=1:height(words39)
if words39.Arousal(i)>4.6
    AroCat(i)={'High'};
    High(i)=1;
else
    AroCat(i)={'Low'};
    High(i)=0;
end
end

words39=[words39 ValCat' AroCat'];
words39.Properties.VariableNames(5:end)={'ValCat' 'AroCat'};
labels=[Positive' High'];

