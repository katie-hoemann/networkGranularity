ID=ones(1,length(eod))';
eod_50=horzcat(ID*50,eod);
xlswrite('eod_50.xlsx',eod_50);
clear eod
clear ID
clc

ARIEOD_100118=vertcat(eod_24,eod_25,eod_26,eod_27,eod_28,eod_29,eod_32,eod_33,eod_34,eod_37,eod_38,eod_40,eod_43,eod_45,eod_46,eod_47,eod_49,eod_50);
xlswrite('ARIEOD_100118.xlsx',ARIEOD_100118);
