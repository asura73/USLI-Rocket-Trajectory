function verticalV(dataFname)  

datalist = xlsread(dataFname); % gets data from the excel file of the alt data  

time = datalist(5:end,1);
altitude = datalist(5:end,2);
altitude = altitude * 0.3048;

v = (altitude(2:end) - altitude(1:end-1))./(time(2:end) - time(1:end-1));
v = [0;v];
data = [time v];

plot(data(:,1),data(:,2))
end
