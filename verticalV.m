dataFname = 'Altimeter 1 Data - Subscale - Dec 13';
datalist = xlsread(dataFname); % gets data from the excel file of the alt data  

time = datalist(5:end,1)';
altitude = datalist(5:end,2)';
altitude = altitude * 0.3048;

v = diff(altitude)./diff(time);
v = [0 v];
data = [time;v];

plot(data(1,:),data(2,:));
