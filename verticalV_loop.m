clear;clc;close all
dataFname = 'Altimeter 1 Data - Subscale - Dec 13';
datalist = xlsread(dataFname); % gets data from the excel file of the alt data  

time = datalist(5:end,1);
altitude = datalist(5:end,2);
altitude = altitude * 0.3048; % m
altitude_filtered = zeros(length(altitude),1);

for i=1:length(altitude)
     % Average buffer 
     if altitude(i) < 0
         altitude_filtered(i) = altitude(i-1);
     else
         altitude_filtered(i) = altitude(i);
     end
end

v = diff(altitude_filtered)./diff(time);
v = [0;v];
velocity_filtered = zeros(length(v),1);
avgvelocity = zeros(length(v) - 1, 1);
time2 = zeros(length(v) - 1, 1);

for j=1:length(v)
     % Average buffer 
     if abs(v(j)) > 500
         velocity_filtered(j) = 0;
     else
         velocity_filtered(j) = v(j);
     end
     if j > 1
         time2(j) = (time(j) + time(j-1))/2;
         avgvelocity(j) = (velocity_filtered(j) + velocity_filtered(j-1))/2;
     else
         time2(j) = 0;
         avgvelocity(j) = 0;
     end
     
    
     
end


data = [time2 avgvelocity];

plot(data(:,1),data(:,2));
