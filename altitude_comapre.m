clear;clc;close all
dataFname = 'Altimeter 1 Data - Subscale - Dec 13';
dataFname2 = 'open_rocket_altitude';
datalist = xlsread(dataFname); % gets data from the excel file of the alt data  

time = datalist(5:end,1);
altitude = datalist(5:end,2);
altitude_filtered = zeros(length(altitude),1);

for i=1:length(altitude)
     % Average buffer 
     if altitude(i) < 0
         altitude_filtered(i) = altitude(i-1);
     else
         altitude_filtered(i) = altitude(i);
     end
     
     if i < length(altitude)
        if abs(altitude_filtered(i+1) - altitude_filtered(i) > 100) 7
            altitude(i+1) = altitude_filtered(i);
        end
     end
end

plot(time, altitude_filtered);

hold on

datalist2 = xlsread(dataFname2);
plot(datalist2(:,1),datalist2(:,2),'k');



