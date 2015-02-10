clear;clc;close all
dataFname = 'Altimeter 1 Data - Subscale - Dec 13';
datalist = xlsread(dataFname); % gets data from the excel file of the alt data  

time = datalist(5:end,1);
altitude = datalist(5:end,2);
altitude = altitude * 0.3048; % m
altitude_filtered = zeros(length(altitude),1);
%% Filter altitude data
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

for j=1:length(altitude)
     % Average buffer 
     if abs(j) > 500
         altitude_filtered(j) = 0;
     else
         altitude_filtered(j) = altitude(j);
     end
     
     
end
    mask3=(abs(v)>500);
    %mask4 = [mask3(2:end);mask3(1)];
    v(mask3)=0;
    %altitudeFit = fit(time, v,'exp2' );
    %plot(altitudeFit);


data = [time v];

plot(data(:,1),data(:,2));
