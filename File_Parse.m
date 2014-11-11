clear
clc


fh = fopen('J150.txt');

time=[];
mass = [];
force = [];

cnt = true;
while cnt == true
    line = fgetl(fh);
    if line==-1
        cnt = false;
    end
    
    if line~=-1
        [~,rest1] = strtok(line,'=');
        [t,~] = strtok(rest1,'=''"');
        t = str2double(t);
        time=[time;t];
        [~,rest2] = strtok(rest1,'='); 
        [f,~] = strtok(rest2,'=''"');
        f = str2double(f);
        force = [force;f];
        [~,rest3] = strtok(rest2,'='); 
        [m,~] = strtok(rest3,'=''"');
        m = str2double(m);
        mass = [mass;m];
    end
    
   
end
        
mass
force
time