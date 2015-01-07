function xvecvec = traject

clear
clc
hold off

d2r = pi / 180;
r2d = 1 / d2r;
in2m = 0.025;
m2in = 1 / in2m;
ft2m = .3048;
m2ft = 1 / ft2m;

% declare global vars

global pin;
global boost;

% atmosphere parameters

Pgl = 101325; %N/m^2
rhogl = 1.225; %kg/m^3

% rocket parameters

dia = 4.03*25.4/1000; % m
area = pi * (dia/2)^2; % m^2
m0 = 10.7 * 0.45; % kg
cd_r = 0.5; % rocket cd
pin = 1*in2m; % pin extention
boost = 1; % start with motor on

% parse engine data

parseThrust('I305.txt')
% state space

rvec0 = [0;0;0]; % position in NWZ
rvecd0 = [0;0;0]; % velocity in NWZ
avec0 = [0;5;0].*d2r; % orientation

xvec0 = [rvec0;rvecd0;avec0;m0];

% setup ode45

odeOptions = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@events);
[TOUT, xvecOUT] = ode45(@ascent, [0 1000], xvec0, odeOptions);

% global h_ideal;
% global time;
% 
% h_ideal = [xvecOUT(:,3), TOUT];
% 
% boost = 1;
% global control;
% global n;
% n = 1;
% control = zeros(1500,3);

% global extraBurn;
% extraBurn = 1;
% [TOUT_controlled, xvecOUT_controlled] = ode45(@controlled_ascent, [0 1000], xvec0, odeOptions);

% boost = 1;
% [TOUT_uncontrolled, xvecOUT_uncontrolled] = ode45(@uncontrolled_ascent, [0 1000], xvec0, odeOptions);

x_r = xvecOUT(:,1);
y_r = xvecOUT(:,2);
z_r = xvecOUT(:,3);

x_r = x_r * m2ft;
y_r = y_r * m2ft;
z_r = z_r * m2ft;
% z_r_c = xvecOUT_controlled(:,3) * m2ft;
% z_r_uc = xvecOUT_uncontrolled(:,3) * m2ft;

apogee = z_r(end)
plot(TOUT,z_r)
% v_r = xvecOUT(:,6);
% 
% v_r = v_r * m2ft;
% 
% m_r = xvecOUT(:,10);
% subplot(2,1,1)
% plot(TOUT,h_ideal(:,1)*m2ft,'k')
% hold on
% plot(TOUT_uncontrolled,z_r_uc,'b')
% plot(TOUT_controlled,z_r_c,'g')
% plot(TOUT_controlled,xvecOUT_controlled(:,6)*m2ft,'r')
% legend('nominal','uncontrolled','controlled','v')
% subplot(2,1,2)
% plot(control(:,3),control(:,1),'k')
% hold on
% plot(control(:,3),control(:,2),'r')
% legend('error','pins')

end

function xvecd = ascent(t, xvec)

global boost;
global pin;
global time;
global mvecd;
global Tvec;

if t > max(time)
    boost = 0;
end

d2r = pi / 180;

% xvec = [rvec ; rvecd ; avec ; m]
% xvecd = [rvecd ; rvecdd ; avecd ; md]

rvec = xvec(1:3);
rvecd = xvec(4:6);
avec = xvec(7:9);
m = xvec(10);


% get atmosphereic conditions

rho = 1.225;

g = 9.81;
% get rocket conditions

cD_r = 0.5;
cD_p = 1.25;

% sum forces on rocket in rocket frame


% drag = .5 * rho * S * V^2 * cD

servo_angle = 45;

area_p = -0.0000001*servo_angle^2 + 0.00002*servo_angle - 0.0002;
drag_r = .5*rho*.0045*norm(rvecd)^2*cD_r;
drag_p = .5*rho*area_p*norm(rvecd)^2*cD_p;


v = norm(rvecd);
if v ~= 0
    rvecdhat = rvecd ./ v;
else
    rvecdhat = [0;0;0];
end

dragvec = -1*(drag_r+drag_p)*rvecdhat;

if boost == 1
    %T = thrust(t);
    %T = 150; % N
    %delm = .951 - .355;
    %6.34
    [~,index] = min(abs(time-t));
    %md = mvecd(index);
    T = Tvec(index);
    md = -0.09400630915;
else
    T = 0;
    md = 0;
end

F_ijk  = [0;0;T] + dragvec;
F_nwz = rotmat(2,5*d2r)*F_ijk - [0;0;m*g];

% xvecd = [rvecd ; rvecdd ; avecd ; md]

% newton's second law
rvecdd = F_nwz ./ m;

xvecd = [rvecd; rvecdd; 0;0;0;md];

end

function xvecd = controlled_ascent(t, xvec)

global boost;
global pin;
global time;
global mvecd;
global Tvec;
global h_ideal;


% BOOST FOR AN EXTRA HALF SECOND... OOPSS 
    
if t > (max(time))
    boost = 0;
end

d2r = pi / 180;

% xvec = [rvec ; rvecd ; avec ; m]
% xvecd = [rvecd ; rvecdd ; avecd ; md]

rvec = xvec(1:3);
rvecd = xvec(4:6);
avec = xvec(7:9);
m = xvec(10);


% get atmosphereic conditions

rho = 1.225;

g = 9.81;
% get rocket conditions

cD_r = 0.5;
cD_p = 1.25;

% sum forces on rocket in rocket frame


% drag = .5 * rho * S * V^2 * cD


% change the servo angle based on the error in servo angle

[~,index] = min(abs(h_ideal(:,2)-t));
z_ideal = h_ideal(index,1);

err = rvec(3) - z_ideal;
gain = 0.05;

servo_angle = err*gain*norm(rvecd) + 45;
servo_angle = 90;
if servo_angle > 90
    servo_angle = 90;
elseif servo_angle < 0
    servo_angle = 0;
end


area_p = -0.0000001*servo_angle^2 + 0.00002*servo_angle - 0.0002;
drag_r = .5*rho*.0045*norm(rvecd)^2*cD_r;
drag_p = .5*rho*area_p*norm(rvecd)^2*cD_p;

v = norm(rvecd);
if v ~= 0
    rvecdhat = rvecd ./ v;
else
    rvecdhat = [0;0;0];
end

dragvec = -1*(drag_r+drag_p)*rvecdhat;

if boost == 1
    %T = thrust(t);
    %T = 150; % N
    %delm = .951 - .355;
    %6.34
    
    [~,index] = min(abs(time-t));
    %md = mvecd(index);
    T = Tvec(index);
    md = -0.09400630915;
else
    T = 0;
    md = 0;
end

% error in burn 
global extraBurn
if t > 4 && t < (4+extraBurn)
    T = T + 40;
end

F_ijk  = [0;0;T] + dragvec;
F_nwz = rotmat(2,5*d2r)*F_ijk - [0;0;m*g];

global n;
global control;
control(n,1) = err;
control(n,2) = servo_angle;
control(n,3) = t;
n = n + 1;

% xvecd = [rvecd ; rvecdd ; avecd ; md]

% newton's second law
rvecdd = F_nwz ./ m;

xvecd = [rvecd; rvecdd; 0;0;0;md];

end

function xvecd = uncontrolled_ascent(t, xvec)

global boost;
global pin;
global time;
global mvecd;
global Tvec;
global h_ideal;


% BOOST FOR AN EXTRA HALF SECOND... OOPSS 
    
if t > (max(time))
    boost = 0;
end

d2r = pi / 180;

% xvec = [rvec ; rvecd ; avec ; m]
% xvecd = [rvecd ; rvecdd ; avecd ; md]

rvec = xvec(1:3);
rvecd = xvec(4:6);
avec = xvec(7:9);
m = xvec(10);


% get atmosphereic conditions

rho = 1.225;

g = 9.81;
% get rocket conditions

cD_r = 0.5;
cD_p = 1.25;

% sum forces on rocket in rocket frame


% drag = .5 * rho * S * V^2 * cD


% change the servo angle based on the error in servo angle

servo_angle = 45;

area_p = -0.0000001*servo_angle^2 + 0.00002*servo_angle - 0.0002;
drag_r = .5*rho*.0045*norm(rvecd)^2*cD_r;
drag_p = .5*rho*area_p*norm(rvecd)^2*cD_p;
v = norm(rvecd);

if v ~= 0
    rvecdhat = rvecd ./ v;
else
    rvecdhat = [0;0;0];
end

dragvec = -1*(drag_r+drag_p)*rvecdhat;

if boost == 1
    %T = thrust(t);
    %T = 150; % N
    %delm = .951 - .355;
    %6.34
    
    [~,index] = min(abs(time-t));
    %md = mvecd(index);
    T = Tvec(index);
    md = -0.09400630915;
else
    T = 0;
    md = 0;
end

% error in burn 
global extraBurn;
if t > 4 && t < (4+extraBurn)
    T = T + 40;
end

F_ijk  = [0;0;T] + dragvec;
F_nwz = rotmat(2,5*d2r)*F_ijk - [0;0;m*g];

% xvecd = [rvecd ; rvecdd ; avecd ; md]

% newton's second law
rvecdd = F_nwz ./ m;

xvecd = [rvecd; rvecdd; 0;0;0;md];

end

function [value,isterminal,direction] = events(t,xvec)

    % end solution when vertical velocity equals zero
    
    value = xvec(6);
    isterminal = 1;
    direction = -1;
end

function parseThrust(filename)
fh = fopen(filename);

global time;
global force;
global mass;
global Tvec;
global mvecd;

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

%define new finer time
time_new = linspace(0,max(time),1000);

%interpolate other variables to fill space of new time vector
Tvec = interp1(time,force, time_new);
mass = interp1(time,mass, time_new);

%need to transpose these, interp1 spits out row vector
%and we need a column vector

time_new = time_new';
Tvec = Tvec';
mass = mass';

%calculate average mass rate in new time space
mass1 = [0;mass(1:end-1)];
time1 = [0;time_new(1:end-1)];

mvecd = (mass-mass1)./(time_new-time1);
mvecd(1) = 0;
mvecd = mvecd/1000;
%assign new time space to time global variable
time = time_new;
        
    % take data file with motor name and populate global thrust 
    % vector
end





