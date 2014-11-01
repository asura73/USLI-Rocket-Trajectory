function xvecvec = traject

clear
clc
close all

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

Pgl = 1;
rhogl = 1;

% rocket parameters

dia = 4.03*25.4/1000; % m
area = pi * (dia/2)^2; % m^2
m0 = 10.7 * 0.45; % kg
cd_r = 0.5; % rocket cd
pin = 1*in2m; % pin extention
boost = 1; % start with motor on

% state space

rvec0 = [0;0;0]; % position in NWZ
rvecd0 = [0;0;0]; % velocity in NWZ
avec0 = [0;5;0].*d2r; % orientation

xvec0 = [rvec0;rvecd0;avec0;m0];

% setup ode45

odeOptions = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@events);
[TOUT, xvecOUT] = ode45(@ascent, [0 1000], xvec0, odeOptions);

x_r = xvecOUT(:,1);
y_r = xvecOUT(:,2);
z_r = xvecOUT(:,3);

x_r = x_r * m2ft;
y_r = y_r * m2ft;
z_r = z_r * m2ft;

apogee = z_r(end)

v_r = xvecOUT(:,6);

v_r = v_r * m2ft;

m_r = xvecOUT(:,10);
%plot(TOUT,z_r)
%plot3(x_r,y_r,z_r)
%axis equal
end

function xvecd = ascent(t, xvec)

global boost;
global pin;

if t > 6.4
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
cD_p = 0.5;

% sum forces on rocket in rocket frame


% drag = .5 * rho * S * V^2 * cD

drag_r = .5*rho*.0082*norm(rvecd)^2*cD_r;
drag_p = 0;

v = norm(rvecd);
if v ~= 0
    rvecdhat = rvecd ./ v;
else
    rvecdhat = [0;0;0];
end

dragvec = -1*(drag_r+drag_p)*rvecdhat;

if boost == 1
    %T = thrust(t);
    T = 150; % N
    %delm = .951 - .355;
    %6.34
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

function [value,isterminal,direction] = events(t,xvec)

    % end solution when vertical velocity equals zero
    
    value = xvec(6);
    isterminal = 1;
    direction = -1;
end

function parseThrust(filename)

    % take data file with motor name and populate global thrust 
    % vector
end

