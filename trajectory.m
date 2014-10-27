function [Time,Height,velocity]=trajectory

n = 1000;

%Initial Condition 

p0 = 2.1162*10^3; % Air Pressure at 0ft (lb/ft^2)
r0 = 2.3769*10^-3; % Air Density at 0ft (slugs/ft^3)
d = 4;            % Diameter of rocket (ft)
A = (d/2)^2*pi;     % Cross Sectional Area of rocket (ft^2)
Mp = 1.25156;     %Propellant Mass (lb)
Mi = 10.7;     %Initial Launch Mass (lb)
It = 213.4;      %Total Impulse (Pound-Seconds)
Tf = 6.40;      %Burn Time (s)
Ft = 33.4;      %Average Thrust (lb)
R = 1545.4;     % Gas constant (ft*lbf/(lbmol*R))
g = 32.2;
theta = 5/180*pi;

initial = [[0;0;0];[0;0;0];0;theta];
time = linspace(0,Tf,n); 


[Tb, Ub] = ode45(@boost, time, initial);

function res=boost(t, U)
    %Differential Equation Function for Solving with ode45
    %Inputs: Time, Input Vector
    %Input Vector: velocity,height,mass
    %Outputs: Output Vector
    %Output Vector: acceleration,velocity,mass rate of change
    
    
    %Unpack Vector
    Vr = U(1:3);  %Velocity
    Hr = U(4:6);  %Height
    Mr = U(7);  %Mass
    Ang = U(8);
    
    h = Hr(3);
    T = ((15-0.0065*h)+273)*1.8;  %temperature(R) at a given altitude from
                                  %Toussaint's Formula
    p = p0*exp(-32.2/(R*T))*h; %Pressure at a given altitude
    %Calculate Acceleration
    
    % Thrust = thrust * ihat
     
    A_rocket = (Ft/Mr -g-.5/Mr*Cd*p*A*Vr(3)^2)*Vr;
    Vr(3)=A_rocket*Tf/n;

    %Calculate Mass Rate of Change as Function of Current Thrust
    dM = -Ft*Mp/It;
    
    %Pack Result Vector
    res = [A_rocket; Vr; dM];
end


function ans = thrust(t)
%Returns Thrust of C6 Rocket at time t after ignition
%Uses Linear Interpolation from data on data sheet
%Inputs: t, Outputs: Thrust

    %Load Data from Data Sheet
    T = [0 0.0310 0.0920 0.1390 0.1920 0.2090 0.2310 0.2480 0.2920 ...
        0.3700 0.4750 0.6710 0.7020 0.7230 0.8500 1.0630 1.2110 1.2420 ...
        1.3030 1.4680 1.6560 1.8210 1.8340 1.8470 1.8600];
    F = [0 0.9460 4.8260 9.9360 14.0900 11.4460 7.3810 6.1510 5.4890 ...
        4.9210 4.4480 4.2580 4.5420 4.1640 4.4480 4.3530 4.3530 4.0690 ...
        4.2580 4.3530 4.4480 4.4480 2.9330 1.3250 0];
    
    %If Not in Range, Return No Thrust
    %If In Range, Return Linear Interpolation based on Data from Data Sheet
    if(t < T(1) || t > T(end))
        ans = 0;
    else
        ans = interp1(T,F,t);
    end
    
end
    
%%%%%%%%%%%%%%%%%%%
%%% COAST PHASE %%%
%%%%%%%%%%%%%%%%%%%

%Setup Initial Conditions and run ODE45 Solver
initial = Ub(end,:);
time = linspace(Tb(end),Tb(end)+10,n);
options = odeset('Events', @events);
[Tc, Uc] = ode45(@coast, time, initial,options);

function res=coast(t, U)
    %Differential Equation Function for Solving with ode45
    %Inputs: Time, Input Vector
    %Input Vector: velocity,height,mass
    %Outputs: Output Vector
    %Output Vector: acceleration,velocity,mass rate of change

    %Unpack Vector
    H_Dot = U(1);  %Velocity
    H     = U(2);  %Height
    M     = U(3);  %Mass

    %Get Current Thrust
    Ft = thrust(t);  %Rocket Thrust

    %Calculate Acceleration
    H_DDot = Ft/M - g - .5/M*Cd*p*A*H_Dot^2;

    %Pack Result Vector
    res = [H_DDot; H_Dot; 0];
end

function [value,isterminal,direction] = events(t,W)
    %Used by odeset to determine end condition for Coast Phase
    %End condition is Velocity = 0
    %Input: Time, [Velocity, Position, Mass]

    value = W(1);
    isterminal = 1;
    direction = -1;
end

%%%%%%%%%%%%%%%%%%%%%%
%%% COMBINE PHASES %%%
%%% & RETURN H(T)  %%%
%%%%%%%%%%%%%%%%%%%%%%

T  = [Tb;Tc];
Ht = [Ub(:,2);Uc(:,2)];

end