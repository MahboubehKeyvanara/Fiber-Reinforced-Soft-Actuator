% Numerical code for the SoRo course - Wageningan University and Research
% version: v01.01.01 230526
%
% Author: Mahboubeh Keyvanara (m.keyvanara@tue.nl) - 
%
% This code studies the relation between physical parameters of a fiber
% reinforced bending actuator and the input pressure of the actuator. The
% final results are presented in a plot. Please change the final plot based
% on the physical parameter you choose to study in line#35.

%% define the parameters of the robot:

% Physical parameters:
% Actuator length (mm):
l=160;
% Chamber radius (mm):
a=8.0;
% Wall thickness (mm):
t=2.0;
b=t;

% Material properties:
% Shear modulus:
mu=0.314*10^3;

% Atmospheric pressure (Kpa):
P_atm= 101;

%% 
% In this section the code finds the relation between inpout pressure and 
% physical properties of the robot using equilibrium equation over 
% the moments acting on the body of the robot.

% Calculating M_theta --> a moment caused by stress and stretch in the
% bottom layer of the actuator. Choose one the loops bellow:

% A loop to change wall thickness (t);
T=linspace(1,3,1000);

% A loop to change the chamber radius (a):
%  A=linspace(6,12,1000);

% A loop to change the actuator length (l):
%  L=linspace(80,180,1000);

% starting pressure for simulations:
P_in=0;

for n=1:length(T)
    % Xhoose the parameter to study vs pressure:
    t=T(n);
    % a=A(n);
    % l=L(n);

    % Choose the bending angle of the actuator:
    theta=pi/2;

    % Part 1 of the integration (do not change):
    M1 = 0;
    beta = linspace(0,b,1000);
    y_beta = beta .* theta/l+1;
    s_beta = mu * (y_beta-1 ./ (y_beta.^3));
    for i=2:length(beta)
        M1 = (s_beta(i) * beta(i) + ...
            s_beta(i-1) * beta(i-1)) * 0.5 * (beta(i) - beta(i-1)) + M1;
    end
    M1 = 2 * (a+t) * M1;

    %part2 of the integration (do not change):
    M2 = 0;
    tau = linspace(0,t,1000);
    MM = 0;
    phi = linspace(0,pi/2,length(tau));
    R = l/theta;
    y_phi = (R+b+sin(phi) .* (a+tau)) ./ R;
    s_phi = mu * (y_phi-1 ./ (y_phi.^3));
    for j=2:length(tau)
        M21 = 0;
        for i=2:length(phi)
            M21 = (s_phi(i) * ((a+tau(j))^2*sin(phi(i))+b*(a+tau(j)))+ ...
                s_phi(i-1) * ((a+tau(j))^2*sin(phi(i-1))+ ...
                b * (a+tau(j))))*0.5*(phi(i)-phi(i-1)) + M21;
        end
        MM(j) = M21; %#ok<SAGROW> 
        M2 = (MM(j)+MM(j-1))*0.5*(tau(j)-tau(j-1)) + M2;
    end
    % Total moment (do not change):
    M = M1 + 2 * M2;
    %%%%input pressure(do not change):
    P_in(n) = 6 * M / (4*a^3+3*a^2*pi*b); %#ok<SAGROW> 
    P1 = P_in + P_atm;
end

figure(1)
grid on
hold on
plot(T,P_in,'-')
xlabel('wall thickness (mm)')
ylabel('pressure (kpa)')
