clc
clear variables
close all
format compact

%% Initial Conditions (Relative State in LVLH frame)
omega_orb_radps = 0.001259;   % [rad/s] orbital angular velocity
m_chaser_kg = 100;            % [kg]

r0_rel_LVLH_m = [  50;   -20;   10 ];      % [m]
v0_rel_LVLH_mps = [ 0.02; 0.00; -0.01 ];   % [m/s]

x0_LVLH = [ r0_rel_LVLH_m;
            v0_rel_LVLH_mps ];

% initial condition

Tsim = 15000;