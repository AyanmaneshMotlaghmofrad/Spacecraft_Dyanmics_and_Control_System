% ES0
clc
clear variables
close all
format compact

%% parameters

I_B_kgm2 = diag([10 15 20]);

% Angular velocity of Body wrt Inertial, expressed in Body frame [rad/s]
omega_BN_B_Init_radps = [0.0; 0.1; 0.2];

q_NB_Init = 0.5*ones(4,1);

M_ext_B_Nm = 0;


% Simulation setting
T_sim = 60;