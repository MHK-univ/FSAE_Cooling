%% 'Start.m' must be executed before starting the simulation. 
% Here, we define the desired variables, which serve as the basis for the simulation.
% I have included an example at the beginning. You can modify it as needed in the future

clear;

%% Setting the simulation time
% Enter the final time recorded in AIM or ECU-Master.

t_sim = 1697; % Unit [s]

%% Setting radiator geometry variables
% Enter the geometry variables of the radiator you want to use

R_h = 352 ; % Radiator Core Height [mm]
R_l = 302 ; % Radiator Core Length [mm]
R_t = 32 ; % Radiator Core Thickness [mm]
Fin_h = 4.85 ; % Fin height [mm]
FPI = 16 ; % Fin per inch [/inch]

save('Radi_parameter.mat', 'R_h', 'R_l', "R_t", "Fin_h", "FPI", '-mat');

