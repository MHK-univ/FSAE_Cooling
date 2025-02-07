%% 

clear;

t_sim = 1697;

%% 

R_h = 352 ; % [mm]
R_l = 302 ; % [mm]
R_t = 32 ; % [mm]
Fin_h = 4.85 ; % [mm]
FPI = 16 ; % [mm]

Matrix_parameter = [R_h; R_l; R_t; Fin_h; FPI;];

save('Radi_parameter.mat', 'Matrix_parameter', '');

