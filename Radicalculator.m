function [Q,U,e,NTU] = Radicalculator(R_l,R_h,R_t,Fin_h,FPI,w_Q,a_Q)
% Heat Transfer Performance of Radiator 
% We are using based on e-NTU method 

% R_l : Radi length [mm]
% R_h : Radi height [mm]
% R_t : Radi thickness [mm]
% FPI : Fin per Inch []
% w_Q : water volumetric flow rate [m^3/s]
% a_Q : air volumetric flow rate [m^3/s]

%% Constant value
a_Ti = 35; % Air inlet temperature [Celsius]
a_row = 1.146; % Air density [kg/(m^3)]
a_cp = 1.006; % Air Specific heat capacity [KJ/(kg*K)]
a_mu = 0.00001884; % Air Dynamic viscosity [Pa/s]
a_k = 0.02699; % Air thermal conductivity [W/(m*K)]

w_Ti = 90; % Water inlet temperature [degree]
w_row = 965.3; % Water density [kg/(m^3)]
w_cp = 4.21; % Water Specific heat capacity [KJ/(kg*K)]
w_mu = 0.000315; % Water Dynamic viscosity [Pa/s]
w_k = 0.67288; % Water thermal conductivity [W/(m*K)]

al_k = 170; % Aluminum 3003 thermal conductivity [W/(m*K)]

%% Definitions of the radiator geometry
B_h = R_h; % Core height [mm]
B_w = R_l; % Core Width [mm]
B_t = R_t; % Core thickness [mm]
N_r = 1; % Number of rows of tubes in the core depth dimension []
N_ct ; % Number of coolant tubes in one row []
N_p ; % Number of profiles (Fins in one row) []
N_f = FPI; % Number of Fins per inch [/inch]
F_t = 0.1; % Fin thickness [mm]
F_h = Fin_h; % Fin height [mm]
F_p ; % Fin pitch [mm]
R_f = 0.08; % Fin end radius [mm]
alp_f ; % Angle of fin [degree]
Y_l = B_w; % Coolant tube length [mm]
Y_cl = B_t; % Coolant tube cross-section length [mm]
Y_cw = 1.8; % Coolant tube cross-section width [mm]
Y_t = 0.127; % Coolant tube thickness [mm]
Y_p = F_h + Y_cw ; % Coolant tube pitch [mm]
R_t = Y_cw/2 ; % Coolant tube end radius [mm]
L_p = 1; % Louvre pitch [mm]
L_l = 4.5; % Louvre length [mm]
L_h = 1; % Louvre height [mm]

%% Adjust Radiator geometry
t_c = (B_h-F_h)/(Y_p); % Number of Tubes count []
f_c = t_c + 1; % Number of Fins count []
N_fm = N_f/25.4*10^(3); % Fin per metre [/m] 
F_p = 25.4/(N_f/2)/2; % Fin pitch [mm]
N_ct = floor(t_c); % Adjust N_ct []
N_p = floor(f_c); % Adjust N_p []
B_h = (F_h+Y_cw)*N_ct+F_h; % Adjust B_h [mm]
alp_f = atand((F_p-2*R_f)/(F_h-2*R_f)); % Adjust alp_f [degree]

%% Calculation of relevant heat transfer areas
F_l = pi*R_f + (F_h-2*R_f)/cosd(alp_f); % Fin length [mm]
A_frr = B_h*B_w; % Radiator core frontal area [mm^2]
A_frt = Y_cw*Y_l*N_ct; % Coolant tube frontal area [mm^2]
A_frf = F_t*F_l*N_fm*Y_l*N_p; % Fin frontal heat ransfer area [mm^2]
A_f = 2*B_t*F_l*N_fm*Y_l*N_p; % Fin heat transfer area [mm^2]
A_a = A_f+2*N_ct*Y_l*N_r*((Y_cl-2*R_t)+(pi*R_t)); % Total heat transfer area on the coolant side [mm^2]
A_c = (2*pi*(R_t-Y_t)+2*(Y_cl-2*R_t))*Y_l*N_ct*N_r; % Total heat transfer area on the coolant side [mm^2]
A_pa = A_frr - A_frf - A_frt; % Total Air pass area [mm^2]
A_pc = (pi*(R_t-Y_t)^2+(Y_cw-2*Y_t)*(Y_cl-2*R_t))*N_ct*N_r; % Total Coolant pass area [mm^2]

%% Calculation of Air side heat transfer coefficients
a_v = a_Q/(A_frr*10^    (-6)) ; % Air velocity in front of Radiator
Pr_a = (a_cp*10^3*a_mu)/a_k; % Prandtl Number of air []
Re_Lp = (a_row*a_v*L_p*10^(-3))/a_mu; % Reynolds Number of air by Louvre pitch []
j = 0.249*Re_Lp^(-0.42)*L_h^(0.33)*(L_l/F_h)^(1.1)*F_h^(0.26); % Colburn modulus j factor correlation []
Nu_avLp = 0.906*Re_Lp^(1/2)*Pr_a^(1/3); % Nussult Number of average louvre pitches []
h_a = (Nu_avLp*a_k)/(L_p*10^(-3)); % heat transfer coefficient on louvre pitch [W/(m^2*K)]

%% Calculation of air-side fin efficiency 
L = F_h /2; % effective fin length [mm]
m = sqrt((2*h_a)/(al_k*F_t)); % Fin efficincy parameter []
eta_f = tanh(m*L)/m*L; % Fin efficiency []
eta_o = 1 - (A_f/A_a)*(1 - eta_f); % Total surface efficiency []

%% Calculation of Water velocity in each tubes
k = N_ct; % Number of coolant tubes in one rows []
m_dot_w = w_row*w_Q; % The coolant mass flowrate [kg/s]
m_dot_ti = m_dot_w/k; % The coolant mass flowrate of each tubes [kg/s]
c_f = 0.09; % Flow non-uniformity coefficient []
m_dot_k = m_dot_ti+sqrt((c_f*m_dot_ti)^2/k); % Adjust m_dot_ti by c_f [kg/s]
Q_dot_k = m_dot_k/w_row; % Adjust coolant mass flowrate of each tubes [m^3/s]
w_v = Q_dot_k/(A_pc/N_ct*10^(-6)); % Water velocity of each tubes [m/s]

%% Calculation of Water side heat transfer coefficients
D_h = 4*A_pc/(A_c/Y_l)*10^(-3); % hydraulic diameter [m]
Pr_w = (w_cp*10^3*w_mu)/w_k; % Prandtl Number of Water []
Re_dh = (w_row*w_v*D_h)/w_mu; % Reynolds Number of water []
f = (0.79*log(Re_dh)-1.64)^(-2); % Friction factor in the tube []

if Re_dh < 2300 % Laminar
    Nu_c = 3.66+((0.0668*(D_h/Y_l)*Re_dh*Pr_w)/(1+0.04*((D_h/Y_l)*Re_dh*Pr_w)^(2/3))); % Nusselt number of water by Hansen equation []
elseif Re_dh >= 2300 && Re_dh <= 10000 % Transition
    Nu_c = ((f/8)*(Re_dh-1000)*Pr_w)/(1+12.7*sqrt(f/8)*(Pr_w^(2/3)-1)); % Nusselt number of water by Gnielinski equation []
else % Turbulent
    Nu_c = ((f/8)*Re_dh*Pr_w)/(1.07+12.7*sqrt(f/8)*(Pr_w^(2/3)-1)); % Nusselt number of water by Petukhow equation []
end

h_c = Nu_c*w_k/D_h; % heat transfer coefficient of coolant side [W/(m^2*K)]

%% Calculation of Overall heat transfer coefficient
U = 1/(((1/(eta_o*h_a*A_a*10^(-6)))+(1/(h_c*A_c*10^(-6))))*(A_frr*10^(-6))); % Overall heat transfer coefficient [W/(m^2*K)]

%% e-NTU Method
m_dot_a = a_row*a_Q; % The air mass flow rate [kg/s]
C_a = m_dot_a*a_cp; % Total heat capacities of air [kW/K]
C_w = m_dot_W*w_cp; % Total heat capacities of water [kW/K]

if C_a < C_w
    C_min = C_a;
    C_max = C_w;
else
    C_min = C_w;
    C_max = C_a;
end

C_r = C_min/C_max; % The heat Capacity rate ratio []

NTU = U*(A_frr*10^(-6))/(C_min*10^(3));
e = 1-exp(NTU^(0.22)/C_r*(exp(-1*C_r*NTU^(0.78))-1)); 

Q = e*C_min*(w_Ti-a_Ti); 

end

