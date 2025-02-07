function [adj_Q] = Q_air(D_l, D_h, R_l, R_h, F_d, v_real)
%% Parameter 
a_row = 1.146; % Air density [kg/m^3]

%% Fan Spec calculator
T = readtable("Fan_spec.xlsx");

x1 = T{:,1};
y1 = T{:,2};

coeffs1 = polyfit(x1, y1, 1);

C_0 = coeffs1(1);
C_1 = coeffs1(2);

%% K_r calculator 
H = readtable("Radi_prediction.xlsx");

x2 = H{:,1};
y2 = H{:,2};

Kr1 = (2*y2)/(a_row*x2^2);

ln_x2 = log(Kr1);

coeffs2 = polyfit(ln_x2, y, 1);

m = coeffs2(1);
n = coeffs2(2);

%% Air mass flow rate calculator
x3 = zeros(1,20);
y3 = zeros(1,20);

A0 = D_l*D_h; % Duct Area [mm^2]
A1 = R_l*R_h; % Radi front Area [mm^2]
A3 = A1; % Radi Back Area [mm^2]
A4 = pi/4*F_d^2; % Fan Front Area [mm^2]
A5 = A4; % Fan Back Area [mm^2]

for i = 1:20
    V0 = i;
    V1 = A0/A1*V0;
    V3 = V1;
    V4 = A3/A4*V3;
    V5 = V4;

    P0 = 101325; % Atmosphere Pressure [Pa]
    P1 = P0 + (1/2)*a_row*(V0^2-V1^2); 
    K_r = m*ln(V1)+n;
    P3 = P1 - K_r*a_row/2*V1^2;
    P4 = P3 + (1/2)*a_row*(V3^2-V4^2);
    P5 = P0;

    x3 = V0*A0*10^(-6);
    y3 = P5-P4;
end

coeffs3 = polyfit(x3, y3, 2);

a = coeffs3(1);
b = coeffs3(2);
c = coeffs3(3);

%% Quadratic Formula calculation
a_cal = a;
b_cal = b-m;
c_cal = c_n;

solQ = (-b_cal+sqrt(b_cal^2-4*a_cal*c_cal))/(2*a_cal);
solP = m*solQ + n;
v_min = solQ/(A1*10^(-6));

%% Adjust to input parameter

Q_real = v_real*A0;

if Q_real<solQ
    adj_Q = solQ;
else
    adj_Q = Q_real;
end

end

