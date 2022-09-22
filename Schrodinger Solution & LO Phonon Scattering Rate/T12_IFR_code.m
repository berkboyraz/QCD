
m0 = 9.10938188E-31; % electron mass [kg]
e = 1.602176487E-19; % electron charge [C]
h = 6.62606896E-34; % Planck constant [J.s]
h_bar = h/(2*pi);
eps_0 = 8.854187817620E-12; % Vaccum dielectric constant [F/m]

m_star_l = 0.07*m0; % conduction band effective mass (material dependent)

% The following three lines are necessary if the wave functions  and
% z-range will be supplied by a file. In the case of using
% "cascade_script.m", these lines are not necessary. load("psic.mat"); %
% wave-functions of states in quantum wells load("z.mat"); % spatial range
% along z-direction load("Ec.mat"); % energy values of states in quantum
% wells

psi_1=psic(:,1); % wave-function of 1st state 
psi_2=psic(:,2); % wave-function of 2nd state  

delta_Ec= (Ec(2)-Ec(1))*e; % energy difference between 1st and 2nd states, J

DELTA_1 = 5.67e-10; % mean height at the 1st interface, m
DELTA_2 = 5.67e-10; % mean height at the 2nd interface, m
LAMBDA_1 = 33e-10; % lateral size of Gaussian fluctiation at the 1st interface, m
LAMBDA_2 = 33e-10; % lateral size of Gaussian fluctiation at the 2nd interface, m

z_1_index = 2002;
z_2_index = 2504;

delta_U_1 = abs(V0(z_1_index)-V0(z_1_index+1))*e; % potential change at the 1st interface, V
delta_U_2 = abs(V0(z_2_index)-V0(z_2_index-1))*e; % potential change at the 2nd interface, V

sum_1 = (DELTA_1^2)*(LAMBDA_1^2)*(delta_U_1^2)*((psi_1(z_1_index)*psi_2(z_1_index))^2)*exp((-1*(LAMBDA_1^2)*m_star_l*delta_Ec)/(2*(h_bar^2)));
sum_2 = (DELTA_2^2)*(LAMBDA_2^2)*(delta_U_2^2)*((psi_1(z_2_index)*psi_2(z_2_index))^2)*exp((-1*(LAMBDA_2^2)*m_star_l*delta_Ec)/(2*(h_bar^2)));

T_12=(h_bar^3)/(pi*m_star_l*(sum_1+sum_2))*1e12;

display("The T_12 is :"+T_12+" ps");

