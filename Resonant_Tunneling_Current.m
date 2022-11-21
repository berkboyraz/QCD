%Constants
e = 1.602176487E-19; % electron charge [C]
h = 6.62606896E-34; % Planck constant [J.s]
h_bar = h/(2*pi);
eps_0 = 8.854187817620E-12; % Vaccum dielectric constant [F/m]

%Parameters
d=0; % Spatial seperation between centroid of wavefunctions [m]
delta_d=0; % Detuning energy [J]
hbar_omega=0; % Coupling energy through barrier [J]
omega=hbar_omega./h_bar;
L_p=0; % Length of single period [m]
gamma=7.5*e; % Some parameter(?) [J]
n5=0; % Population of 5th state
n4=0; % Population of 4th state
T=0; %Temperature[K]
kB= 1.3806488E-23;% Boltzmann's constant [J/K]


%Tunneling current formula
J_tunnel= (e*d*(omega^2)*2*gamma*h_bar*(n5-n4*(exp(-delta_d/(kB*T)))))...
    /(4*gamma^2*L_p+delta_d^2);
