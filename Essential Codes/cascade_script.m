clear all
close all
clc

T=300; % [kelvin]
P_in_eV=1e-3; % [eV]
[psic, Ec, z] = SchrodingerPoisson1D_CB_Kane_Main(300, 0.01);
[ni_matrix, Tijp_matrix, Tipj_matrix] = rate_equation_solver_function(psic, Ec, z,P_in_eV);

J_photo = J_photo_function(ni_matrix, Tijp_matrix, Tipj_matrix);

responsivity = responsivity_function(P_in_eV,J_photo);

