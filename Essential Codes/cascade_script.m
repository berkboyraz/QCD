clear all
close all
clc

T=300; % [kelvin]
P_in_eV=1e-3; % [eV]
[psic, Ec, z] = SchrodingerPoisson1D_CB_Kane_Main(300);
[ni_matrix, Tijp_matrix, Tipj_matrix] = rate_equation_solver(psic, Ec, z,P_in_eV);

J_photo = J_photo_solver(ni_matrix, Tijp_matrix, Tipj_matrix);

responsivity = responsivity_calculator(P_in_eV,J_photo);
