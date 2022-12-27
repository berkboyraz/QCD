clear all
close all
clc

T = 300; % [kelvin]
c = 3e8; % speed of light in free space [m/s]    
P_in_W_m2 = 1e-12;
% [psic, Ec, z] = SchrodingerPoisson1D_CB_Kane_Main(300, 0.01);
wavelength = 9e-6;
w_pht = 2*pi*c./wavelength;

load("schrodinger_solver_output_300K_0V.mat");


[ni_matrix, Tijp_matrix, Tipj_matrix, Tij_matrix] = rate_equation_solver_function(psic, Ec, z, P_in_W_m2, w_pht, T);

J_photo = J_photo_function(ni_matrix, Tijp_matrix, Tipj_matrix);

responsivity = responsivity_function(P_in_W_m2,J_photo);

