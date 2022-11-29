clear all
close all
clc

load("schrodinger_solver_output_300K_0V.mat");
c = 3e8; % [m/s]
% lambda=linspace(8e-6,12e-6, 10);
lambda = 8e-6;
w_pht = c/lambda*2*pi;
P_in_W_m2 = 1e7; % [W/m2] 

[ni_matrix, Tijp_matrix, Tipj_matrix] = rate_equation_solver_function(psic, Ec, z,P_in_W_m2, w_pht);
photocurrent_density = J_photo_function(ni_matrix, Tijp_matrix, Tipj_matrix);
responsivity = responsivity_function(P_in_W_m2,photocurrent_density);
