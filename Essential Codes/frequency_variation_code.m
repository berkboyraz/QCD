clear all
close all
clc

load("schrodinger_solver_output_300K_0V.mat");
c = 3e8; % [m/s]
% lambda=linspace(8e-6,12e-6, 10);
lambda = [8e-6 8.5e-6 9e-6 9.5e-6 10e-6 11e-6];
w_pht = (c*2*pi)./lambda;
P_in_W_m2 = 10; % [W/m2] 

photocurrent_density=zeros(length(lambda),1);
responsivity=zeros(length(lambda),1);
i=1;
for w_pht_ = w_pht
[ni_matrix, Tijp_matrix, Tipj_matrix] = rate_equation_solver_function(psic, Ec, z,P_in_W_m2, w_pht_);
photocurrent_density(i) = J_photo_function(ni_matrix, Tijp_matrix, Tipj_matrix);
responsivity(i) = responsivity_function(P_in_W_m2,photocurrent_density(i));
i=i+1;
end

