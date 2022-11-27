clear all
close all
clc

T_array= linspace(10,290, 9); % [kelvin]
P_in_eV=1e-3; % [eV]

responsivity_array = zeros(1,length(T_array));

i=1;
while i<length(T_array)
    [psic, Ec, z] = SchrodingerPoisson1D_CB_Kane_Main(T_array(i));
    [ni_matrix, Tijp_matrix, Tipj_matrix] = rate_equation_solver_function(psic, Ec, z,P_in_eV);
    J_photo = J_photo_function(ni_matrix, Tijp_matrix, Tipj_matrix);
    responsivity_array(i) = responsivity_function(P_in_eV,J_photo);
    i=i+1;
end

figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(T_array, responsivity_array*1e3, LineWidth=5);
grid("minor");
xlabel("Temperature (Kelvin)")
ylabel("Responsivity (mA/W)");



