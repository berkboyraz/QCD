clear all
close all
clc

T_array= linspace(100,300, 2); % [kelvin]
P_in_eV=1e-3; % [eV]
applied_voltage=0; % [V]

responsivity_array = zeros(1,length(T_array));

iii=1;

while iii<length(T_array)+1
    display("Calculating for T = "+num2str(T_array(iii))+" Kelvin");
    [psic, Ec, z] = SchrodingerPoisson1D_CB_Kane_Main(T_array(iii), applied_voltage);
    [ni_matrix, Tijp_matrix, Tipj_matrix] = rate_equation_solver_function(psic, Ec, z,P_in_eV);
    J_photo = J_photo_function(ni_matrix, Tijp_matrix, Tipj_matrix);
    responsivity_array(iii) = responsivity_function(P_in_eV,J_photo);
    display(num2str(iii/length(T_array)*100)+"% of simulation is completed.");
    iii=iii+1;
end

figure('units','normalized','outerposition',[0 0 1 1]);
plot(T_array, abs(responsivity_array).*1e3, LineWidth=5);
grid("minor");
xlabel("Temperature (Kelvin)")
ylabel("Responsivity (mA/W)");



