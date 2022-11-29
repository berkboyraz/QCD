clear all
close all
clc

T= 300; % [kelvin]
P_in_J_array = [0.002 0.004 0.006]; % [J/cm^2]
applied_voltage=0; % [V]

photocurrent_density_array = zeros(1,length(P_in_J_array));

[psic, Ec, z] = SchrodingerPoisson1D_CB_Kane_Main(T, applied_voltage);

iii=1;

while iii<length(P_in_J_array)+1
    display("Calculating for P_in = "+num2str(P_in_J_array(iii))+" eV");
    [ni_matrix, Tijp_matrix, Tipj_matrix] = rate_equation_solver_function(psic, Ec, z,P_in_J_array(iii));
    photocurrent_density_array(iii) = J_photo_function(ni_matrix, Tijp_matrix, Tipj_matrix);
    display(num2str(iii/length(P_in_J_array)*100)+"% of simulation is completed.");
    iii=iii+1;
end

figure('units','normalized','outerposition',[0 0 1 1]);
plot(P_in_J_array, abs(photocurrent_density_array), LineWidth=5);
grid("minor");
xlabel("Pin(W/m^2)")
ylabel("J-photo (A/m^2)");



