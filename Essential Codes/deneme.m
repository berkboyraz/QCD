
load("schrodinger_solver_output_300K_0V.mat");

e = 1.602176487E-19; % electron charge [C]
h = 6.62606896E-34; % Planck constant [J.s]
h_bar = h/(2*pi);
c= 3e8; % speed of light in free space [m/s]
eps_0 = 8.854187817620E-12; % Vaccum dielectric constant [F/m]
n_eff=3.5;

psi_1=psic(:,1);
psi_5=psic(:,10);

% term_1 = 15e-3*e;
% term_2 = (Ec(10)-Ec(1))*e-h_bar*w_pht;
% 
% term_3 = term_1/(term_2^2-(term_1/2)^2);
% 
% term_4 = (e^2*w_pht*fi*(cos(45/180*pi))^2)/(2*c*n_eff*eps_0);
% 
% integral_result_sqroot=0;
% 
% i=1;
% while i < length(z)+1
%     integral_result_sqroot=integral_result_sqroot+(1e-11)*psi_1(i)*z(i)*psi_5(i);
%     % Here, dz (infinitesimal differantial variable) is represented with 1e-11.
%     i=i+1;
% end
% 
% term_5=integral_result_sqroot^2;
% 
% W_rad=term_4*term_5*term_3;
% 
% T_rad=1/W_rad;

wavelength= linspace(6e-6,13e-6,1e3);
wavelength_in_um=wavelength*10^6;
w_pht=2*pi*c./wavelength;

for i=1:length(wavelength)
rate_matrix(i)=radiative_transition_rate_function(45,3.5,0.015,Ec,w_pht(i),10, psic,z,1,10);
end

figure('units','normalized','outerposition',[0 0 1 1]);
plot(wavelength_in_um, rate_matrix, LineWidth=5);
xlabel('Wavelength in um','FontSize',15)
ylabel('Radiative Transition Rate in 1/s','FontSize',15)
grid("minor");

