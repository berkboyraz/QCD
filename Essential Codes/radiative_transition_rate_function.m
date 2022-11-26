function radiative_transmission_rate = radiative_transition_rate_function(theta_degree,n_eff,gamma_eV,Ec,w_pht, P_in_cm_2, psic,z,state_i_index,state_f_index)

e = 1.602176487E-19; % electron charge [C]
% electromagnetic radiation and the  z’ direction
h = 6.62606896E-34; % Planck constant [J.s]
h_bar = h/(2*pi);
c= 3e8; % speed of light in free space [m/s]
eps_0 = 8.854187817620E-12; % Vaccum dielectric constant [F/m]

theta_radian = theta_degree/180*pi; % angle between electric field of the 
gamma_J=gamma_eV*e; % [J]

Ef_J = Ec(state_f_index); % [J]
Ei_J = Ec(state_i_index); % [J]
Ef_Ei_diff_J = abs(Ef_J-Ei_J) ; % [J]
psi_i = psic(:,1);
psi_f = psic(:,10);

% temporary parameters
w_pht = (Ef_Ei_diff_J)/h_bar;

P_in_m2 = P_in_cm_2*1e4; 
fi = P_in_m2/Ef_Ei_diff_J; % photon flux, number of incoming photon per second [1/s]

part_1= ((e^2)*gamma_J*w_pht*fi*(cos(theta_radian))^2)/(2*c*n_eff*eps_0);
part_2= ((Ef_Ei_diff_J)-h_bar*w_pht)^2+(0.5*gamma_J)^2;

integral_result=0;
i=1;
while i < length(z)+1
    integral_result=integral_result+abs((1e-11)*psi_i(i)*z(i)*psi_f(i))^2;
    % Here, dz (infinitesimal differantial variables) is represented with 1e-11.
    i=i+1;
end

radiative_transmission_rate = part_1*integral_result/part_2;

end