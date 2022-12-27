function T_12= T_LO_function(state_1_index, state_2_index,psic,Ec,z,T)

m0 = 9.10938188E-31; % electron mass [kg]
e = 1.602176487E-19; % electron charge [C]
h = 6.62606896E-34; % Planck constant [J.s]
h_bar = h/(2*pi);
eps_0 = 8.854187817620E-12; % Vaccum dielectric constant [F/m]
kB   = 1.3806488E-23;               %% Boltzmann's constant [J/K]

m_star = 0.07*m0; % conduction band effective mass (material dependent)
h_bar_w_lo = e*36e-3; % optical phonon enery, J
eps_inf = 10.6; % high frequency relative permittivity (material dependent)
eps_static = 12.9; % relative static permittivity (material dependent)

psi_1=psic(:,state_1_index); % wave-function of 1st state
psi_2=psic(:,state_2_index); % wave-function of 2nd state

w_lo=h_bar_w_lo/h_bar; % optical phonon frequency

eps_p=1/((1/eps_inf)-(1/eps_static));
% eps_p=64;

Ec_1 = Ec(state_1_index)*e;
Ec_2 = Ec(state_2_index)*e;

% The case of phonon absorption
if Ec_1 - Ec_2 + h_bar_w_lo > 0
    Q = abs(((2*m_star/(h_bar^2))*(Ec_1-Ec_2+h_bar*w_lo))^(1/2));
    % The following is the integral to calculate I_12(k_l).
    integral_result=0;
    i=1;
    j=1;
    while i < length(z)+1
        while j < length(z)+1
            integral_result=integral_result+(1e-22)*psi_1(i)*psi_2(i)*psi_1(j)*psi_2(j)*exp(-Q*abs(z(i)-z(j)));
            % Here, dz*dz' (infinitesimal differantial variables) is represented with 1e-22.
            j=j+1;
        end
        j=1;
        i=i+1;
    end
    I_12=(1/Q)*integral_result;
    n = 1/((exp(h_bar_w_lo/(kB*T)))-1);
    T_12_abs=((4.*(h_bar^2).*eps_p*eps_0)./(m_star.*(e^2).*w_lo.*I_12))/n; % s
else
    T_12_abs=1e20;
end

% The case of phonon emission
if Ec_1 - Ec_2 - h_bar_w_lo > 0 
    Q = ((2*m_star/(h_bar^2))*(Ec_1-Ec_2-h_bar*w_lo))^(1/2);
    % The following is the integral to calculate I_12(k_l).
    integral_result=0;
    i=1;
    j=1;
    while i < length(z)+1
        while j < length(z)+1
            integral_result=integral_result+(1e-22)*psi_1(i)*psi_2(i)*psi_1(j)*psi_2(j)*exp(-Q*abs(z(i)-z(j)));
            % Here, dz*dz' (infinitesimal differantial variables) is represented with 1e-22.
            j=j+1;
        end
        j=1;
        i=i+1;
    end
    I_12=(1/Q)*integral_result;
    n = 1/((exp(h_bar_w_lo/(kB*T)))-1);
    T_12_emi=((4.*(h_bar^2).*eps_p*eps_0)./(m_star.*(e^2).*w_lo.*I_12))/(1+n); % s
else
    T_12_emi=1e20;
end

% LO phonon lifetime    
T_12= 1/(1/T_12_abs+1/T_12_emi);

end


