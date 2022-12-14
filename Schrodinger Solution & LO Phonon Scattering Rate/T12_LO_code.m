
m0 = 9.10938188E-31; % electron mass [kg]
e = 1.602176487E-19; % electron charge [C]
h = 6.62606896E-34; % Planck constant [J.s]
h_bar = h/(2*pi);
eps_0 = 8.854187817620E-12; % Vaccum dielectric constant [F/m]

m_star_l = 0.07*m0; % conduction band effective mass (material dependent)
h_bar_w_lo = e*36e-3; % optical phonon enery, J
eps_inf = 10.6; % high frequency relative permittivity (material dependent)
eps_static = 12.9; % relative static permittivity (material dependent)

% The following three lines are necessary if the wave functions  and z-range
% will be supplied by a file. In the case of using "cascade_script.m",
% these lines are not necessary. 
% load("psic.mat"); % wave-functions of states in quantum wells
% load("z.mat"); % spatial range along z-direction
% load("Ec.mat"); % energy values of states in quantum wells

psi_1=psic(:,1); % wave-function of 1st state 
psi_2=psic(:,2); % wave-function of 2nd state  

delta_Ec= (Ec(2)-Ec(1))*e; % energy difference between 1st and 2nd states, J

w_lo=h_bar_w_lo/h_bar; % optical phonon frequency

eps_p=1/((1/eps_inf)-(1/eps_static)); 
% eps_p=64; 

% The lower state electron momentum, based on k_u=0 assumption. 
k_l=(abs(w_lo/h_bar-delta_Ec/(h_bar^2))*2*m_star_l)^(1/2); 

% The following is the integral to calculate I_12(k_l).
integral_result=0;
i=1;
j=1;
while i < length(z)+1
    while j < length(z)+1
        integral_result=integral_result+(1e-22)*psi_1(i)*psi_2(i)*psi_1(j)*psi_2(j)*exp(-k_l*abs(z(i)-z(j)));
        % Here, dz*dz' (infinitesimal differantial variables) is represented with 1e-22.
        j=j+1;
    end
    j=1;
    i=i+1;
end

% integral_result = 0;
% for i = 1:length(z)
%     for j = 1:length(z)
%         dIij = psi_1(j) * psi_2(j) * exp(-k_l*abs(z(i)-z(j))) * psi_1(i) * psi_2(j) * (dz).^2;
%         integral_result = integral_result + dIij;
%     end
% end


% sum1=sum(psi_1.*psi_1)*1e-11
% sum2=sum(psi_2.*psi_2)*1e-11

I_12=(1/k_l)*integral_result;

T_12=((4.*(h_bar^2).*eps_p*eps_0)./(m_star_l.*(e^2).*w_lo.*I_12))*1e12; % ps
display("The T_12 is :"+T_12+" ps");

% PsiTPsii=0;
% PsiTPsij=0;
%  for q = 1:length(z)
%      PsiTPsii = PsiTPsii + psi_1(q) * ((Ec(2)-V0(q))/(Ec(2)-V0(q)+Eg(q))) * psi_1(q) * 1e-11;
%      PsiTPsij = PsiTPsij + psi_1(q) * ((Ec(1)-V0(q))/(Ec(1)-V0(q)+Eg(q))) * psi_1(q) * 1e-11;
%  end
% 
%  T_12_new=((4.*(h_bar^2).*eps_p*eps_0)./(m_star_l.*(e^2).*w_lo.*I_12))*1e12*(1+PsiTPsii)*(1+PsiTPsij); % ps
% 
% display("The T_12_new is :"+T_12_new+" ps");




