m0 = 9.10938188E-31; % electron mass [kg]
e0 = 1.602176487E-19; % electron charge [C]
McE=0.07;
h = 6.62606896E-34; % Planck constant [J.s]
hbar = h/(2*pi);
eps_0 = 8.854187817620E-12; % Vaccum dielectric constant [F/m]
hwlo = e*36e-3; 

Psii = psic(:,2);
Psij = psic(:,1);
dZ=
Qp = sqrt(2*McE(state2)*m0*e0 /hbar.^2 * (Ec(2)-Ec(1)-hwlo));

Iij = 0;
for n = 1:length(z)
     for m = 1:length(z)
         dIij = Psii(m) * Psij(m) * exp(-Qp*abs(z(n)-z(m))*1e-10) * Psii(n) * Psij(n) * (dZ*1e-10).^2;
         Iij = Iij + dIij;
     end
 end
 
PsiTPsii=0;
PsiTPsij=0;
 for q = 1:length(z)
     PsiTPsii = PsiTPsii + Psii(q) * ((Ec(2)-V0(q))/(Ec(2)-V0(q)+Eg(q))) * Psii(q) * dZ*1e-10;
     PsiTPsij = PsiTPsij + Psij(q) * ((Ec(1)-V0(q))/(Ec(1)-V0(q)+Eg(q))) * Psij(q) * dZ*1e-10;
 end
 
inversetauij = (McE(state2)*m0*(e0/hbar).^3 * hwlo/(4* eps0 * kp0))*(1e-10*Iij/(1e-10*Qp)) * (1+PsiTPsii)*(1+PsiTPsij);

tauij = 1e12/inversetauij;
