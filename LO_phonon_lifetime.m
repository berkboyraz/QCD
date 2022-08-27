m_0=9.10938356*10e-31; %Unit: kg
m_l=                 ; %Unit: kg
e= 1.60217662*10e-19 ; %Unit: Coulomb
hw_lo=35*10e-3*e ; % Unit: Joule
h_bar= 1.0545718*10e-34 ; % Unit: Joule.s
w_lo=hw_lo/h_bar;
eps_rho=69; % Unit: F.m^-1, Bundan emin değilim


k_l= sqrt(((hw_lo-(E_u-E_l))*2*m_l)/(h_bar^2));

I=;

t_ij_inverse= (m_l*e^2*w_lo*I)/(4*h_bar^2*eps_rho);
