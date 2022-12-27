function [ni_matrix, Tijp_matrix, Tipj_matrix, Tij_matrix] =rate_equation_solver_function(psic, Ec, z, P_in_W_m2,w_pht,T)

% A script to solve rate equations for a system with N=5 state using 4th
% order Runge Kutta method. I will express the problem using matrices. 

% Useful source (starts from basics to final 4th order Runge Kutta method):
% URL: https://lpsa.swarthmore.edu/NumInt/NumIntIntro.html

% Lifetimes, should be calculated previously:

% For the current settings (27.11.2022), i.e. number of Poisson loops is
% NLoops=4, and input structure is as in "input_file_Saha_Kumar.m", the
% state labels and corresponding wavefunctions are as follows:

% psi(:,1)-->1
% psi(:,2)-->1++
% psi(:,3)-->1+
% psi(:,4)-->2+
% psi(:,5)-->2
% psi(:,6)-->3
% psi(:,7)-->3+
% psi(:,8)-->4
% psi(:,9)-->4+
% psi(:,10)-->5
% psi(:,11)-->5++
% psi(:,12)-->5+

% Considering the labels above, the lifetimes are calculated as:

T11p= T_LO_function(1,3,psic,Ec,z,T);
T1p1= T_LO_function(3,1,psic,Ec,z,T);

T12= T_LO_function(1,5,psic,Ec,z,T);
T12p= T_LO_function(1,4,psic,Ec,z,T);
T1p2= T_LO_function(3,5,psic,Ec,z,T);

T21= T_LO_function(5,1,psic,Ec,z,T);
T21p= T_LO_function(5,3,psic,Ec,z,T);
T2p1= T_LO_function(4,1,psic,Ec,z,T);

T13= T_LO_function(1,6,psic,Ec,z,T);
T13p= T_LO_function(1,7,psic,Ec,z,T);
T1p3= T_LO_function(3,6,psic,Ec,z,T);

T31= T_LO_function(6,1,psic,Ec,z,T);
T31p= T_LO_function(6,3,psic,Ec,z,T);
T3p1= T_LO_function(7,1,psic,Ec,z,T);

T14= T_LO_function(1,8,psic,Ec,z,T);
T14p= T_LO_function(1,9,psic,Ec,z,T);
T1p4= T_LO_function(3,8,psic,Ec,z,T);

T41= T_LO_function(8,1,psic,Ec,z,T);
T41p= T_LO_function(8,3,psic,Ec,z,T);
T4p1= T_LO_function(9,1,psic,Ec,z,T);

T15= 1/(1/T_LO_function(1,10,psic,Ec,z,T)+radiative_transition_rate_function(45,3.5,15e-3,Ec,w_pht,P_in_W_m2,psic,z,1,10));
% T15= T_LO_function(1,10,psic,Ec,z,T);
T15p= T_LO_function(1,12,psic,Ec,z,T);
T1p5= T_LO_function(3,10,psic,Ec,z,T);

% T51= 1/(T_LO_function(10,1,psic,Ec,z,T)+radiative_transition_rate_function(45,3.5,15e-3,Ec,w_pht,P_in_W_m2,psic,z,10,1));
T51= T_LO_function(10,1,psic,Ec,z,T);
T51p= T_LO_function(10,3,psic,Ec,z,T);
T5p1= T_LO_function(12,1,psic,Ec,z,T);

T22p= T_LO_function(5,4,psic,Ec,z,T);
T2p2= T_LO_function(4,5,psic,Ec,z,T);

T23= T_LO_function(5,6,psic,Ec,z,T);
T23p= T_LO_function(5,7,psic,Ec,z,T);
T2p3= T_LO_function(4,6,psic,Ec,z,T);

T32= T_LO_function(6,5,psic,Ec,z,T);
T32p= T_LO_function(6,4,psic,Ec,z,T);
T3p2= T_LO_function(7,5,psic,Ec,z,T);

T24= T_LO_function(5,8,psic,Ec,z,T);
T24p= T_LO_function(5,9,psic,Ec,z,T);
T2p4= T_LO_function(4,8,psic,Ec,z,T);

T42= T_LO_function(8,5,psic,Ec,z,T);
T42p= T_LO_function(8,4,psic,Ec,z,T);
T4p2= T_LO_function(9,5,psic,Ec,z,T);

T25= T_LO_function(5,10,psic,Ec,z,T);
T25p= T_LO_function(5,12,psic,Ec,z,T);
T2p5= T_LO_function(4,10,psic,Ec,z,T);

T52= T_LO_function(10,5,psic,Ec,z,T);
T52p= T_LO_function(10,4,psic,Ec,z,T);
T5p2= T_LO_function(12,5,psic,Ec,z,T);

T33p= T_LO_function(6,7,psic,Ec,z,T);
T3p3= T_LO_function(7,6,psic,Ec,z,T);

T34= T_LO_function(6,8,psic,Ec,z,T);
T34p= T_LO_function(6,9,psic,Ec,z,T);
T3p4= T_LO_function(7,8,psic,Ec,z,T);

T43= T_LO_function(8,6,psic,Ec,z,T);
T43p= T_LO_function(8,7,psic,Ec,z,T);
T4p3= T_LO_function(9,6,psic,Ec,z,T);

T35= T_LO_function(6,10,psic,Ec,z,T);
T35p= T_LO_function(6,12,psic,Ec,z,T);
T3p5= T_LO_function(7,10,psic,Ec,z,T);

T53= T_LO_function(10,6,psic,Ec,z,T);
T53p= T_LO_function(10,7,psic,Ec,z,T);
T5p3= T_LO_function(12,6,psic,Ec,z,T);

T44p= T_LO_function(8,9,psic,Ec,z,T);
T4p4= T_LO_function(9,8,psic,Ec,z,T);

T45= T_LO_function(8,10,psic,Ec,z,T);
T45p= T_LO_function(8,12,psic,Ec,z,T);
T4p5= T_LO_function(9,10,psic,Ec,z,T);

T54= T_LO_function(10,8,psic,Ec,z,T);
T54p= T_LO_function(10,9,psic,Ec,z,T);
T5p4= T_LO_function(12,8,psic,Ec,z,T);

T55p= T_LO_function(10,12,psic,Ec,z,T);
T5p5= T_LO_function(12,10,psic,Ec,z,T);

Tijp_matrix= [T11p T12p T13p T14p T15p; T21p T22p T23p T24p T25p; T31p T32p T33p T34p T35p; T41p T42p T43p T44p T45p; T51p T52p T53p T54p T55p];
Tipj_matrix= [T1p1 T1p2 T1p3 T1p4 T1p5; T2p1 T2p2 T2p3 T2p4 T2p5; T3p1 T3p2 T3p3 T3p4 T3p5; T4p1 T4p2 T4p3 T4p4 T4p5; T5p1 T5p2 T5p3 T5p4 T5p5];
Tij_matrix= [0 T12 T13 T14 T15; T21 0 T23 T24 T25; T31 T32 0 T34 T35; T41 T42 T43 0 T45; T51 T52 T53 T54 0];

% Formation of coefficient matrix:

a11= -(1/T12+1/T12p+1/T1p2+1/T13+1/T13p+1/T1p3+1/T14+1/T14p+1/T1p4+1/T15+1/T15p+1/T1p5);
a12= 1/T21+1/T21p+1/T2p1;
a13= 1/T31+1/T31p+1/T3p1;
a14= 1/T41+1/T41p+1/T4p1;
a15= 1/T51+1/T51p+1/T5p1;

a21= 1/T12+1/T12p+1/T1p2;
a22= -(1/T21+1/T21p+1/T2p1+1/T23+1/T23p+1/T2p3+1/T24+1/T24p+1/T2p4+1/T25+1/T25p+1/T2p5);
a23= 1/T32+1/T32p+1/T3p2;
a24= 1/T42+1/T42p+1/T4p2;
a25= 1/T52+1/T52p+1/T5p2;

a31= 1/T13+1/T13p+1/T1p3;
a32= 1/T23+1/T23p+1/T2p3;
a33= -(1/T31+1/T31p+1/T3p1+1/T32+1/T32p+1/T3p2+1/T34+1/T34p+1/T3p4+1/T35+1/T35p+1/T3p5);
a34= 1/T43+1/T43p+1/T4p3;
a35= 1/T53+1/T53p+1/T5p3;

a41= 1/T14+1/T14p+1/T1p4;
a42= 1/T24+1/T24p+1/T2p4;
a43= 1/T34+1/T34p+1/T3p4;
a44= -(1/T41+1/T41p+1/T4p1+1/T42+1/T42p+1/T4p2+1/T43+1/T43p+1/T4p3+1/T45+1/T45p+1/T4p5);
a45= 1/T54+1/T54p+1/T5p4;

a51= 1/T15+1/T15p+1/T1p5;
a52= 1/T25+1/T25p+1/T2p5;
a53= 1/T35+1/T35p+1/T3p5;
a54= 1/T45+1/T45p+1/T4p5;
a55= -(1/T51+1/T51p+1/T5p1+1/T52+1/T52p+1/T5p2+1/T53+1/T53p+1/T5p3+1/T54+1/T54p+1/T5p4);

A = [ a11 a12 a13 a14 a15 ; a21 a22 a23 a24 a25 ; a31 a32 a33 a34 a35 ;
    a41 a42 a43 a44 a45 ; a51 a52 a53 a54 a55 ];

% % % % Time step:
% % % 
% % % h=1e-25;
% % % t=linspace(0,2e-18,2e7+1);

% % % % Densities of states:
% % % 
% % % n1=zeros(1, length(t));
% % % n2=zeros(1, length(t));
% % % n3=zeros(1, length(t));
% % % n4=zeros(1, length(t));
% % % n5=zeros(1, length(t));

% % % % % Initial values (setting 1):
% % % % 
% % % % n1(1)=2e17;
% % % % n2(1)=1e7;
% % % % n3(1)=1e4;
% % % % n4(1)=1e1;
% % % % n5(1)=1e-1;
% % % 
% % % % Initial values (setting 2):
% % % 
% % % n1(1)=2e17;
% % % n2(1)=1e3;
% % % n3(1)=1e3;
% % % n4(1)=1e3;
% % % n5(1)=1e3;

% % % % Solution:
% % % 
% % % i=1;
% % % 
% % % Ah1 = h*A;
% % % Ah2 = Ah1 * Ah1;
% % % Ah3 = Ah2 * Ah1;
% % % Ah4 = Ah3 * Ah1;

% % % while i<length(t)
% % %     
% % % %     t0=t(i);
% % % 
% % %     Q=[n1(i) n2(i) n3(i) n4(i) n5(i)]';
% % % 
% % % %     k1=A*Q;
% % % %     Q1=Q+k1*h/2;
% % % %     k2=A*Q1;
% % % %     Q2=Q+k2*h/2;
% % % %     k3=A*Q2;
% % % %     Q3=Q+k3*h;
% % % %     k4=A*Q3;
% % % %     Q_next=Q+((k1+2.*k2+2.*k3+k4)./6)*h;
% % % 
% % % %     Q_next=Q+(((A*Q)+2.*(A*(Q+(A*Q)*h/2))+2.*(A*(Q+(A*(Q+(A*Q)*h/2))*h/2))+(A*(Q+(A*(Q+(A*(Q+(A*Q)*h/2))*h/2))*h)))./6)*h;
% % %     
% % %     Q_next = Q + Ah1*Q + (1/2)*Ah2*Q + (1/6)*Ah3*Q + (1/24)*Ah4*Q;
% % % 
% % %     i=i+1;
% % %     n1(i)=Q_next(1); 
% % %     n2(i)=Q_next(2); 
% % %     n3(i)=Q_next(3);
% % %     n4(i)=Q_next(4);
% % %     n5(i)=Q_next(5);
% % % end


% Adaptive Runge-Kutta-Fehlberg Method
% A lecture on the method: https://www.youtube.com/watch?v=9g_geBN08IY
% Initial parameters

tf=1e-4;
h=1e-15;
n_max=5e4;
e_min=1e0;
e_max=1e2;
h_min=1e-25;
h_max=1e-5;
t=zeros(1 ,n_max);
n1=zeros(1, n_max);
n2=zeros(1, n_max);
n3=zeros(1, n_max);
n4=zeros(1, n_max);
n5=zeros(1, n_max);
t(1)=0;
n1(1)=2e17;
n2(1)=1e3;
n3(1)=1e3;
n4(1)=1e3;
n5(1)=1e3;
i=1;

while i<n_max && t(i)<tf
    if h<h_min
        h=h_min;
    elseif h>h_max
        h=h_max;
    end

    Q=[n1(i) n2(i) n3(i) n4(i) n5(i)]';

    k_1 = A*Q;
    Q_1 = Q+k_1*h/4;
    k_2 = A*Q_1;
    Q_2 = Q+k_1*h*3/32+k_2*h*9/32;
    k_3 = A*Q_2;
    Q_3 = Q+k_1*h*1932/2197-k_2*h*7200/2197+k_3*h*7296/2197;
    k_4 = A*Q_3;
    Q_4 = Q+k_1*h*439/216-k_2*h*8+k_3*h*3680/513-k_4*h*845/4104;
    k_5 = A*Q_4;
    Q_next_5 = Q+k_1*h*25/216+k_3*h*1408/2565+k_4*h*2197/4104-k_5*h*1/5;

    Q_5 = Q-k_1*h*8/27+k_2*h*2-k_3*h*3544/2565+k_4*h*1859/4104-k_5*h*11/40;
    k_6 = A*Q_5;

    Q_next_6 = Q+k_1*h*16/135+k_3*h*6656/12825+k_4*h*28561/56430-k_5*h*9/50+k_6*h*2/55;

    error = max(abs(Q_next_5-Q_next_6));

    if error>e_max && h>h_min
        h=h/2;
    else 
        i=i+1;
        t(i)=t(i-1)+h;
        n1(i)=Q_next_6(1);
        n2(i)=Q_next_6(2);
        n3(i)=Q_next_6(3);
        n4(i)=Q_next_6(4);
        n5(i)=Q_next_6(5);

        if error<e_min
            h=h*2;
        end
    end
end

figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(t, n1, t, n2, t, n3, t, n4, t, n5, LineWidth=5);
grid("minor");
legend("n1", "n2", "n3", "n4", "n5",Location="best");

ni_matrix= [n1(i-1) n2(i-1) n3(i-1) n4(i-1) n5(i-1)];

end

