% A script to solve rate equations for a system with N=5 state using 4th
% order Runge Kutta method. I will express the problem using matrices. 

% Useful source (starts from basics to final 4th order Runge Kutta method):
% URL: https://lpsa.swarthmore.edu/NumInt/NumIntIntro.html

% Lifetimes, should be calculated previously:

% For the current settings (8.11.2022), i.e. number of Poisson loops is
% NLoops=10, and input structure is as in "input_file_Saha_Kumar.m", the
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

T12= T_LO_code_for_rate_equation_solver(1,5,psic,Ec,z);
T12p= T_LO_code_for_rate_equation_solver(1,4,psic,Ec,z);
T1p2= T_LO_code_for_rate_equation_solver(3,5,psic,Ec,z);

T21= T_LO_code_for_rate_equation_solver(5,1,psic,Ec,z);
T21p= T_LO_code_for_rate_equation_solver(5,3,psic,Ec,z);
T2p1= T_LO_code_for_rate_equation_solver(4,1,psic,Ec,z);

T13= T_LO_code_for_rate_equation_solver(1,6,psic,Ec,z);
T13p= T_LO_code_for_rate_equation_solver(1,7,psic,Ec,z);
T1p3= T_LO_code_for_rate_equation_solver(3,6,psic,Ec,z);

T31= T_LO_code_for_rate_equation_solver(6,1,psic,Ec,z);
T31p= T_LO_code_for_rate_equation_solver(6,3,psic,Ec,z);
T3p1= T_LO_code_for_rate_equation_solver(7,1,psic,Ec,z);

T14= T_LO_code_for_rate_equation_solver(1,8,psic,Ec,z);
T14p= T_LO_code_for_rate_equation_solver(1,9,psic,Ec,z);
T1p4= T_LO_code_for_rate_equation_solver(3,8,psic,Ec,z);

T41= T_LO_code_for_rate_equation_solver(8,1,psic,Ec,z);
T41p= T_LO_code_for_rate_equation_solver(8,3,psic,Ec,z);
T4p1= T_LO_code_for_rate_equation_solver(9,1,psic,Ec,z);

T15= T_LO_code_for_rate_equation_solver(1,10,psic,Ec,z);
T15p= T_LO_code_for_rate_equation_solver(1,12,psic,Ec,z);
T1p5= T_LO_code_for_rate_equation_solver(3,10,psic,Ec,z);

T51= T_LO_code_for_rate_equation_solver(10,1,psic,Ec,z);
T51p= T_LO_code_for_rate_equation_solver(10,3,psic,Ec,z);
T5p1= T_LO_code_for_rate_equation_solver(12,1,psic,Ec,z);



T23= T_LO_code_for_rate_equation_solver(5,6,psic,Ec,z);
T23p= T_LO_code_for_rate_equation_solver(5,7,psic,Ec,z);
T2p3= T_LO_code_for_rate_equation_solver(4,6,psic,Ec,z);

T32= T_LO_code_for_rate_equation_solver(6,5,psic,Ec,z);
T32p= T_LO_code_for_rate_equation_solver(6,4,psic,Ec,z);
T3p2= T_LO_code_for_rate_equation_solver(7,5,psic,Ec,z);

T24= T_LO_code_for_rate_equation_solver(5,8,psic,Ec,z);
T24p= T_LO_code_for_rate_equation_solver(5,9,psic,Ec,z);
T2p4= T_LO_code_for_rate_equation_solver(4,8,psic,Ec,z);

T42= T_LO_code_for_rate_equation_solver(8,5,psic,Ec,z);
T42p= T_LO_code_for_rate_equation_solver(8,4,psic,Ec,z);
T4p2= T_LO_code_for_rate_equation_solver(9,5,psic,Ec,z);

T25= T_LO_code_for_rate_equation_solver(5,10,psic,Ec,z);
T25p= T_LO_code_for_rate_equation_solver(5,12,psic,Ec,z);
T2p5= T_LO_code_for_rate_equation_solver(4,10,psic,Ec,z);

T52= T_LO_code_for_rate_equation_solver(10,5,psic,Ec,z);
T52p= T_LO_code_for_rate_equation_solver(10,4,psic,Ec,z);
T5p2= T_LO_code_for_rate_equation_solver(12,5,psic,Ec,z);



T34= T_LO_code_for_rate_equation_solver(6,8,psic,Ec,z);
T34p= T_LO_code_for_rate_equation_solver(6,9,psic,Ec,z);
T3p4= T_LO_code_for_rate_equation_solver(7,8,psic,Ec,z);

T43= T_LO_code_for_rate_equation_solver(8,6,psic,Ec,z);
T43p= T_LO_code_for_rate_equation_solver(8,7,psic,Ec,z);
T4p3= T_LO_code_for_rate_equation_solver(9,6,psic,Ec,z);

T35= T_LO_code_for_rate_equation_solver(6,10,psic,Ec,z);
T35p= T_LO_code_for_rate_equation_solver(6,12,psic,Ec,z);
T3p5= T_LO_code_for_rate_equation_solver(7,10,psic,Ec,z);

T53= T_LO_code_for_rate_equation_solver(10,6,psic,Ec,z);
T53p= T_LO_code_for_rate_equation_solver(10,7,psic,Ec,z);
T5p3= T_LO_code_for_rate_equation_solver(12,6,psic,Ec,z);



T45= T_LO_code_for_rate_equation_solver(8,10,psic,Ec,z);
T45p= T_LO_code_for_rate_equation_solver(8,12,psic,Ec,z);
T4p5= T_LO_code_for_rate_equation_solver(9,10,psic,Ec,z);

T54= T_LO_code_for_rate_equation_solver(10,8,psic,Ec,z);
T54p= T_LO_code_for_rate_equation_solver(10,9,psic,Ec,z);
T5p4= T_LO_code_for_rate_equation_solver(12,8,psic,Ec,z);


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

% Time step:

h=1e-6;
t=linspace(0,1e2,1e8+1);

% Densities of states:

n1=zeros(1, length(t));
n2=zeros(1, length(t));
n3=zeros(1, length(t));
n4=zeros(1, length(t));
n5=zeros(1, length(t));

% % Initial values (setting 1):
% 
% n1(1)=2e17;
% n2(1)=1e7;
% n3(1)=1e4;
% n4(1)=1e1;
% n5(1)=1e-1;

% Initial values (setting 2):

n1(1)=2e17;
n2(1)=1e14;
n3(1)=1e11;
n4(1)=1e8;
n5(1)=1e5;

% Solution:

i=1;
while i<length(t)
    t0=t(i);
    Q=[n1(i) n2(i) n3(i) n4(i) n5(i)]';
    k1=A*Q;
    Q1=Q+k1*h/2;
    k2=A*Q1;
    Q2=Q+k2*h/2;
    k3=A*Q2;
    Q3=Q+k3*h;
    k4=A*Q3;
    Q_next=Q+((k1+2.*k2+2.*k3+k4)./6)*h;
    i=i+1;
    n1(i)=Q_next(1); 
    n2(i)=Q_next(2); 
    n3(i)=Q_next(3);
    n4(i)=Q_next(4);
    n5(i)=Q_next(5);
end

figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(t, n1, t, n2, t, n3, t, n4, t, n5, LineWidth=5);
grid("minor");
legend("n1", "n2", "n3", "n4", "n5",Location="best");




