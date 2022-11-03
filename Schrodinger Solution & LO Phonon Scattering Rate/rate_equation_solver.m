% A script to solve rate equations for a system with N=5 state using 4th
% order Runge Kutta method. I will expEcress the problem using matrices. 

% Useful source (starts from basics to finak 4th order Runge Kutta method):
% URL: https://lpsa.swarthmore.edu/NumInt/NumIntIntro.html

% Lifetimes, should be calculated previously:

% T12=;
% T12p=;
% T1p2=;
% 
% T21=;
% T21p=;
% T2p1=;
% 
% T13=;
% T13p=;
% T1p3=;
% 
% T31=;
% T31p=;
% T3p1=;
% 
% T14=;
% T14p=;
% T1p4=;
% 
% T41=;
% T41p=;
% T4p1=;
% 
% T15=;
% T15p=;
% T1p5=;
% 
% T51=;
% T51p=;
% T5p1=;
% 
% 
% 
% T23=;
% T23p=;
% T2p3=;
% 
% T32=;
% T32p=;
% T3p2=;
% 
% T24=;
% T24p=;
% T2p4=;
% 
% T42=;
% T42p=;
% T4p2=;
% 
% T25=;
% T25p=;
% T2p5=;
% 
% T52=;
% T52p=;
% T5p2=;
% 
% 
% 
% T34=;
% T34p=;
% T3p4=;
% 
% T43=;
% T43p=;
% T4p3=;
% 
% T35=;
% T35p=;
% T3p5=;
% 
% T53=;
% T53p=;
% T5p3=;
% 
% 
% 
% T45=;
% T45p=;
% T4p5=;
% 
% T54=;
% T54p=;
% T5p4=;


% Formation of coefficient matrix:

% a11= -(1/T12+1/T12p+1/T1p2+1/T13+1/T13p+1/T1p3+1/T14+1/T14p+1/T1p4+1/T15+1/T15p+1/T1p5);
% a12= 1/T21+1/T21p+1/T2p1;
% a13= 1/T31+1/T31p+1/T3p1;
% a14= 1/T41+1/T41p+1/T4p1;
% a15= 1/T51+1/T51p+1/T5p1;
% 
% a21= 1/T12+1/T12p+1/T1p2;
% a22= -(1/T21+1/T21p+1/T2p1+1/T23+1/T23p+1/T2p3+1/T24+1/T24p+1/T2p4+1/T25+1/T25p+1/T2p5);
% a23= 1/T32+1/T32p+1/T3p2;
% a24= 1/T42+1/T42p+1/T4p2;
% a25= 1/T52+1/T52p+1/T5p2;
% 
% a31= 1/T13+1/T13p+1/T1p3;
% a32= 1/T23+1/T23p+1/T2p3;
% a33= -(1/T31+1/T31p+1/T3p1+1/T32+1/T32p+1/T3p2+1/T34+1/T34p+1/T3p4+1/T35+1/T35p+1/T3p5);
% a34= 1/T43+1/T43p+1/T4p3;
% a35= 1/T53+1/T53p+1/T5p3;
% 
% a41= 1/T14+1/T14p+1/T1p4;
% a42= 1/T24+1/T24p+1/T2p4;
% a43= 1/T34+1/T34p+1/T3p4;
% a44= -(1/T41+1/T41p+1/T4p1+1/T42+1/T42p+1/T4p2+1/T43+1/T43p+1/T4p3+1/T45+1/T45p+1/T4p5);
% a45= 1/T54+1/T54p+1/T5p4;
% 
% a51= 1/T15+1/T15p+1/T1p5;
% a52= 1/T25+1/T25p+1/T2p5;
% a53= 1/T35+1/T35p+1/T3p5;
% a54= 1/T45+1/T45p+1/T4p5;
% a55= -(1/T51+1/T51p+1/T5p1+1/T52+1/T52p+1/T5p2+1/T53+1/T53p+1/T5p3+1/T54+1/T54p+1/T5p4);
% 
% A = [ a11 a12 a13 a14 a15 ; a21 a22 a23 a24 a25 ; a31 a32 a33 a34 a35 ;
%     a41 a42 a43 a44 a45 ; a51 a52 a53 a54 a55 ];

% Time step:

h=1e-9;
t=linspace(0,1e-3,1e6+1);

% Densities of states:

n1=zeros(1, length(t));
n2=zeros(1, length(t));
n3=zeros(1, length(t));
n4=zeros(1, length(t));
n5=zeros(1, length(t));

% Initial values:

% n1(1)=1e10;
% n2(1)=1e8;
% n3(1)=1e6;
% n4(1)=1e4;
% n5(1)=1e2;

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




