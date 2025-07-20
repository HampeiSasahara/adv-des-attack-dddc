%% From Control Tutorials for MATLAB & SIMULINK 
% https://ctms.engin.umich.edu/CTMS/index.php?aux=Home
% https://www.amazon.com/Control-Tutorials-MATLAB-Simulink-Web-Based/dp/0201477009

%% Cruise Control (CC)
%  state: velocity
%  input: accelaration
%   m: vehicle mass 1000kg
%   b: damping coefficient 50 N.s/m
%  control objective: velocity tracking

%% Motor Position (MP)
%  state: motor position, motor speed, armature current
%  input: voltage source
%   J: moment of inertia of the rotor 3.2284e-6 kg.m^2
%   b: motor viscous riction constant 3.5077e-6 N.m.s
%   Kb: electromotive force constant 0.0274 V/rad/sec
%   Kt: motor torque constant 0.0274 N.m/Amp
%   R: electric resistance 4 Ohm
%   L: electric inductance 2.75e-6 H
%  control objective: motor position tracking (to, e.g., 1 rad)

%% Suspension System (SS)
%  state: body mass distance (from natural length), body mass velocity, suspension mass distance (from natural length), suspension mass velocity
%  input: active suspension between body mass and suspension mass
%   M1: 1/4 bus body mass 2500kg
%   M2: suspension mass 320kg
%   K1: spring constant of suspension system 80,000 N/m
%   K2: spring constant of wheel and tire 500,000 N/m
%   b1: damping constant of suspension system 350 N.s/m
%   b2: damping constant of wheel and tire 15,020 N.s/m
%  control objective: stabilization to origin

%% Inverted Pendulum (IP) mounted to motorized cart -> rocket attitude control
%  state: cart position, cart velocity, pendulum angle, pendulum angular velocity
%  input: force to cart
%   M: cart mass 0.5kg
%   m: pendulum mass 0.2kg
%   b: friction coefficient for cart 0.1N/m/sec
%   l: length to pendulum center to mass 0.3m
%   I: mass moment of inertia of pendulum 0.006kg.m^2
%  control objective: stabilization to origin

%% Aircraft Pitch (AP)
%  state: three dimensional angle deviation
%  input: elevator deflection angle
%  control objective: pitch angle tracking (to, e.g., 0.2rad)
%

%% Ball and Beam (BB) -> robot balancing
%  state: ball position, ball velocity, beam angular velocity, beam angular acceleration
%  input: torque to beam
%   m: ball mass 0.11kg
%   r: ball radius 0.015m
%   d: lever arm offset 0.03m
%   g: gravitational acceleration 9.8 m/s^2
%   L: beam length 1.0m
%   J: ball's moment of inertia 9.99e-6 kg.m^2
%  control objective: ball position tracking


% appropriate eps?
% BB: 80% destab. at eps=1e-5 with gamma=1e-3

Ts=0.1; % sampling period
switch sys_flag
    case "CC"
        m_CC=1000;
        b_CC=50;
        A_c=-b_CC/m_CC;
        B_c=1/m_CC;
    case "MP"
        J_MP=3.2284e-6;
        b_MP=3.5077e-6;
        K_MP=0.0274;
        R_MP=4;
        L_MP=2.75e-6;
        A_c = [0 1 0;
            0 -b_MP/J_MP K_MP/J_MP;
            0 -K_MP/L_MP -R_MP/L_MP];
        B_c = [0 ; 0 ; 1/L_MP];
    case "SS"
        m1_SS=2500;
        m2_SS=320;
        k1_SS=80000;
        k2_SS=500000;
        b1_SS=350;
        b2_SS=15020;
        A_c=[0                 1   0                                              0;
          -(b1_SS*b2_SS)/(m1_SS*m2_SS)   0   ((b1_SS/m1_SS)*((b1_SS/m1_SS)+(b1_SS/m2_SS)+(b2_SS/m2_SS)))-(k1_SS/m1_SS) -(b1_SS/m1_SS);
           b2_SS/m2_SS             0  -((b1_SS/m1_SS)+(b1_SS/m2_SS)+(b2_SS/m2_SS))                      1;
           k2_SS/m2_SS             0  -((k1_SS/m1_SS)+(k1_SS/m2_SS)+(k2_SS/m2_SS))                      0];
        B_c=[0                 0
           1/m1_SS              (b1_SS*b2_SS)/(m1_SS*m2_SS)
           0                -(b2_SS/m2_SS)
           (1/m1_SS)+(1/m2_SS)    -(k2_SS/m2_SS)];
    case "IP"
        M_IP=0.5;
        m_IP=0.2;
        b_IP=0.1;
        I_IP=0.006;
        g_IP=9.8;
        l_IP=0.3;
        p_IP=I_IP*(M_IP+m_IP)+M_IP*m_IP*l_IP^2; %denominator for the A and B matrices
        A_c = [0      1              0           0;
             0 -(I_IP+m_IP*l_IP^2)*b_IP/p_IP  (m_IP^2*g_IP*l_IP^2)/p_IP   0;
             0      0              0           1;
             0 -(m_IP*l_IP*b_IP)/p_IP       m_IP*g_IP*l_IP*(M_IP+m_IP)/p_IP  0];
        B_c = [     0;
             (I_IP+m_IP*l_IP^2)/p_IP;
                  0;
                m_IP*l_IP/p_IP];
    case "AP"
        A_c = [-0.313 56.7 0; -0.0139 -0.426 0; 0 56.7 0];
        B_c = [0.232; 0.0203; 0];
    case "BB"
        m_BB=0.111;
        R_BB=0.015;
        g_BB=-9.8;
        L_BB=1.0;
        d_BB=0.03;
        J_BB=9.99e-6;
        H_BB=-m_BB*g_BB/(J_BB/(R_BB^2)+m_BB);
        A_c = [0 1 0 0;
           0 0 H_BB 0;
           0 0 0 1;
           0 0 0 0];
        B_c = [0 0 0 1]';
    case "TT"
        A=[0.96 0 0;
            0.04 0.97 0;
            -0.04 0 0.9];
        B=[8.8 -2.3 0;
            0.2 2.2 4.9;
            -0.21 -2.2 1.9];
        n=size(B,1); m=size(B,2);
    case "QT"
        As=[28 32 28 32];
        as=[.071 .057 .071 .057];
        grav=981;
        hs=[12.6 13 4.8 4.9];
        ks=[3.14 3.29];
        %ks=[2 5];
        gams=[0.43 0.34];
        tmp=As./as;
        Tss=tmp.*sqrt(2*hs/grav);
        
        A_c=[-1/Tss(1) 0 As(3)/(As(1)*Tss(3)) 0;
            0 -1/Tss(2) 0 As(4)/(As(2)*Tss(4));
            0 0 -1/Tss(3) 0;
            0 0 0 -1/Tss(4)];
        B_c=[gams(1)*ks(1)/As(1) 0;
            0 gams(2)*ks(2)/As(2);
            0 (1-gams(2))*ks(2)/As(3);
            (1-gams(1))*ks(1)/As(4) 0];
end

if sys_flag~="TT"
    [n,m]=size(B_c);
    sys_c=ss(A_c,B_c,eye(n),0);
    sys_d=c2d(sys_c,Ts);
    A=sys_d.a; B=sys_d.b;
end


% Data-driven Control Parameters
Q=q*eye(n);
R=r*eye(m); V=sqrtm(R);
[K_ana,P_ana]=dlqr(A,-B,Q,R);

% Data Length
if strcmp(T_str,'short')
    T=n+m+1;
elseif strcmp(T_str,'medium')
    T=2*(n+m);
elseif strcmp(T_str,'long')
    T=3*(n+m);
end

% auxiliary function
if n==1
    rho=@(K) abs(A+B*K);
else
    rho=@(K) spec_rad_clsys(A,B,K);
end

% for analytical expression of gradients
nb=2*n+m;
anal_mats.E_Z=[eye(n) zeros(n,n) zeros(n,m)];
anal_mats.E_X=[zeros(n,n) eye(n) zeros(n,m)];
anal_mats.E_U=[zeros(m,n) zeros(m,n) eye(m)];
E_Z=anal_mats.E_Z;
E_X=anal_mats.E_X;
E_U=anal_mats.E_U;
barE1=[eye(n+m);zeros(2*n,n+m)];
barE2=[zeros(n+m,2*n);eye(2*n)];
anal_mats.E_1=kron(barE1,barE1);
anal_mats.E_2=kron(barE2,barE2);
anal_mats.kE_X=kron(eye(T),E_X);
anal_mats.kE_U=kron(eye(T),E_U);


tmp=[eye(m);zeros(n,m)];
anal_mats.dF1dS=kron(tmp,tmp);
anal_mats.dF2dS=zeros(4*n^2,m^2);
C_T_mpn=Create_ComMat(T,m+n);
C_T_2n=Create_ComMat(T,2*n);
C_n_mpn=Create_ComMat(n,m+n);
C_n_2n=Create_ComMat(n,2*n);
C_1=kron(kron(eye(n),C_T_mpn),eye(m+n));
C_1_tilde=kron(kron(eye(T),C_n_mpn),eye(m+n));
C_2=kron(kron(eye(n),C_T_2n),eye(2*n));
C_2_tilde=kron(kron(eye(T),C_n_2n),eye(2*n));

tmp1=[zeros(m,n);eye(n)]; tmp1=tmp1(:);
tmp2=[V*E_U;zeros(n,nb)];
tmp3=[zeros(m,nb);E_X];
C=Create_ComMat(T,n);
tmp4=kron(tmp1,eye((m+n)*T));
tmp5=kron(eye((m+n)*T),tmp1);
tmp6=kron(eye(T),tmp2);
ddF1dDdL_1=C_1*tmp4*tmp6;
ddF1dDdL_2=C_1*tmp4*(kron(eye(T),tmp3));
ddF1dDdL_3=kron(C',eye((m+n)^2))*C_1_tilde*tmp5*tmp6;
anal_mats.ddF1dDdL=ddF1dDdL_1+ddF1dDdL_2+ddF1dDdL_3;

tmp1=[eye(n);zeros(n,n)]; tmp1=tmp1(:);
tmp2=[E_X;zeros(n,nb)];
tmp3=[zeros(n,n);eye(n)]; tmp3=tmp3(:);
tmp4=[E_Z;zeros(n,nb)];
tmp5=kron(tmp3,eye(2*n*T));
tmp6=kron(eye(T),tmp4);
ddF2dDdL_1=C_2*kron(tmp1,eye(2*n*T))*kron(eye(T),tmp2);
ddF2dDdL_2=C_2*tmp5*tmp6;
ddF2dDdL_3=kron(C',eye(4*n^2))*C_2_tilde*kron(eye(2*n*T),tmp3)*tmp6;
ddF2dDdL_4=C_2*tmp5*kron(eye(T),[zeros(n,nb);E_X]);
anal_mats.ddF2dDdL=ddF2dDdL_1+ddF2dDdL_2+ddF2dDdL_3+ddF2dDdL_4;

anal_mats.ddFdDdL_1=kron(eye(n*T),anal_mats.E_1);
anal_mats.ddFdDdL_2=kron(eye(n*T),anal_mats.E_2);

anal_mats.dGamdX=kron(eye(T),[zeros(m,n);eye(n)]);
anal_mats.dGamdU=kron(eye(T),[eye(m);zeros(n,m)]);