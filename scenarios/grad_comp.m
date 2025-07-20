
N_trials=10;
result=zeros(N_trials,3);
% 1: numerical differentiation
% 2: implicit differentiation (numdif for KKT)
% 3: implicit differentiation (analdif for KKT)

if n==1
    rho=@(K) abs(A+B*K);
else
    rho=@(K) spec_rad_clsys(A,B,K);
end

for i=1:N_trials
    disp(i);
    % data generation
    data_gen;

    % numerical differentiation
    rhoK_D=@(D) rho(ctrl_design_LQR(D,n,m,Q,R,gamma));
    tic
    drhodD_1=AutoDiff(rhoK_D,D);
    gradD_1=reshape(drhodD_1,[2*n+m,T]);
    result(i,1)=toc

    % implicit differentiation (numdif for KKT)
    tic
    [K_dd,~,S_dd,L_dd,Lam1_dd,Lam2_dd]=ctrl_design_LQR(D,n,m,Q,R,gamma);
    %gradD_2=gradD_func_LQR(K_dd,S_dd,L_dd,Lam1_dd,Lam2_dd,D,Q,R,gamma,n,m,T,rho);
    gradD_2=grad_func_LQR_KKTnum(K_dd,S_dd,L_dd,Lam1_dd,Lam2_dd,D,Q,R,gamma,n,m,T,rho,A,B,anal_mats);
    result(i,2)=toc

    % implicit differentiation (analdif for KKT)
    tic
    [K_dd,~,S_dd,L_dd,Lam1_dd,Lam2_dd]=ctrl_design_LQR(D,n,m,Q,R,gamma);
    gradD_3=gradD_func_LQR_anal(K_dd,S_dd,L_dd,Lam1_dd,Lam2_dd,D,Q,R,gamma,n,m,T,rho,A,B,anal_mats);
    result(i,3)=toc
end

M=mean(result);
SD=std(result);