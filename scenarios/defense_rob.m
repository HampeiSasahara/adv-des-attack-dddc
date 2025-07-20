
if n==1
    rho=@(K) abs(A+B*K);
else
    rho=@(K) spec_rad_clsys(A,B,K);
end

epsrob_array=10.^eps_exp_array;
eps_length=length(epsrob_array);
J_array=zeros(1,eps_length);
data_gen;
gamma=1e-8;
[~,~,~,~,~,~,J_LQR]=ctrl_design_LQR(D,n,m,Q,R,gamma);

for i_eps=1:eps_length
    epsrob=epsrob_array(i_eps);
    pert_size=epsrob;
    data_gen;
    step_size=max(D_norms)*pert_size;
    %OptDelta_GradProj;
    %flag_des
    [~,~,~,~,~,J_robust]=ctrl_design_Robust(D,n,m,Q,R,epsrob,D_norms);
    J_array(i_eps)=J_robust/J_LQR;
end
J_array