

flag_des=false;

[K_Del,~,S_Del,L_Del,Lam1_Del,Lam2_Del]=ctrl_design_LQR(D,n,m,Q,R,gamma);
gradD=gradD_func_LQR_anal(K_Del,S_Del,L_Del,Lam1_Del,Lam2_Del,D,Q,R,gamma,n,m,T,rho,A,B,anal_mats);
Delta=step_size*gradD;
Z_nor=D_norms(1); X_nor=D_norms(2); U_nor=D_norms(3);
Del_Z=E_Z*Delta; Del_X=E_X*Delta; Del_U=E_U*Delta;
Del_Z=sign(Del_Z).*Z_nor*pert_size;
Del_X=sign(Del_X).*X_nor*pert_size;
Del_U=sign(Del_U).*U_nor*pert_size;
Delta_DGSM=[Del_Z;Del_X;Del_U];
D_Del_DGSM=D+Delta_DGSM;

[Z_Del,X_Del,U_Del]=sep_ZXU(D_Del_DGSM,n,m);
[Del_Z,Del_X,Del_U]=sep_ZXU(D-D_Del_DGSM,n,m);

[K_Del_DGSM,~,~,~,~,~]=ctrl_design_LQR(D_Del_DGSM,n,m,Q,R,gamma);
if rho(K_Del_DGSM)>1
    flag_des=true;
end