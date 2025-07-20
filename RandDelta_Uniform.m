
if strcmp(pert_norm,'Fro')
    Del_Zvec_rand=unit_dist_samples_in_d_ball(n*T,1);
    Del_Xvec_rand=unit_dist_samples_in_d_ball(n*T,1);
    Del_Uvec_rand=unit_dist_samples_in_d_ball(m*T,1);
else
    Del_Zvec_rand=2*rand(n*T,1)-1;
    Del_Xvec_rand=2*rand(n*T,1)-1;
    Del_Uvec_rand=2*rand(m*T,1)-1;
end

Del_Z_rand=reshape(Del_Zvec_rand,[n,T])*D_norms(1)*pert_size;
Del_X_rand=reshape(Del_Xvec_rand,[n,T])*D_norms(2)*pert_size;
Del_U_rand=reshape(Del_Uvec_rand,[m,T])*D_norms(3)*pert_size;
    
Delta_rand=[Del_Z_rand; Del_X_rand; Del_U_rand];