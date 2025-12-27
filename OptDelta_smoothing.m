

flag_des=false;
Delta=zeros(size(D));
D_Del=D+Delta;
flag_end=false; % infeasibility
i_nan_end=10;
rho_pre=rho(K_dd);
eps_rho=0*1e-3;
step_size_init=step_size;
eps_Delta=max(D_norms)*pert_size;

for i_ite=1:Ite_num
    i_ite
    [K_Del,~,S_Del,L_Del,Lam1_Del,Lam2_Del]=ctrl_design_LQR(D_Del,n,m,Q,R,gamma);

    if ~all(isfinite(K_Del))
        scaling_tmp=1;
        D_Del=D_Del-Delta;
        for i_nan=1:i_nan_end
            scaling_tmp=scaling_tmp/2;
            D_Del_tmp=D_Del+scaling_tmp*Delta;
            [K_Del,~,S_Del,L_Del,Lam1_Del,Lam2_Del]=ctrl_design_LQR(D_Del,n,m,Q,R,gamma);
            if all(isfinite(K_Del))
                break;
            end
            if i_nan==i_nan_end
                flag_end=true;
            end
        end
    end

    if flag_end
        break;
    end

    if rho(K_Del)>1
        flag_des=true;
        break;
    end
    
    if i_ite ~= 1 && abs(rho_pre-rho(K_Del))<eps_rho
        break;
    else
        rho_pre=rho(K_Del);
    end

    % gradient step
    gradD_adj=gradD_func_LQR_adjoint_smoothing(K_Del,S_Del,L_Del,Lam1_Del,Lam2_Del,D_Del,Q,R,gamma,n,m,T,rho,A,B,anal_mats,i_ite*1e2);
    gradD_np=gradD_adj.detach().numpy();
    gradD=reshape(double(gradD_np),[nb,T]);

    if ~all(isfinite(gradD))
        scaling_tmp=1;
        D_Del=D_Del-Delta;
        for i_nan=1:i_nan_end
            scaling_tmp=scaling_tmp/2;
            D_Del_tmp=D_Del+scaling_tmp*Delta;
            [K_Del,~,S_Del,L_Del,Lam1_Del,Lam2_Del]=ctrl_design_LQR(D_Del,n,m,Q,R,gamma);
            gradD_adj=gradD_func_LQR_adjoint(K_Del,S_Del,L_Del,Lam1_Del,Lam2_Del,D_Del,Q,R,gamma,n,m,T,rho,A,B,anal_mats);
            gradD_np=gradD_adj.detach().numpy();
            gradD=reshape(double(gradD_np),[nb,T]);
            if all(isfinite(gradD))
                break;
            end
            if i_nan==i_nan_end
                flag_end=true;
            end
        end
    end

    step_size = 5*step_size_init/sqrt(i_ite);
    %step_size = 1e1*step_size_init/i_ite;
    Delta=Delta+step_size*gradD;

    % projection step
    Delta=proj_D(Delta,D_norms,pert_size,pert_norm,n,m);

    D_Del=D+Delta;

    if norm(Delta,Inf)<eps_Delta
        break;
    end
end

[Z_Del,X_Del,U_Del]=sep_ZXU(D_Del,n,m);
[Del_Z,Del_X,Del_U]=sep_ZXU(D-D_Del,n,m);
