

flag_des=false;
Delta=zeros(size(D));
D_Del=D+Delta;
flag_end=false; % infeasibility
i_nan_end=10;

for i_ite=1:Ite_num
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

    % gradient step
    %gradD=gradD_func_LQR(K_Del,S_Del,L_Del,Lam1_Del,Lam2_Del,D_Del,Q,R,gamma,n,m,T,rho);
    gradD=gradD_func_LQR_anal(K_Del,S_Del,L_Del,Lam1_Del,Lam2_Del,D_Del,Q,R,gamma,n,m,T,rho,A,B,anal_mats);

    if ~all(isfinite(gradD))
        scaling_tmp=1;
        D_Del=D_Del-Delta;
        for i_nan=1:i_nan_end
            scaling_tmp=scaling_tmp/2;
            D_Del_tmp=D_Del+scaling_tmp*Delta;
            [K_Del,~,S_Del,L_Del,Lam1_Del,Lam2_Del]=ctrl_design_LQR(D_Del,n,m,Q,R,gamma);
            gradD=gradD_func_LQR_anal(K_Del,S_Del,L_Del,Lam1_Del,Lam2_Del,D_Del,Q,R,gamma,n,m,T,rho,A,B,anal_mats);
            if all(isfinite(gradD))
                break;
            end
            if i_nan==i_nan_end
                flag_end=true;
            end
        end
    end

    Delta=Delta+step_size*gradD;

    % projection step
    Delta=proj_D(Delta,D_norms,pert_size,pert_norm,n,m);

    D_Del=D+Delta;
end

[Z_Del,X_Del,U_Del]=sep_ZXU(D_Del,n,m);
[Del_Z,Del_X,Del_U]=sep_ZXU(D-D_Del,n,m);
