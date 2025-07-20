
% scenario parameters
eps_array=10.^eps_exp_array;
eps_length=length(eps_array);
des_result_DGSM=zeros(1,eps_length);
des_result=zeros(1,eps_length);
des_result_rand=zeros(1,eps_length);

for i_eps=1:eps_length
    % Perturbation Size
    i_eps_rev=eps_length+1-i_eps;
    pert_size=eps_array(i_eps_rev);
    if i_eps>1 && des_result(i_eps_rev+1)==0
        break;
    end

    for i_sample=1:num_sample
        % data generation & step size
        data_gen;
        step_size=max(D_norms)*pert_size;
        
        % Controller Design
        [K_dd,P_dd,S_dd,L_dd,Lam1_dd,Lam2_dd]=ctrl_design_LQR(D,n,m,Q,R,gamma);

        % DGSM
        DGSM_Delta_GradProj;
        if flag_des
            des_result_DGSM(i_eps_rev)=des_result_DGSM(i_eps_rev)+1;
        end
        
        % I-DGSM
        OptDelta_GradProj;
        if flag_des
            des_result(i_eps_rev)=des_result(i_eps_rev)+1;
        end
    
        % Random Perturbation
        RandDelta_Uniform;
        K_Del_rand=ctrl_design_LQR(D+Delta_rand,n,m,Q,R,gamma);
        if all(isfinite(K_Del_rand))
            if rho(K_Del_rand)>1
                des_result_rand(i_eps_rev)=des_result_rand(i_eps_rev)+1;
            end
        end
    
        fprintf('i_eps = %d/%d, i_sample = %d/%d\n',i_eps,eps_length,i_sample,num_sample);
        toc
    end
end
des_result_DGSM=des_result_DGSM/num_sample;
des_result=des_result/num_sample;
des_result_rand=des_result_rand/num_sample;

%% Plot

% Define the colors
color1 = [0, 0.4470, 0.7410];    % dark blue
color2 = [0.8500, 0.3250, 0.0980]; % vermilion
color3 = [0.4660, 0.6740, 0.1880]; % dark green

figure;
semilogx(eps_array,des_result_rand,':','LineWidth',p.line_width_thick,'Color', color3); hold on;
semilogx(eps_array,des_result_DGSM,'--','LineWidth',p.line_width_thick, 'Color', color2);
semilogx(eps_array,des_result,'-','LineWidth',p.line_width_thick, 'Color', color1);
legend({'Baseline','DGSM','I-DGSM'});
ylabel('ASR','FontSize',p.fs,'interpreter','latex');
xlabel('$\epsilon$','FontSize',p.fs,'interpreter','latex');
set(gca,'Fontsize',p.label_size);
if p.save
    p.pl_save.save_png('deseps');
end
%close;
