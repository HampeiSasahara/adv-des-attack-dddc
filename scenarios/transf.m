
eps_array=10.^eps_exp_array;
eps_length=length(eps_array);
des_result_transf=zeros(1,eps_length);

if n==1
    rho=@(K) abs(A+B*K);
else
    rho=@(K) spec_rad_clsys(A,B,K);
end

sys_ASR=strcat('/graphs/ASR/ASR_',sys_flag,'.fig');
fig = openfig(sys_ASR);
lines = findall(fig, 'Type', 'line');
des_result_full = get(lines(1), 'YData');
des_result_rand = get(lines(end), 'YData');
close(fig);

for i_eps=1:eps_length
    i_eps_rev=eps_length+1-i_eps;
    pert_size=eps_array(i_eps_rev);
    if i_eps>1 && des_result_transf(i_eps_rev+1)==0
        break;
    end

    for i_sample=1:num_sample
        % data generation & step size
        data_gen;
        D_ori=D;

        % I-DGSM without knoeledge of data
        data_gen;
        step_size=max(D_norms)*pert_size;
        OptDelta_GradProj;
        D_hyp_Del=D_ori+Delta;
        K_hyp=ctrl_design_LQR(D_hyp_Del,n,m,Q,R,gamma);
        if rho(K_hyp)>=1
            des_result_transf(i_eps_rev)=des_result_transf(i_eps_rev)+1;
        end
    
        fprintf('i_eps = %d/%d, i_sample = %d/%d\n',i_eps,eps_length,i_sample,num_sample);
        toc
    end

end

des_result_transf=des_result_transf/num_sample;


%% Plot

% Define the colors
color1 = [0, 0.4470, 0.7410];    % dark blue
color2 = [0.8500, 0.3250, 0.0980]; % vermilion
color3 = [0.4660, 0.6740, 0.1880]; % dark green
color4 = [0.4940, 0.1840, 0.5560];  % dark purple

figure;
semilogx(eps_array,des_result_rand,':','LineWidth',p.line_width_thick, 'Color', color3); hold on; grid on;
semilogx(eps_array,des_result_full,'-','LineWidth',p.line_width_thick, 'Color', color1);
semilogx(eps_array,des_result_transf,'-.','LineWidth',p.line_width_thick, 'Color', color4);
legend({'Baseline', 'Full Knowledge', 'Partial Knowledge'});
ylabel('ASR','FontSize',p.fs,'interpreter','latex');
xlabel('$\epsilon$','FontSize',p.fs,'interpreter','latex');
set(gca,'Fontsize',p.label_size);
if p.save
    p.pl_save.save_png('transf');
end