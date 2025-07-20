

% scenario parameters
gam_array=10.^gam_exp_array;
gam_length=length(gam_array);
des_result=zeros(1,gam_length);
LQRper_result=zeros(num_sample,gam_length);
data_gen;
gamma=1e-8;
[~,~,~,~,~,~,J_LQR]=ctrl_design_LQR(D,n,m,Q,R,gamma);

for i_gam=1:gam_length
    % Perturbation Size
    gamma=gam_array(i_gam);
    %if i_gam>1 && des_result(i_gam-1)==0
    %    break;
    %end

    for i_sample=1:num_sample
        % data generation & step size
        data_gen;
        step_size=max(D_norms)*pert_size;
        
        % Controller Design
        [K_dd,P_dd,S_dd,L_dd,Lam1_dd,Lam2_dd]=ctrl_design_LQR(D,n,m,Q,R,gamma);
        
        % I-DGSM
        OptDelta_GradProj;
        if flag_des
            des_result(i_gam)=des_result(i_gam)+1;
        end

        
        % Performance Evaluation
        D_rand=D+pert_size*(2*rand(size(D))-1);
        [K_reg,~,~,~,~,~]=ctrl_design_LQR(D_rand,n,m,Q,R,gamma);
        while any(abs(eig(A+B*K_reg))>1)
            D_rand=D+pert_size*(2*rand(size(D))-1);
            [K_reg,~,~,~,~,~]=ctrl_design_LQR(D_rand,n,m,Q,R,gamma);
        end
        CGram=dlyap(A+B*K_reg,eye(n));
        LQRperformance=trace(Q*CGram)+trace(K_reg'*R*K_reg*CGram);
        LQRper_result(i_sample,i_gam)=LQRperformance/J_LQR;
    
        fprintf('i_gam = %d/%d, i_sample = %d/%d\n',i_gam,gam_length,i_sample,num_sample);
        toc
    end
end
des_result=des_result/num_sample;
LQRper_means=mean(LQRper_result);
LQRper_std=std(LQRper_result);

%% Plot

% Define the colors
color1 = [0, 0.4470, 0.7410];    % dark blue
color2 = [0.8500, 0.3250, 0.0980]; % vermilion
color3 = [0.4660, 0.6740, 0.1880]; % dark green

%figure;
%semilogx(gam_array,des_result,'-','LineWidth',p.line_width_thick, 'Color', color1);
%ylabel('ASR','FontSize',p.fs,'interpreter','latex');
%xlabel('$\gamma$','FontSize',p.fs,'interpreter','latex');
%set(gca,'Fontsize',p.label_size);


figure;
subplot(2,1,1);
semilogx(gam_array,des_result,'-','LineWidth',p.line_width_thick, 'Color', color1);
%legend({'Baseline','DGSM','I-DGSM'});
ylabel('ASR','FontSize',p.fs,'interpreter','latex');
%xlabel('$\gamma$','FontSize',p.fs,'interpreter','latex');
set(gca,'Fontsize',p.label_size); grid on;
subplot(2,1,2);
%tmpp=LQRper_means+LQRper_std; tmpn=LQRper_means-LQRper_std;
%semilogx(gam_array,tmpp,'-','LineWidth',p.line_width_thick, 'Color', color1); hold on;
%semilogx(gam_array,tmpn,'-','LineWidth',p.line_width_thick, 'Color', color1);
%patch([gam_array fliplr(gam_array)],[tmpp fliplr(tmpn)],'g');
semilogx(gam_array,LQRper_means,'-','LineWidth',p.line_width_thick, 'Color', color1); hold on;
%legend({'Baseline','DGSM','I-DGSM'});
ylabel('RCP','FontSize',p.fs,'interpreter','latex');
xlabel('$\gamma$','FontSize',p.fs+6,'interpreter','latex');
set(gca,'Fontsize',p.label_size); grid on;
if p.save
    p.pl_save.save_png('defense_reg');
end
%close;

[gam_exp_array; des_result; LQRper_means]
