

% data generation
data_gen;

% Controller Design
[K_dd,P_dd,S_dd,L_dd,Lam1_dd,Lam2_dd]=ctrl_design_LQR(D,n,m,Q,R,gamma);

step_size=max(D_norms)*pert_size; % gradient step size
if n==1
    rho=@(K) abs(A+B*K);
else
    rho=@(K) spec_rad_clsys(A,B,K);
end

% Adversarial Perturbation Generation
OptDelta_GradProj;

% Random Perturbation Generation
RandDelta_Uniform;
D_Del_rand=D+Delta_rand;
[K_Del_rand,~,S_Del_rand,L_Del_rand]=ctrl_design_LQR(D_Del_rand,n,m,Q,R,gamma);
sol_cons_mat=[zeros(n); eye(n); zeros(m,n)];

%% Plot
% timeseries
figure;
subplot(2,1,1); hold on; grid on; set(gca,'Fontsize',p.label_size);
plot([1:T],U(1,:),'-','LineWidth',p.line_width_thick,'Color',[0,0,0]);
yl=strcat('$u_',num2str(1),'$');
ylabel(yl,'FontSize',p.fs,'interpreter','latex');
ylim([-2,2]);
subplot(2,1,2); hold on; grid on; set(gca,'Fontsize',p.label_size);
plot([1:T],X(1,:),'-','LineWidth',p.line_width_thick,'Color',[0,0,0]);
yl=strcat('$x_',num2str(1),'$');
ylabel(yl,'FontSize',p.fs,'interpreter','latex');
ylim([-2,2]);
xlabel('time','FontSize',p.fs,'interpreter','latex');
if p.save
    p.pl_save.save_png('UX');
end
%close;

figure;
subplot(2,1,1); hold on; grid on; set(gca,'Fontsize',p.label_size);
plot([1:T],Del_U(1,:),'-','LineWidth',p.line_width_thick,'Color',[0,0,0]);
yl=strcat('$\Delta u_',num2str(1),'$');
ylabel(yl,'FontSize',p.fs,'interpreter','latex');
ylim([-2,2]);
%xlabel('time','FontSize',p.fs,'interpreter','latex');
subplot(2,1,2); hold on; grid on; set(gca,'Fontsize',p.label_size);
plot([1:T],Del_X(1,:),'-','LineWidth',p.line_width_thick,'Color',[0,0,0]);
yl=strcat('$\Delta x_',num2str(1),'$');
ylabel(yl,'FontSize',p.fs,'interpreter','latex');
ylim([-2,2]);
xlabel('time','FontSize',p.fs,'interpreter','latex');
if p.save
    p.pl_save.save_png('dUX');
end
%close;

figure;
subplot(2,1,1); hold on; grid on; set(gca,'Fontsize',p.label_size);
plot([1:T],U_Del(1,:),'-','LineWidth',p.line_width_thick,'Color',[0,0,0]);
yl=strcat('$u_{\Delta',num2str(1),'}$');
ylabel(yl,'FontSize',p.fs,'interpreter','latex');
ylim([-2,2]);
subplot(2,1,2); hold on; grid on; set(gca,'Fontsize',p.label_size);
plot([1:T],X_Del(1,:),'-','LineWidth',p.line_width_thick,'Color',[0,0,0]);
yl=strcat('$x_{\Delta',num2str(1),'}$');
ylabel(yl,'FontSize',p.fs,'interpreter','latex');
ylim([-2,2]);
xlabel('time','FontSize',p.fs,'interpreter','latex');
if p.save
    p.pl_save.save_png('UXD');
end
%close;


circle_lim=1.2;

% eigenvalues
dLam=eig(A+B*K_dd);
figure; hold on;
t_circle = linspace(0,2*pi,100);
plot(sin(t_circle),cos(t_circle),'k--');
xline(0); yline(0);
%scatter(real(dLam),imag(dLam),p.msize,'x','MarkerEdgeColor',[0 0.4470 0.7410],'LineWidth',p.mline_width);
scatter(real(dLam),imag(dLam),p.msize,'x','MarkerEdgeColor',[0 0 0],'LineWidth',p.mline_width);
xlim([-circle_lim,circle_lim]); ylim([-circle_lim,circle_lim]); axis square;
xlabel('Re$(\lambda)$','FontSize',p.fs,'interpreter','latex');
ylabel('Im$(\lambda)$','FontSize',p.fs,'interpreter','latex');
set(gca,'Fontsize',p.label_size);
if p.save
    p.pl_save.save_png('eigs_ori');
end
%close;

dLam_Del=eig(A+B*K_Del);
figure; hold on;
t_circle = linspace(0,2*pi,100);
plot(sin(t_circle),cos(t_circle),'k--');
xline(0); yline(0);
%scatter(real(dLam_Del),imag(dLam_Del),p.msize,'x','MarkerEdgeColor',[0.4940 0.1840 0.5560],'LineWidth',p.mline_width);
scatter(real(dLam_Del),imag(dLam_Del),p.msize,'x','MarkerEdgeColor',[0,0,0],'LineWidth',p.mline_width);
xlim([-circle_lim,circle_lim]); ylim([-circle_lim,circle_lim]); axis square;
xlabel('Re$(\lambda)$','FontSize',p.fs,'interpreter','latex');
ylabel('Im$(\lambda)$','FontSize',p.fs,'interpreter','latex');
set(gca,'Fontsize',p.label_size);
if p.save
    p.pl_save.save_png('eigs_pert');
end
%close;

dLam_Del_rand=eig(A+B*K_Del_rand);
figure; hold on;
t_circle = linspace(0,2*pi,100);
plot(sin(t_circle),cos(t_circle),'k--');
xline(0); yline(0);
%scatter(real(dLam_Del_rand),imag(dLam_Del_rand),p.msize,'x','MarkerEdgeColor',[0.4940 0.1840 0.5560],'LineWidth',p.mline_width);
scatter(real(dLam_Del_rand),imag(dLam_Del_rand),p.msize,'x','MarkerEdgeColor',[0,0,0],'LineWidth',p.mline_width);
xlim([-circle_lim,circle_lim]); ylim([-circle_lim,circle_lim]); axis square;
xlabel('Re$(\lambda)$','FontSize',p.fs,'interpreter','latex');
ylabel('Im$(\lambda)$','FontSize',p.fs,'interpreter','latex');
set(gca,'Fontsize',p.label_size);
if p.save
    %p.pl_save.save_png('eigs_pert_rand');
end
close;

cond_nums=cond_nums(cond_nums~=0);
mean_demo=mean(cond_nums);
std_demo=std(cond_nums);

%open_sys=ss(A,eye(n),eye(n),0,Ts);
%cl_sys=ss(A+B*K_dd,eye(n),eye(n),0,Ts);
%cl_sys_per=ss(A+B*K_Del,eye(n),eye(n),0,Ts);