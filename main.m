%%
%% Adversarial Destabilization Attacks to Direct Data-Driven Control
%%

tic
clear all;
addpath("tools","systems","scenarios","LQRtools","RobustControlTools","graphs","graphs/ASR");

%% Seed
rand_seed=3;
s=RandStream('mt19937ar','seed',rand_seed);
RandStream.setGlobalStream(s);
stream=RandStream.getGlobalStream;

%% Plot Parameters
p.label_size=20;
p.label_size2=16;
p.label_size_large=40;
p.fs=50;
p.fs2=30;
p.fs_l=20;
p.fs_large=66;
p.x_width=11;
p.y_width=8;
p.y_width2=5;
p.line_width=1;
p.line_width_thick=2;
p.mline_width=2;
p.msize=4000;
p.msize2=75;
p.pl_save=plot_save(p.x_width,p.y_width,p.y_width2);
p.save=false;

%% General Parameter Setting
cvx_precision default
pert_norm='max'; %'Fro': Frobenius, 'max': element-wise max
w=0; % disturbance size (variance: w^2)
v=0; % measurement noise size
q=1e0; r=1e-2; % LQR weights
%gamma=1e0; % regularization parameter
%gamma=1e-2;
gamma=1e-4;
Ite_num=50; % Max iteration number
T_str='medium'; % short, medium, long

%% Scenario
% 1. demonstration: instance for demonstration
% 2. deseps: destabilization ratio vs epsilon for particular system
% 3. grad_comp: gradient computation time comparison
% 4. defense_reg: regularization-based defense method
% 5. defense_rob: robust control-based defense method
% 6. transf: transferability across data 
scenario_flag=4;
switch scenario_flag
    case 1
        sys_flag="TT";
        build_sys;
        pert_size=5e-3; % perturbation size
        demonstration;
    case 2
        eps_exp_array=-7:1:0;
        num_sample=10;
        sys_flag="TT";
        build_sys;
        deseps;
    case 3
        sys_flag="TT";
        build_sys;
        grad_comp;
    case 4
        %gam_exp_array=-8:1:-4;
        gam_exp_array=[-5,-3];
        num_sample=20;
        pert_size=1e-2;
        sys_flag="TT";
        build_sys;
        defense_reg;
    case 5
        %pert_exp_array=-6:1:-2;
        %num_sample=2;
        eps_exp_array=-6:1:0;
        sys_flag="TT";
        build_sys;
        defense_rob;
    case 6
        eps_exp_array=-5:.5:0;
        num_sample=20;
        sys_flag="TT";
        build_sys;
        transf;
end

toc