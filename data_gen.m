
% Signal Generation
[U,X,Z,~,~]=gen_sig_exp_rep(n,m,T,w,v,A,B);
D=[Z;X;U];
if strcmp(pert_norm,'Fro')
    Z_nor=norm(Z,'fro'); X_nor=norm(X,'fro'); U_nor=norm(U,'fro');
elseif strcmp(pert_norm,'max')
    Z_nor=max(abs(Z),[],"all"); X_nor=max(abs(X),[],"all"); U_nor=max(abs(U),[],"all");
end
D_norms=[Z_nor,X_nor,U_nor];