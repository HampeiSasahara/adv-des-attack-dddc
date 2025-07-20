function F = KKTeq_LQR2(K,S,L,Lam,D,Q,R,gamma,n,m,T,rho,A,B)

V=sqrtm(R);
E_Z=[eye(n) zeros(n,n) zeros(n,m)];
E_X=[zeros(n,n) eye(n) zeros(n,m)];
E_U=[zeros(m,n) zeros(m,n) eye(m)];
Z=E_Z*D; X=E_X*D; U=E_U*D;
Gam=[U;X];
Pi=eye(T)-pinv(Gam)*Gam;
F=blkdiag([S V*E_U*D*L; (V*E_U*D*L)' E_X*D*L],[E_X*D*L-eye(n) E_Z*D*L; (E_Z*D*L)' E_X*D*L]);

J_func=@(L,S,D) trace(Q*E_X*D*L)+trace(S)+gamma*norm(Pi*L,'fro')^2;
F1_func=@(L,S,D) [S V*E_U*D*L; (V*E_U*D*L)' E_X*D*L];
F2_func=@(L,S,D) [E_X*D*L-eye(n) E_Z*D*L; (E_Z*D*L)' E_X*D*L];
F_func=@(L,S,D) [F1_func(L,S,D) zeros(m+n,2*n); zeros(2*n,m+n) F2_func(L,S,D)];
Lag_func=@(L,S,Lam,D) J_func(L,S,D)-trace(F_func(L,S,D)*Lam');

Lag_func_L=@(L) Lag_func(L,S,Lam,D);
Lag_func_S=@(S) Lag_func(L,S,Lam,D);
dLagdL_num=AutoDiff(Lag_func_L,L);
dLagdS_num=AutoDiff(Lag_func_S,S);
G1=dLagdL_num;
G2=dLagdS_num;
G3=F*Lam';
G4=Lam-Lam';
G5=Lam(1:n+m,n+m+1:end);
G6=Lam(n+m+1:end,1:n+m);
%F=[G1(:);G2(:);G3(:);G4(:)];
F=[G1(:);G2(:);G3(:);G4(:);G5(:);G6(:)];

end

