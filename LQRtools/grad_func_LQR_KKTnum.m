function [gradD,dGdvar,dGdD,LinSol] = grad_func_LQR_KKTnum(K,S,L,Lam1,Lam2,D,Q,R,gamma,n,m,T,rho,A,B,anal_mats)

V=sqrtm(R);
E_Z=anal_mats.E_Z;
E_X=anal_mats.E_X;
E_U=anal_mats.E_U;
Z=E_Z*D; X=E_X*D; U=E_U*D;
Gam=[U;X];
Pi=eye(T)-pinv(Gam)*Gam;
Lam=blkdiag(Lam1,Lam2);
F=blkdiag([S V*E_U*D*L; (V*E_U*D*L)' E_X*D*L],[E_X*D*L-eye(n) E_Z*D*L; (E_Z*D*L)' E_X*D*L]);

G=@(L,S,Lam,D) KKTeq_LQR2(K,S,L,Lam,D,Q,R,gamma,n,m,T,rho,A,B);
GL=@(L) G(L,S,Lam,D);
GS=@(S) G(L,S,Lam,D);
GLam=@(Lam) G(L,S,Lam,D);
GD=@(D) G(L,S,Lam,D);

dGdL=AutoDiff(GL,L);
dGdS=AutoDiff(GS,S);
dGdLam=AutoDiff(GLam,Lam);
dGdD=AutoDiff(GD,D);

dGdvar=[dGdL dGdS dGdLam];
LinSol=linsolve(dGdvar,-dGdD);
dLdD=LinSol(1:n*T,:);
dSdD=LinSol(n*T+1:n*T+m^2,:);
dLamdD=LinSol(n*T+m^2+1:end,:);

K_func_LQR_D=@(D) K_func_LQR(D,L);
K_func_LQR_L=@(L) K_func_LQR(D,L);
drKdrD=AutoDiff(K_func_LQR_D,D);
drKdrL=AutoDiff(K_func_LQR_L,L);
drhodK=AutoDiff(rho,K);
dKdD=drKdrD+drKdrL*dLdD;
drhodD=drhodK*dKdD;
gradD=reshape(drhodD,[2*n+m,T]);

end

