function gradD = gradD_func_LQR(K_Del,S_Del,L_Del,Lam1_Del,Lam2_Del,D_Del,Q,R,gamma,n,m,T,rho)

% dLdD by implicit differentiation
HS=@(S) KKTeq_LQR(S,L_Del,Lam1_Del,Lam2_Del,D_Del,Q,R,gamma,n);
HL=@(L) KKTeq_LQR(S_Del,L,Lam1_Del,Lam2_Del,D_Del,Q,R,gamma,n);
HLam1=@(Lam1) KKTeq_LQR(S_Del,L_Del,Lam1,Lam2_Del,D_Del,Q,R,gamma,n);
HLam2=@(Lam2) KKTeq_LQR(S_Del,L_Del,Lam1_Del,Lam2,D_Del,Q,R,gamma,n);
HD=@(D) KKTeq_LQR(S_Del,L_Del,Lam1_Del,Lam2_Del,D,Q,R,gamma,n);
dHdS=AutoDiff(HS,S_Del);
dHdL=AutoDiff(HL,L_Del);
dHdLam1=AutoDiff(HLam1,Lam1_Del);
dHdLam2=AutoDiff(HLam2,Lam2_Del);
dHdD=AutoDiff(HD,D_Del);
Amat=[dHdS dHdL dHdLam1 dHdLam2]; Bmat=-dHdD;
LinSol=linsolve(Amat,Bmat);
dLdD=LinSol(m^2+1:m^2+n*T,:);

% others

Eu=[zeros(m,2*n) eye(m)]; Ex=[zeros(n) eye(n) zeros(n,m)];
Kd=Ex*D_Del*L_Del;
tmp1=inv(Kd)'; tmp2=Eu-K_Del*Ex;
drKdrD=kron(tmp1*L_Del',tmp2);
drKdrL=kron(tmp1,tmp2*D_Del);


%{
K_func_LQR_D=@(D) K_func_LQR(D,L_Del);
K_func_LQR_L=@(L) K_func_LQR(D_Del,L);
drKdrD=AutoDiff(K_func_LQR_D,D_Del);
drKdrL=AutoDiff(K_func_LQR_L,L_Del);
%}

drhodK=AutoDiff(rho,K_Del);

dKdD=drKdrD+drKdrL*dLdD;
drhodD=drhodK*dKdD;

gradZ=reshape(drhodD(1:n*T),[n,T]);
gradX=reshape(drhodD(n*T+1:2*n*T),[n,T]);
gradU=reshape(drhodD(2*n*T+1:end),[m,T]);
gradD_old=[gradZ;gradX;gradU];

gradD=reshape(drhodD,[2*n+m,T]);

end

