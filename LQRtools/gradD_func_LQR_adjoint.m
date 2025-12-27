function gradD = gradD_func_LQR_adjoint(K,S,L,Lam1,Lam2,D,Q,R,gamma,n,m,T,rho,A,B,anal_mats)

V=sqrtm(R);
E_Z=anal_mats.E_Z;
E_X=anal_mats.E_X;
E_U=anal_mats.E_U;
Z=E_Z*D; X=E_X*D; U=E_U*D;
Gam=[U;X];
Pi=eye(T)-pinv(Gam)*Gam;
Lam=blkdiag(Lam1,Lam2);
F=blkdiag([S V*E_U*D*L; (V*E_U*D*L)' E_X*D*L],[E_X*D*L-eye(n) E_Z*D*L; (E_Z*D*L)' E_X*D*L]);

E_1=anal_mats.E_1;
E_2=anal_mats.E_2;
kE_X=anal_mats.kE_X;
kE_U=anal_mats.kE_U;

%% Proposition3
drhodK=drhodK_func(A,B,K);
tmp=transpose(inv(X*L));
FK=E_U*D*L/(E_X*D*L);
dFKdL=kron(tmp,U-FK*X);
tmp=tmp*transpose(L);
dFKdD=kron(tmp,E_U-FK*E_X);

%% Ajoint Method
WT=drhodK*dFKdL;
WTdLdD=py.interface.WTdLdD(L, S, Lam, D, Q, R, gamma, n, m, T, V, Pi, WT);
gradD=WTdLdD+drhodK*dFKdD;

end


function drdK = drhodK_func(A,B,K)
[V,D,W] = eig(A+B*K);
[~,idx] = max(abs(diag(D)));
lam = D(idx, idx);
v = V(:, idx);
v = v / norm(v);
w = W(:, idx);
w = w / norm(w);
tmp=w.'*(kron(v.',B))/(w.'*v);
tmp=lam*tmp/abs(lam);
drdK=real(tmp);
end
