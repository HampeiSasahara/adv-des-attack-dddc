function [gradD,dGdvar,dGdD,LinSol] = gradD_func_LQR_anal(K,S,L,Lam1,Lam2,D,Q,R,gamma,n,m,T,rho,A,B,anal_mats)

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

%% Proposition4
tmp1=[zeros(m,n); eye(n)]; tmp2=[V*U;zeros(n,T)];
C=Create_ComMat(T,n);
dF1dL=kron(tmp1,tmp2)+(kron(tmp2,tmp1))*C+kron(tmp1,[zeros(m,T);X]);
tmp1=[eye(n);zeros(n,n)]; tmp2=[zeros(n,n);eye(n)]; tmp3=[Z; zeros(n,T)];
dF2dL=kron(tmp1,[X;zeros(n,T)])+kron(tmp2,tmp3)+(kron(tmp3,tmp2))*C+kron(tmp2,[zeros(n,T);X]);
dF1dS=anal_mats.dF1dS;
dF2dS=anal_mats.dF2dS;
dFdL=E_1*dF1dL+E_2*dF2dL;
dFdS=E_1*dF1dS+E_2*dF2dS;
%tmp1=(Q*X)'; tmp1=tmp1(:);
%tmp2=Pi*L; tmp2=tmp2(:);
%tmp3=Lam(:);
%dLagdL=tmp1'+2*gamma*tmp2'-tmp3'*dFdL;
%tmp1=eye(m); tmp1=tmp1(:);
%dLagdS=tmp1'-tmp3'*dFdS;

%% Proposition5
nb=2*n+m;
tmp1=[zeros(m,T);L']; tmp2=[V*E_U;zeros(n,nb)]; tmp3=[zeros(m,nb);E_X];
C=Create_ComMat(nb,T);
dF1dD=kron(tmp1,tmp2)+kron(tmp2,tmp1)*C+kron(tmp1,tmp3);
tmp1=[L';zeros(n,T)]; tmp2=[E_X;zeros(n,nb)]; tmp3=[zeros(n,T);L']; tmp4=[E_Z;zeros(n,nb)]; tmp5=[zeros(n,nb);E_X];
dF2dD=kron(tmp1,tmp2)+kron(tmp3,tmp4)+(kron(tmp4,tmp3))*C+kron(tmp3,tmp5);

ddF1dDdL=anal_mats.ddF1dDdL;
ddF2dDdL=anal_mats.ddF2dDdL;

dFdD=E_1*dF1dD+E_2*dF2dD;
ddFdDdL=anal_mats.ddFdDdL_1*ddF1dDdL+anal_mats.ddFdDdL_2*ddF2dDdL;


% Gs
Gam=[eye(m); zeros(n,m)]*U+[zeros(m,n); eye(n)]*X;
C=Create_ComMat(m+n,T);
dGam4dGam=kron(Gam,eye(m+n))+kron(eye(m+n),Gam)*C;
Gam4=Gam*Gam';
iGam4=inv(Gam4);
dGam2dGam=-kron(transpose(iGam4),iGam4)*dGam4dGam;
dGamdagdGam=kron(iGam4',eye(T))*C+kron(eye(m+n),Gam')*dGam2dGam;

Gamdag=Gam'*iGam4;
dGamdagdX=dGamdagdGam*anal_mats.dGamdX;
dGamdagdU=dGamdagdGam*anal_mats.dGamdU;
tmp1=-kron(Gam',eye(T));
tmp2=-kron(eye(T),Gamdag);
dPidX=tmp1*dGamdagdX+tmp2*anal_mats.dGamdX;
dPidU=tmp1*dGamdagdU+tmp2*anal_mats.dGamdU;
dPidD=dPidX*kE_X+dPidU*kE_U;
dPiLdD=kron(L',eye(T))*dPidD;


tmp=Lam(:); tmp=tmp';
dG1dL=2*gamma*kron(eye(n),Pi);
dG1dS=zeros(n*T,m^2);
dG1dLam=-dFdL';
C=Create_ComMat(n,T);
dG1dD=C*kron(eye(T),Q*E_X)+2*gamma*dPiLdD-kron(eye(n*T),tmp)*ddFdDdL;
dG2dL=zeros(m^2,n*T);
dG2dS=zeros(m^2,m^2);
dG2dLam=-dFdS';
%dG2dD=-kron(eye(m^2),tmp)*ddFdDdS;
dG2dD=zeros(m^2,nb*T);
dG3dL=kron(Lam,eye(3*n+m))*dFdL;
dG3dS=kron(Lam,eye(3*n+m))*dFdS;
C=Create_ComMat(3*n+m,3*n+m);
dG3dLam=kron(eye(3*n+m),F)*C;
dG3dD=kron(Lam,eye(3*n+m))*dFdD;
dG4dL=zeros((3*n+m)^2,n*T);
dG4dS=zeros((3*n+m)^2,m^2);
dG4dLam=eye((3*n+m)^2)-C;
dG4dD=zeros((3*n+m)^2,nb*T);

nbb=2*n*(m+n);
dG5dL=zeros(nbb,n*T);
dG5dS=zeros(nbb,m^2);
dG5dLam=kron([zeros(2*n,m+n) eye(2*n)],[eye(m+n) zeros(m+n,2*n)]);
dG5dD=zeros(nbb,nb*T);
dG6dL=zeros(nbb,n*T);
dG6dS=zeros(nbb,m^2);
dG6dLam=kron([eye(m+n) zeros(m+n,2*n)],[zeros(2*n,m+n) eye(2*n)]);
dG6dD=zeros(nbb,nb*T);

dGdvar=[dG1dL dG1dS dG1dLam;
    dG2dL dG2dS dG2dLam;
    dG3dL dG3dS dG3dLam;
    dG4dL dG4dS dG4dLam;
    dG5dL dG5dS dG5dLam;
    dG6dL dG6dS dG6dLam];
dGdD=[dG1dD; dG2dD; dG3dD; dG4dD; dG5dD; dG6dD];
LinSol=linsolve(dGdvar,-dGdD);
dLdD=LinSol(1:n*T,:);
dSdD=LinSol(n*T+1:n*T+m^2,:);
dLamdD=LinSol(n*T+m^2+1:end,:);

drhodD=drhodK*(dFKdL*dLdD+dFKdD);
gradD=reshape(drhodD,[nb,T]);

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

function K = Create_ComMat(p,q)
I = reshape(1:p*q, [p, q]);
I = I';
I = I(:);
Y = eye(p*q);
K = Y(I,:);
end
