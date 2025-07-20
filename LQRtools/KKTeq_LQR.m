function out = KKTeq_LQR(S,L,Lam1,Lam2,D,Q,R,gamma,n)

T=size(D,2);
Z=D(1:n,:);
X=D(n+1:2*n,:);
U=D(2*n+1:end,:);
W0=[X;U];
Rs=sqrtm(R);
Pi=eye(T)-pinv(W0)*W0;

F1=F1_func(S,L,Rs,U,X);
F2=F2_func(L,Z,n,X);

LagS=@(S) Lagrangian_LQR(S,L,Lam1,Lam2,Z,X,U,Q,Rs,gamma,Pi,n);
LagL=@(L) Lagrangian_LQR(S,L,Lam1,Lam2,Z,X,U,Q,Rs,gamma,Pi,n);
drLagdrS=AutoDiff(LagS,S);
drLagdrL=AutoDiff(LagL,L);
H3=F1*Lam1; H4=F2*Lam2;
%H3=trace(F1*Lam1); H4=trace(F2*Lam2);
H5=Lam1-Lam1'; H6=Lam2-Lam2';

out=[drLagdrS(:);
    drLagdrL(:);
    H3(:);
    H4(:);
    H5(:);
    H6(:);
    ];
end

function Lag=Lagrangian_LQR(S,L,Lam1,Lam2,Z,X,U,Q,Rs,gamma,Pi,n)
J=trace(Q*X*L)+trace(S)+gamma*norm(Pi*L);
F1=F1_func(S,L,Rs,U,X);
F2=F2_func(L,Z,n,X);
Lag=J-trace(F1*Lam1)-trace(F2*Lam2);
end

function out = F1_func(S,L,Rs,U,X)
out=[S Rs*U*L;
    (Rs*U*L)' X*L];
end

function out = F2_func(L,Z,n,X)
out=[X*L-eye(n) Z*L;
    (Z*L)' X*L];
end