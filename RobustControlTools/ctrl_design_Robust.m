function [K,P,S,L,Lams,J] = ctrl_design_Robust(D,n,m,Q,R,eps,D_norms)
%CTRL_DESIGN generate the robust controller

T=size(D,2);
Z=D(1:n,:);
X=D(n+1:2*n,:);
U=D(2*n+1:end,:);
Rs=sqrtm(R);
Dbar=[eye(2*n+m) [Z; -X; -U]]';
eps=eps*max(D_norms);
ebar=eps^2*(2*n+m)*T;
N=Dbar'*[ebar*eye(2*n+m) zeros(2*n+m,T);
    zeros(T,2*n+m) -eye(T)]*Dbar;

cvx_begin sdp quiet
    variable S(m,m) semidefinite
    variable L(T,n)
    variable P(n,n) semidefinite
    variables alp bet
    dual variable Lam1
    dual variable Lam2
    dual variable Lam3
    dual variable Lam4
    dual variable Lam5
    minimize trace(Q*X*L) + trace(S)
    subject to
        [S Rs*U*L;
            (Rs*U*L)' P] >= 0 : Lam1;
        [P-eye(n) Z*L;
            (Z*L)' P] >= 0 : Lam2;
        P-X*L == 0;
        [P-bet*eye(n) zeros(n) zeros(n,m) zeros(n);
            zeros(n) -P -L'*U' zeros(n);
            zeros(m,n) -U*L zeros(m) U*L;
            zeros(n) zeros(n) L'*U' P] - alp*blkdiag(N,zeros(n)) >= 0 : Lam3;
        alp>=0 : Lam4;
        bet>=1 : Lam5;
cvx_end

%P=X*L;
%K=U*L/P;
K=K_func_LQR(D,L);
Lams=cell(1,5);
Lams{1}=Lam1; Lams{2}=Lam2; Lams{3}=Lam3; Lams{4}=Lam4; Lams{5}=Lam5;
J=trace(Q*P)+trace(S);

end

