function [K,P,S,L,Lam1,Lam2,J] = ctrl_design_LQR(D,n,m,Q,R,gamma)
%CTRL_DESIGN generate the optimal LQR controller

T=size(D,2);
Z=D(1:n,:);
X=D(n+1:2*n,:);
U=D(2*n+1:end,:);
W0=[X;U];
Rs=sqrtm(R);
Pi=eye(T)-pinv(W0)*W0;

cvx_begin sdp quiet
    variable S(m,m) semidefinite
    variable L(T,n)
    variable P(n,n) semidefinite
    dual variable Lam1
    dual variable Lam2
    %minimize trace(Q*X*L) + trace(S) + gamma*norm(Pi*L)
    minimize trace(Q*X*L) + trace(S) + gamma*sum_square(vec(Pi*L))
    subject to
        [S Rs*U*L;
            (Rs*U*L)' P] >= 0 : Lam1;
        [P-eye(n) Z*L;
            (Z*L)' P] >= 0 : Lam2;
        P-X*L == 0;
cvx_end

%P=X*L;
%K=U*L/P;
K=K_func_LQR(D,L);
J=trace(Q*P)+trace(S);

end

