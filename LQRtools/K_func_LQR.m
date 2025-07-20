function K = K_func_LQR(D,L)
n=size(L,2);
X=D(n+1:2*n,:);
U=D(2*n+1:end,:);
K=U*L/(X*L);
end

