function [Z,X,U] = sep_ZXU(D,n,m)
Z=D(1:n,:);
X=D(n+1:2*n,:);
U=D(2*n+1:end,:);
end

