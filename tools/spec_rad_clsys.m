function out = spec_rad_clsys(A,B,K)
out=abs(eigs(A+B*K,1));
end

