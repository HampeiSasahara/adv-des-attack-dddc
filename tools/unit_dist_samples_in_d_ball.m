function s = unit_dist_samples_in_d_ball(d,n)
% https://jp.mathworks.com/matlabcentral/answers/439205-generate-n-random-uniformly-distributed-points-in-the-d-ball

s = randn(d,n);
r = rand(1,n).^(1/d);
c = r./sqrt(sum(s.^2,1));
s = bsxfun(@times, s, c);

end

