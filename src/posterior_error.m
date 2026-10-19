function [e_mean, e_cov] = posterior_error(m, S, ref, L)
e_mean = norm(L*(m - ref.m)) / norm(L*ref.m);
e_cov  = norm(L*(S - ref.S)*L', 'fro') / norm(L*(eye(size(S,1)) - ref.S)*L', 'fro');
end
