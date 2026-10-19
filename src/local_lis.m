% Local LIS (problem 1 Cui) 
function [lambda, V] = local_lis(theta, cfg, C, f, L, sigma_obs)
J = adjoint_jacobian(theta, cfg, C, f);
M = (J * L) / sigma_obs;
[~, S, V] = svd(M, 'econ');
lambda = diag(S).^2;
keep   = lambda >= cfg.tau_loc;
lambda = lambda(keep);
V      = V(:, keep);
end
