% objective in whitened coordinaates 
function [F, Jw, res] = map_objective(v, cfg, C, f, L, y_obs, sigma_obs)

theta = cfg.mu_pr + L*v;
y     = fem_solve(theta, cfg, C, f);
res   = (y - y_obs) / sigma_obs;
F     = 0.5*(res'*res) + 0.5*(v'*v);

if nargout > 1
    J  = adjoint_jacobian(theta, cfg, C, f);
    Jw = (J*L) / sigma_obs;
end
end
