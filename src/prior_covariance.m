% Prior covariance and its cholesky factor L 
function [Gamma_pr, L] = prior_covariance(cfg)

x = cfg.x_mid;
D = abs(x - x');

Gamma_pr = cfg.sigma_theta^2 * exp(-D / cfg.ell_pr);
Gamma_pr = (Gamma_pr + Gamma_pr') / 2;

L = chol(Gamma_pr, 'lower');
end
