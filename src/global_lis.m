function [Phi, Xi, gamma] = global_lis(Theta, cfg, C, f, L, sigma_obs)
N = size(Theta, 2);
S = zeros(cfg.n, cfg.n);
for k = 1:N
    [lambda, V] = local_lis(Theta(:,k), cfg, C, f, L, sigma_obs);
    S = S + (V .* lambda') * V';
end
[Psi, gamma] = lis_eig(S / N);
Phi = L  * Psi;
Xi  = L' \ Psi;
end
