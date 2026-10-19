function [Phi, Xi, gamma, info] = adaptive_lis(cfg, C, f, L, y_obs, sigma_obs, v_map, H_map)
n = cfg.n;

[lambda, V] = local_lis(cfg.mu_pr + L*v_map, cfg, C, f, L, sigma_obs);
S = (V .* lambda') * V';
N = 1;
[Psi, gamma] = lis_eig(S / N);

r = max(1, sum(gamma >= cfg.tau_g));
x = Psi(:,1:r)' * v_map;

rng(cfg.lis_seed, 'twister');
info.dist = zeros(1, cfg.lis_iters);
info.acc  = zeros(1, cfg.lis_iters);
info.r    = zeros(1, cfg.lis_iters);
info.Theta = zeros(n, cfg.lis_iters + 1);
info.Theta(:,1) = cfg.mu_pr + L*v_map;

for k = 1:cfg.lis_iters
    Psi_r = Psi(:,1:r);
    Phi_r = L * Psi_r;

    H_red = Psi_r' * H_map * Psi_r;
    R     = chol(inv((H_red + H_red') / 2));
    logpi = @(z) -0.5*sum(((fem_solve(cfg.mu_pr + Phi_r*z, cfg, C, f) - y_obs)/sigma_obs).^2) ...
                 - 0.5*(z'*z);

    [x, info.acc(k)] = subchain(logpi, x, R, 2.38/sqrt(r), cfg.lis_L);

    xi    = randn(n, 1);
    xi    = xi - Psi_r*(Psi_r'*xi);
    theta = cfg.mu_pr + Phi_r*x + L*xi;
    info.Theta(:,k+1) = theta;

    [lambda, V] = local_lis(theta, cfg, C, f, L, sigma_obs);
    S = S + (V .* lambda') * V';
    N = N + 1;

    [Psi_new, gamma_new] = lis_eig(S / N);
    info.dist(k) = weighted_subspace_distance(Psi, gamma, Psi_new, gamma_new, cfg.tau_g);

    r_new = max(1, sum(gamma_new >= cfg.tau_g));
    x     = Psi_new(:,1:r_new)' * (Psi_r * x);
    Psi   = Psi_new;
    gamma = gamma_new;
    r     = r_new;
    info.r(k) = r;
end

Phi = L  * Psi;
Xi  = L' \ Psi;
info.N_hessians = N;
end

function [x, acc] = subchain(logpi, x, R, s, n_steps)
d  = numel(x);
lp = logpi(x);
n_acc = 0;
for k = 1:n_steps
    xp  = x + s * (R' * randn(d,1));
    lpp = logpi(xp);
    if log(rand) < lpp - lp
        x  = xp;
        lp = lpp;
        n_acc = n_acc + 1;
    end
end
acc = n_acc / n_steps;
end
