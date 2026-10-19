function post = posterior_sweep(cfg, C, f, L, Phi, Xi, gamma, v_map, H_map, y_obs, sigma_obs, r)
ntr   = cfg.ntr_list;
nB    = numel(ntr);
n_max = max(ntr);
n     = cfg.n;

Phi_r = Phi(:, 1:r);
Psi_r = L' * Xi(:, 1:r);
sc    = sqrt(gamma(1:r));

Hs     = (H_map + H_map') / 2;
R_full = chol(inv(Hs));
Hr     = Psi_r' * Hs * Psi_r;
R_red  = chol(inv((Hr + Hr') / 2));

misfit = @(G) 0.5*sum(((G - y_obs)/sigma_obs).^2);

logpi = @(v) -misfit(fem_solve(cfg.mu_pr + L*v, cfg, C, f)) - 0.5*(v'*v);
V = run_chain(logpi, v_map, R_full, cfg, cfg.mcmc_seed + 1);
post.full.m = mean(V, 2);
post.full.S = cov(V');

logpi = @(x) -misfit(fem_solve(cfg.mu_pr + Phi_r*x, cfg, C, f)) - 0.5*(x'*x);
X = run_chain(logpi, Psi_r' * v_map, R_red, cfg, cfg.mcmc_seed + 2);
[post.lis.m, post.lis.S] = lis_moments(X, Psi_r, n);

post.r = r;
post.lisgp = cell(cfg.n_rep, nB);
post.alone = cell(cfg.n_rep, nB);

for rep = 1:cfg.n_rep
    Theta = sample_prior(n_max, cfg, L, cfg.design_seed + rep);
    Xr = Xi(:, 1:r)' * (Theta - cfg.mu_pr);
    Y  = zeros(cfg.m, n_max);
    Yr = zeros(cfg.m, n_max);
    for i = 1:n_max
        Y(:,i)  = fem_solve(Theta(:,i), cfg, C, f);
        Yr(:,i) = fem_solve(cfg.mu_pr + Phi_r*Xr(:,i), cfg, C, f);
    end

    for b = 1:nB
        nt   = ntr(b);
        seed = cfg.mcmc_seed + 100*rep + b;

        gp = gp_fit(sc .* Xr(:,1:nt), Yr(:,1:nt), cfg.gp_jitter);
        logpi = @(x) -misfit(gp_predict(gp, sc .* x)) - 0.5*(x'*x);
        X = run_chain(logpi, Psi_r' * v_map, R_red, cfg, seed);
        [post.lisgp{rep,b}.m, post.lisgp{rep,b}.S] = lis_moments(X, Psi_r, n);

        gp = gp_fit(Theta(:,1:nt) - cfg.mu_pr, Y(:,1:nt), cfg.gp_jitter);
        logpi = @(v) -misfit(gp_predict(gp, L*v)) - 0.5*(v'*v);
        V = run_chain(logpi, v_map, R_full, cfg, seed);
        post.alone{rep,b}.m = mean(V, 2);
        post.alone{rep,b}.S = cov(V');
    end
    fprintf('  design %d of %d done\n', rep, cfg.n_rep);
end
end

function X = run_chain(logpi, x0, R, cfg, seed)
opts = struct('s0', 2.38/sqrt(numel(x0)), 'n_pilot', cfg.mcmc_pilot, ...
              'acc_target', cfg.mcmc_acc_target, 'seed', seed);
X = mcmc_rw(logpi, x0, R, cfg.chain_steps, opts);
X = X(:, round(cfg.mcmc_burnin * cfg.chain_steps) + 1 : end);
end

function [m, S] = lis_moments(X, Psi_r, n)
r = size(Psi_r, 2);
m = Psi_r * mean(X, 2);
S = eye(n) + Psi_r * (cov(X') - eye(r)) * Psi_r';
end
