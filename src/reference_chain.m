function Vf = reference_chain(cfg, C, f, L, y_obs, sigma_obs, v_map, H_map)
cache = fullfile(pwd, cfg.ref_cache);
if isfile(cache)
    S = load(cache);
    if S.steps == cfg.mcmc_steps_full && S.seed == cfg.mcmc_seed && S.thin == cfg.ref_thin
        Vf = S.Vf;
        fprintf('reference chain loaded from %s\n', cache);
        return;
    end
end

logpi = @(v) -0.5*sum(((fem_solve(cfg.mu_pr + L*v, cfg, C, f) - y_obs)/sigma_obs).^2) - 0.5*(v'*v);
R_full = chol(inv((H_map + H_map')/2));
opts = struct('s0', 2.38/sqrt(cfg.n), 'n_pilot', cfg.mcmc_pilot, ...
              'acc_target', cfg.mcmc_acc_target, 'seed', cfg.mcmc_seed);
tic;
[V, acc] = mcmc_rw(logpi, v_map, R_full, cfg.mcmc_steps_full, opts);
fprintf('reference chain: %d steps, %.1f s, acceptance %.3f\n', cfg.mcmc_steps_full, toc, acc);

burn = round(cfg.mcmc_burnin * cfg.mcmc_steps_full);
Vf   = V(:, burn+1:cfg.ref_thin:end);
steps = cfg.mcmc_steps_full;
seed  = cfg.mcmc_seed;
thin  = cfg.ref_thin;
save(cache, 'Vf', 'steps', 'seed', 'thin');
end
