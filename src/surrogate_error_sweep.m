function res = surrogate_error_sweep(cfg, C, f, L, Phi, Xi, gamma, Theta_te, G_te, sigma_obs)
ntr   = cfg.ntr_list;
nB    = numel(ntr);
nR    = numel(cfg.r_list);
n_max = max(ntr);

res.alone = zeros(cfg.n_rep, nB);
res.lis   = zeros(cfg.n_rep, nB, nR);
res.floor = zeros(1, nR);

Xte_alone = Theta_te - cfg.mu_pr;
for a = 1:nR
    r = cfg.r_list(a);
    res.floor(a) = truncation_error(Phi(:,1:r), Xi(:,1:r), Theta_te, G_te, cfg, C, f) / sigma_obs;
end

for rep = 1:cfg.n_rep
    Theta = sample_prior(n_max, cfg, L, cfg.design_seed + rep);
    Y = zeros(cfg.m, n_max);
    for i = 1:n_max
        Y(:,i) = fem_solve(Theta(:,i), cfg, C, f);
    end

    Yl = cell(1, nR);
    for a = 1:nR
        r  = cfg.r_list(a);
        Xr = Xi(:,1:r)' * (Theta - cfg.mu_pr);
        Yl{a} = zeros(cfg.m, n_max);
        for i = 1:n_max
            Yl{a}(:,i) = fem_solve(cfg.mu_pr + Phi(:,1:r)*Xr(:,i), cfg, C, f);
        end
    end

    for b = 1:nB
        nt = ntr(b);
        gp = gp_fit(Theta(:,1:nt) - cfg.mu_pr, Y(:,1:nt), cfg.gp_jitter);
        res.alone(rep,b) = rms_error(gp_predict(gp, Xte_alone), G_te) / sigma_obs;

        for a = 1:nR
            r  = cfg.r_list(a);
            sc = sqrt(gamma(1:r));
            Xtr = sc .* (Xi(:,1:r)' * (Theta(:,1:nt) - cfg.mu_pr));
            Xte = sc .* (Xi(:,1:r)' * (Theta_te - cfg.mu_pr));
            gp  = gp_fit(Xtr, Yl{a}(:,1:nt), cfg.gp_jitter);
            res.lis(rep,b,a) = rms_error(gp_predict(gp, Xte), G_te) / sigma_obs;
        end
    end
    fprintf('  design %d of %d done\n', rep, cfg.n_rep);
end
end

function e = rms_error(G_hat, G)
e = sqrt(mean(sum((G_hat - G).^2, 1)));
end

function e = truncation_error(Phi_r, Xi_r, Theta_te, G_te, cfg, C, f)
Tr = cfg.mu_pr + Phi_r * (Xi_r' * (Theta_te - cfg.mu_pr));
G  = zeros(size(G_te));
for k = 1:size(Tr, 2)
    G(:,k) = fem_solve(Tr(:,k), cfg, C, f);
end
e = rms_error(G, G_te);
end
