clear; clc; close all;
addpath(genpath('src'));

cfg = config();
f             = load_vector(cfg);
[C, ~]        = sensor_matrix(cfg);
[~, L]        = prior_covariance(cfg);
[theta_true, sigma_obs] = ground_truth(cfg, C, f);
y_obs = generate_synthetic_data(cfg, C, f, theta_true, sigma_obs);
[~, v_map] = compute_map(cfg, C, f, L, y_obs, sigma_obs);
[~, Jw_map] = map_objective(v_map, cfg, C, f, L, y_obs, sigma_obs);
H_map = Jw_map'*Jw_map + eye(cfg.n);

%% 1. Jacobian against central finite differences
theta = cfg.mu_pr + L*randn(cfg.n, 1);
J  = adjoint_jacobian(theta, cfg, C, f);
Jn = zeros(cfg.m, cfg.n);
h  = 1e-6;
for e = 1:cfg.n
    d = zeros(cfg.n, 1); d(e) = h;
    Jn(:,e) = (fem_solve(theta + d, cfg, C, f) - fem_solve(theta - d, cfg, C, f)) / (2*h);
end
fprintf('1. Jacobian relative error: %.2e\n', norm(J - Jn, 'fro') / norm(Jn, 'fro'));

%% 2. Projector algebra
Theta = sample_prior(200, cfg, L, cfg.sample_seed);
[Phi_t, Xi_t, gamma_t] = global_lis(Theta, cfg, C, f, L, sigma_obs);
r = 8;
Pi_r = Phi_t(:,1:r) * Xi_t(:,1:r)';
fprintf('2. ||Xi''*Phi - I|| = %.2e, ||Pi^2 - Pi|| = %.2e\n', ...
        norm(Xi_t'*Phi_t - eye(cfg.n), 'fro'), norm(Pi_r*Pi_r - Pi_r, 'fro'));

%% 3. Algorithm 1 against a reference built from the full-space chain
logpi_full = @(v) -0.5*sum(((fem_solve(cfg.mu_pr + L*v, cfg, C, f) - y_obs)/sigma_obs).^2) - 0.5*(v'*v);
R_full = chol(inv((H_map + H_map')/2));
opts   = struct('s0', 2.38/sqrt(cfg.n), 'n_pilot', cfg.mcmc_pilot, ...
                'acc_target', cfg.mcmc_acc_target, 'seed', cfg.mcmc_seed);
tic;
Vf = mcmc_rw(logpi_full, v_map, R_full, cfg.mcmc_steps_full, opts);
Vf = Vf(:, round(cfg.mcmc_burnin*cfg.mcmc_steps_full)+1:end);
fprintf('3. full-space chain: %.1f s\n', toc);

idx_ref   = round(linspace(1, size(Vf,2), 1000));
Theta_ref = cfg.mu_pr + L*Vf(:, idx_ref);
[Phi_ref, Xi_ref, gamma_ref] = global_lis(Theta_ref, cfg, C, f, L, sigma_obs);

[Phi_ada, Xi_ada, gamma_ada, info] = adaptive_lis(cfg, C, f, L, y_obs, sigma_obs, v_map, H_map);
fprintf('   Algorithm 1: %d Hessians, mean acceptance %.3f\n', info.N_hessians, mean(info.acc));
fprintf('   directions with gamma >= %.2f: reference %d, Algorithm 1 %d\n', ...
        cfg.tau_g, sum(gamma_ref >= cfg.tau_g), sum(gamma_ada >= cfg.tau_g));
fprintf('   spectrum reference   : %s\n', sprintf('%.3e ', gamma_ref(1:10)));
fprintf('   spectrum Algorithm 1 : %s\n', sprintf('%.3e ', gamma_ada(1:10)));

idx_te = round(linspace(1, size(Vf,2), 200));
Theta_te = cfg.mu_pr + L*Vf(:, idx_te);
G_te = zeros(cfg.m, 200);
for k = 1:200
    G_te(:,k) = fem_solve(Theta_te(:,k), cfg, C, f);
end

fprintf('\n   %3s %14s %14s\n', 'r', 'floor Alg. 1', 'floor ref');
for r = [4 6 8]
    fl_a  = truncation_error(Phi_ada, Xi_ada, r, Theta_te, G_te, cfg, C, f);
    fl_r  = truncation_error(Phi_ref, Xi_ref, r, Theta_te, G_te, cfg, C, f);
    fprintf('   %3d %14.3e %14.3e\n', r, fl_a, fl_r);
end

function e = truncation_error(Phi, Xi, r, Theta_te, G_te, cfg, C, f)
Tr = cfg.mu_pr + Phi(:,1:r) * (Xi(:,1:r)' * (Theta_te - cfg.mu_pr));
Gr = zeros(size(G_te));
for k = 1:size(Tr,2)
    Gr(:,k) = fem_solve(Tr(:,k), cfg, C, f);
end
e = sqrt(mean(sum((Gr - G_te).^2, 1))) / sqrt(mean(sum(G_te.^2, 1)));
end
