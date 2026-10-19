clear; clc; close all;
addpath(genpath('src'));

%% Setup
cfg = config();
f             = load_vector(cfg);
[C, sensors]  = sensor_matrix(cfg);
[Gamma_pr, L] = prior_covariance(cfg);
[theta_true, sigma_obs] = ground_truth(cfg, C, f);
y_obs = generate_synthetic_data(cfg, C, f, theta_true, sigma_obs);

[theta_map, v_map, map_hist] = compute_map(cfg, C, f, L, y_obs, sigma_obs);
[~, Jw_map] = map_objective(v_map, cfg, C, f, L, y_obs, sigma_obs);
H_map = Jw_map'*Jw_map + eye(cfg.n);

fprintf('n = %d, m = %d, sigma_obs = %.4e\n', cfg.n, cfg.m, sigma_obs);
fprintf('MAP: %d iterations, ||theta_map - theta_true|| = %.4e\n', ...
        map_hist.iters, norm(theta_map - theta_true));

%% Global LIS from prior samples
Theta_prior = sample_prior(cfg.N_prior, cfg, L, cfg.sample_seed);
[Phi_pri, Xi_pri, gamma_pri] = global_lis(Theta_prior, cfg, C, f, L, sigma_obs);
fprintf('prior samples: %d directions with gamma >= %.2f\n', sum(gamma_pri >= cfg.tau_g), cfg.tau_g);

%% Global LIS from posterior samples (Algorithm 1)
tic;
[Phi, Xi, gamma, lis_info] = adaptive_lis(cfg, C, f, L, y_obs, sigma_obs, v_map, H_map);
fprintf('Algorithm 1: %d iterations, %d Hessians, %.1f s\n', cfg.lis_iters, lis_info.N_hessians, toc);
fprintf('mean acceptance of the subchains: %.3f\n', mean(lis_info.acc));
fprintf('r(k): first %d, last %d\n', lis_info.r(1), lis_info.r(end));
fprintf('weighted subspace distance at k = 10, 50, 100, 200, 500: %s\n', ...
        sprintf('%.2e ', lis_info.dist([10 50 100 200 500])));
fprintf('posterior samples: %d directions with gamma >= %.2f\n', sum(gamma >= cfg.tau_g), cfg.tau_g);
fprintf('spectrum (leading 12): %s\n', sprintf('%.3e ', gamma(1:12)));

%% Spectra figure
n_show = 200;
lam_pri = zeros(cfg.m, n_show);
lam_pos = zeros(cfg.m, n_show);
Theta_pos = lis_info.Theta(:, end-n_show+1:end);
for k = 1:n_show
    lam_pri(:,k) = local_spectrum(Theta_prior(:,k), cfg, C, f, L, sigma_obs);
    lam_pos(:,k) = local_spectrum(Theta_pos(:,k),   cfg, C, f, L, sigma_obs);
end
fprintf('local directions with lambda >= %.2f (mean over %d samples): prior %.2f, posterior %.2f, ceiling m = %d\n', ...
        cfg.tau_loc, n_show, mean(sum(lam_pri >= cfg.tau_loc, 1)), mean(sum(lam_pos >= cfg.tau_loc, 1)), cfg.m);

figdir = fullfile(pwd, 'figures');
if ~isfolder(figdir), mkdir(figdir); end
plot_spectra(lam_pri, lam_pos, gamma_pri, gamma, cfg, fullfile(figdir, 'Figure_A_spectra'));

%% Reference chain and posterior test set
Vf = reference_chain(cfg, C, f, L, y_obs, sigma_obs, v_map, H_map);
idx_te   = round(linspace(1, size(Vf,2), cfg.n_test));
Theta_te = cfg.mu_pr + L*Vf(:, idx_te);
G_te     = zeros(cfg.m, cfg.n_test);
for k = 1:cfg.n_test
    G_te(:,k) = fem_solve(Theta_te(:,k), cfg, C, f);
end

%% Surrogate error (Figure 1)
tic;
res = surrogate_error_sweep(cfg, C, f, L, Phi, Xi, gamma, Theta_te, G_te, sigma_obs);
fprintf('surrogate sweep: %.1f s\n', toc);
fprintf('RMSE(G_hat)/sigma_obs on posterior test points, median over %d designs\n', cfg.n_rep);
hdr = arrayfun(@(r) sprintf('LIS+GP r=%d', r), cfg.r_list, 'UniformOutput', false);
fprintf('%6s %12s', 'n_tr', 'GP alone');
fprintf(' %12s', hdr{:});
fprintf('\n');
for b = 1:numel(cfg.ntr_list)
    fprintf('%6d %12.2f', cfg.ntr_list(b), median(res.alone(:,b)));
    fprintf(' %12.2f', squeeze(median(res.lis(:,b,:), 1))');
    fprintf('\n');
end
fprintf('%6s %12s', 'floor', '-');
fprintf(' %12.2f', res.floor);
fprintf('\n');
res.r_main = sum(gamma >= cfg.tau_g);
plot_surrogate_error(res, cfg, fullfile(figdir, 'Figure_1_surrogate_error'));
n0 = budget_at(cfg.ntr_list, median(res.alone, 1), 1);
fprintf('n_tr for RMSE = sigma_obs: GP alone %.0f', n0);
for k = 1:numel(cfg.r_list)
    nk = budget_at(cfg.ntr_list, median(res.lis(:,:,k), 1), 1);
    fprintf(', LIS+GP r=%d %.0f (x%.1f)', cfg.r_list(k), nk, n0/nk);
end
fprintf('\n');

%% Posterior mean and covariance (Figure 2)
tic;
post = posterior_sweep(cfg, C, f, L, Phi, Xi, gamma, v_map, H_map, y_obs, sigma_obs, res.r_main);
fprintf('posterior sweep: %.1f s\n', toc);

ref.m = mean(Vf, 2);
ref.S = cov(Vf');
nB = numel(cfg.ntr_list);
E = struct('mean_lisgp', zeros(cfg.n_rep, nB), 'cov_lisgp', zeros(cfg.n_rep, nB), ...
           'mean_alone', zeros(cfg.n_rep, nB), 'cov_alone', zeros(cfg.n_rep, nB));
for rep = 1:cfg.n_rep
    for b = 1:nB
        [E.mean_lisgp(rep,b), E.cov_lisgp(rep,b)] = posterior_error(post.lisgp{rep,b}.m, post.lisgp{rep,b}.S, ref, L);
        [E.mean_alone(rep,b), E.cov_alone(rep,b)] = posterior_error(post.alone{rep,b}.m, post.alone{rep,b}.S, ref, L);
    end
end
[E.mean_lis,  E.cov_lis]  = posterior_error(post.lis.m,  post.lis.S,  ref, L);
[E.mean_full, E.cov_full] = posterior_error(post.full.m, post.full.S, ref, L);

fprintf('posterior error vs reference chain, median over %d designs\n', cfg.n_rep);
fprintf('%6s %12s %12s %12s %12s\n', 'n_tr', 'e_mu alone', 'e_mu LIS+GP', 'e_G alone', 'e_G LIS+GP');
for b = 1:nB
    fprintf('%6d %12.3f %12.3f %12.3f %12.3f\n', cfg.ntr_list(b), median(E.mean_alone(:,b)), ...
            median(E.mean_lisgp(:,b)), median(E.cov_alone(:,b)), median(E.cov_lisgp(:,b)));
end
fprintf('%6s %12.3f %12.3f %12.3f %12.3f\n', 'exact', E.mean_full, E.mean_lis, E.cov_full, E.cov_lis);
plot_posterior_error(E, cfg, post.r, fullfile(figdir, 'Figure_2_posterior_error'));

%% Covariance maps (Figure 3 and Figure 4)
b = find(cfg.ntr_list == cfg.ntr_map, 1);
Gm.pr    = L*L';
Gm.ref   = L*ref.S*L';
Gm.lis   = L*post.lis.S*L';
Gm.lisgp = L*post.lisgp{1,b}.S*L';
Gm.alone = L*post.alone{1,b}.S*L';
plot_covariance_maps(Gm, cfg.ntr_map, fullfile(figdir, 'Figure_3_covariances'), ...
                     fullfile(figdir, 'Figure_4_covariance_errors'));
