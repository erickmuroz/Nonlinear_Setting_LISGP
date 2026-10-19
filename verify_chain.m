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

%% 1. Cost of one forward solve
th = cfg.mu_pr + L*randn(cfg.n, 1);
tic;
for k = 1:1e5
    fem_solve(th, cfg, C, f);
end
t_solve = toc / 1e5;
fprintf('1. one forward solve: %.2f us, expected chain time for %d steps: %.1f s\n', ...
        1e6*t_solve, cfg.mcmc_steps_full, t_solve*cfg.mcmc_steps_full);

%% 2. Two independent chains: one from the MAP, one from a prior draw
logpi = @(v) -0.5*sum(((fem_solve(cfg.mu_pr + L*v, cfg, C, f) - y_obs)/sigma_obs).^2) - 0.5*(v'*v);
R_full = chol(inv((H_map + H_map')/2));
opts = struct('s0', 2.38/sqrt(cfg.n), 'n_pilot', cfg.mcmc_pilot, ...
              'acc_target', cfg.mcmc_acc_target, 'seed', cfg.mcmc_seed);
burn = round(cfg.mcmc_burnin * cfg.mcmc_steps_full);

tic;
[VA, accA, infoA] = mcmc_rw(logpi, v_map, R_full, cfg.mcmc_steps_full, opts);
tA = toc;
opts.seed = cfg.mcmc_seed + 1;
rng(123, 'twister');
v0 = randn(cfg.n, 1);
tic;
[VB, accB, infoB] = mcmc_rw(logpi, v0, R_full, cfg.mcmc_steps_full, opts);
tB = toc;
fprintf('2. chain A (from MAP): %.1f s, step %.3f, acceptance %.3f\n', tA, infoA.s, accA);
fprintf('   chain B (from prior draw): %.1f s, step %.3f, acceptance %.3f\n', tB, infoB.s, accB);

fprintf('   distinct consecutive states: A %.3f, B %.3f (should equal acceptance)\n', ...
        mean(any(diff(VA,1,2) ~= 0, 1)), mean(any(diff(VB,1,2) ~= 0, 1)));

VA = VA(:, burn+1:end);
VB = VB(:, burn+1:end);

%% 3. Convergence on the 8 most data-informed whitened directions
[Q, D] = eig(Jw_map'*Jw_map);
[~, ord] = sort(real(diag(D)), 'descend');
Q8 = Q(:, ord(1:8));
ZA = (Q8' * VA)';
ZB = (Q8' * VB)';

N  = size(ZA, 1);
h  = floor(N/2);
fprintf('\n3. direction   R-hat    ESS(A)   ESS(B)   mean A    mean B    std A     std B\n');
for j = 1:8
    S = [ZA(1:h,j), ZA(h+1:2*h,j), ZB(1:h,j), ZB(h+1:2*h,j)];
    fprintf('   %5d     %7.4f  %7.0f  %7.0f  %8.4f  %8.4f  %8.4f  %8.4f\n', j, rhat(S), ...
            ess(ZA(:,j)), ess(ZB(:,j)), mean(ZA(:,j)), mean(ZB(:,j)), std(ZA(:,j)), std(ZB(:,j)));
end

%% 4. Data misfit along the chains
res2 = @(V) mean(arrayfun(@(k) sum(((fem_solve(cfg.mu_pr + L*V(:,k), cfg, C, f) - y_obs)/sigma_obs).^2), ...
                          round(linspace(1, size(V,2), 2000))));
fprintf('\n4. mean ||(G(theta) - y_obs)/sigma||^2: A %.3f, B %.3f, at MAP %.3f (m = %d)\n', ...
        res2(VA), res2(VB), sum(((fem_solve(cfg.mu_pr + L*v_map, cfg, C, f) - y_obs)/sigma_obs).^2), cfg.m);

%% 5. Posterior mean and covariance of the two chains in theta
muA = cfg.mu_pr + L*mean(VA, 2);
muB = cfg.mu_pr + L*mean(VB, 2);
fprintf('5. relative difference of the posterior means (A vs B): %.3e\n', norm(muA - muB)/norm(muA - cfg.mu_pr));
CA = L*cov(VA')*L';
CB = L*cov(VB')*L';
fprintf('   relative difference of the posterior covariances (A vs B): %.3e\n', norm(CA - CB, 'fro')/norm(CA, 'fro'));

function r = rhat(S)
[n, m] = size(S);
W = mean(var(S, 0, 1));
B = n * var(mean(S, 1));
r = sqrt(((n-1)/n * W + B/n) / W);
end

function e = ess(x)
x = x(:) - mean(x);
N = numel(x);
F = fft(x, 2^nextpow2(2*N));
ac = real(ifft(F .* conj(F)));
ac = ac(1:N) / ac(1);
tau = 1;
for k = 1:2:N-2
    p = ac(k+1) + ac(k+2);
    if p < 0, break; end
    tau = tau + 2*p;
end
e = N / tau;
end
