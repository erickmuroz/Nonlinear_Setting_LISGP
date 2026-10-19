% Creates the measured data
% used for MAP, the reference chain, the LIS, and surrogate errors 
% fixed seed: conclusions are for a single noise realization
function y_obs = generate_synthetic_data(cfg, C, f, theta_true, sigma_obs)

rng(cfg.data_seed, 'twister');
y_clean = fem_solve(theta_true, cfg, C, f);
y_obs   = y_clean + sigma_obs * randn(cfg.m, 1);
end
