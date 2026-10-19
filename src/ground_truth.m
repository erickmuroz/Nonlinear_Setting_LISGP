% Ground truth - stiffness dip at midspan + noise level sigma_obs
% what we would like to recover, not exactly reproducible 
% theta_true is the log-stiffnes field we pretend generated the data
% see MAP || theta_map - theta_true || = 0.87

% SNR =  50 : Is saying noise is 2% of the maximum signal, following Cui et
% al 4.1. It survives any change of load, stiffness or units 

% it is computes once from the noise-free data G(theta_true) and then
% treated as a constant
function [theta_true, sigma_obs] = ground_truth(cfg, C, f)

theta_true = cfg.mu_pr - cfg.dip_amplitude * ...
    exp(-(cfg.x_mid - cfg.dip_centre).^2 / (2 * cfg.dip_width^2));

G_true    = fem_solve(theta_true, cfg, C, f);
sigma_obs = max(abs(G_true)) / cfg.SNR;
end
