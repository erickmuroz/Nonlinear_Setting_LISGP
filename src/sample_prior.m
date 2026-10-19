% how prior draws are built theta = mu_pr + L z
function Theta = sample_prior(N, cfg, L, seed)
rng(seed, 'twister');
Theta = cfg.mu_pr + L * randn(cfg.n, N);
end
