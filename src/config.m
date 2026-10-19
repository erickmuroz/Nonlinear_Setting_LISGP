% All parameters (mesh, load, prior, sensors, seeds, budgets)
function cfg = config()

cfg.L = 2.0;
cfg.n = 50;
cfg.l = cfg.L / cfg.n;
cfg.x_mid = (cfg.l/2 : cfg.l : cfg.L - cfg.l/2)';

cfg.q0          = 4e6;
cfg.q_bump      = 2e6;
cfg.bump_width  = cfg.L / 8;
cfg.bump_centre = cfg.L / 2;

cfg.m           = 6;
cfg.sensor_seed = 41;
cfg.SNR         = 50;

cfg.EA_ref      = 1e7;
cfg.mu_pr       = log(cfg.EA_ref);
cfg.sigma_theta = 0.60;
cfg.ell_pr      = 0.50;

cfg.dip_amplitude = 0.50;
cfg.dip_width     = 0.10;
cfg.dip_centre    = cfg.L / 2;

cfg.data_seed = 2718;
cfg.map_tol   = 1e-10;
cfg.map_maxit = 50;

cfg.tau_loc     = 0.1;
cfg.tau_g       = 0.1;
cfg.N_prior     = 2000;
cfg.sample_seed = 7;

cfg.lis_iters = 500;
cfg.lis_L     = 200;
cfg.lis_seed  = 555;

cfg.mcmc_steps_full = 300000;
cfg.mcmc_pilot      = 2000;
cfg.mcmc_burnin     = 0.20;
cfg.mcmc_acc_target = 0.30;
cfg.mcmc_seed       = 4242;

cfg.n_test    = 300;
cfg.ntr_list  = [20 40 80 160 320 640 1280 2560];
cfg.r_list    = [4 8 10];
cfg.n_rep     = 3;
cfg.design_seed = 1000;
cfg.gp_jitter = 1e-8;
cfg.ref_thin  = 10;
cfg.ref_cache = 'reference_chain.mat';
cfg.chain_steps = 200000;
cfg.ntr_map   = 160;
end
