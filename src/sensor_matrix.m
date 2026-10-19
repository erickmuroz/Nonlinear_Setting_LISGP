% Sensor operator m x n , m =6 chosen at random (seed 41)
function [C, idx] = sensor_matrix(cfg)

rng(cfg.sensor_seed, 'twister');
perm = randperm(cfg.n);
idx  = sort(perm(1:cfg.m))';

C = zeros(cfg.m, cfg.n);
C(sub2ind(size(C), (1:cfg.m)', idx)) = 1;
end
