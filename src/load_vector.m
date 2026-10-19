% Nodal vector f from the distributed q(z) 
% constant + Gaussian bump
function f = load_vector(cfg)

q  = cfg.q0 + cfg.q_bump * exp(-(cfg.x_mid - cfg.bump_centre).^2 ...
                               / (2 * cfg.bump_width^2));
fe = q * cfg.l / 2;

f = zeros(cfg.n, 1);
f(1) = fe(1);
f(1:cfg.n-1) = f(1:cfg.n-1) + fe(2:cfg.n);
f(2:cfg.n)   = f(2:cfg.n)   + fe(2:cfg.n);
end
