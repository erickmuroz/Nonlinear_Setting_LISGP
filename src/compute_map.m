% compute and store for later 
function [theta_map, v_map, hist] = compute_map(cfg, C, f, L, y_obs, sigma_obs)

v = zeros(cfg.n, 1);
[F, Jw, res] = map_objective(v, cfg, C, f, L, y_obs, sigma_obs);

hist.F        = F;
hist.step     = 0;
hist.n_backtrack = 0;

for it = 1:cfg.map_maxit
    g  = Jw'*res + v;
    Hn = Jw'*Jw + eye(cfg.n);
    dv = -(Hn \ g);

    alpha = 1;
    nbt   = 0;
    while true
        v_try = v + alpha*dv;
        F_try = map_objective(v_try, cfg, C, f, L, y_obs, sigma_obs);
        if F_try < F || nbt >= 30
            break;
        end
        alpha = alpha/2;
        nbt   = nbt + 1;
    end

    F_old = F;
    v = v + alpha*dv;
    [F, Jw, res] = map_objective(v, cfg, C, f, L, y_obs, sigma_obs);

    hist.F(end+1)    = F;
    hist.step(end+1) = norm(alpha*dv);
    hist.n_backtrack(end+1) = nbt;

    if abs(F_old - F) <= cfg.map_tol * max(abs(F_old), 1)
        break;
    end
end

hist.iters = it;
v_map      = v;
theta_map  = cfg.mu_pr + L*v_map;
end
