function [X, acc, info] = mcmc_rw(logpi, x0, R_prop, n_steps, opts)

d = numel(x0);

rng(opts.seed, 'twister');

[Xp, acc_pilot] = rw_core(logpi, x0(:), R_prop, opts.s0, opts.n_pilot);

ratio = acc_pilot / opts.acc_target;
s     = opts.s0 * min(2, max(0.5, sqrt(ratio)));
x_start = Xp(:, end);

[X, acc] = rw_core(logpi, x_start, R_prop, s, n_steps);

info = struct('s0', opts.s0, 's', s, 'acc_pilot', acc_pilot, ...
              'n_pilot', opts.n_pilot, 'd', d);
end

function [X, acc] = rw_core(logpi, x, R_prop, s, n_steps)

d  = numel(x);
X  = zeros(d, n_steps);
lp = logpi(x);
n_acc = 0;

for k = 1:n_steps
    xp  = x + s * (R_prop' * randn(d,1));
    lpp = logpi(xp);
    if log(rand) < lpp - lp
        x  = xp;
        lp = lpp;
        n_acc = n_acc + 1;
    end
    X(:,k) = x;
end

acc = n_acc / n_steps;
end
