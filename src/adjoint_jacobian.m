
function J = adjoint_jacobian(theta, cfg, C, f)

n = cfg.n;
K = assemble_stiffness(theta, cfg);
u = K \ f;
A = K \ C';
ke = exp(theta(:)) / cfg.l;

uext = [0; u];
Aext = [zeros(1, cfg.m); A];

delta = uext(1:n) - uext(2:n+1);
Dmat  = Aext(1:n, :) - Aext(2:n+1, :);

J = -(Dmat') .* ((ke .* delta)');
end
