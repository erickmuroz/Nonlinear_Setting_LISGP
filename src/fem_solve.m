% Solve of full system
function [y, u] = fem_solve(theta, cfg, C, f)

K = assemble_stiffness(theta, cfg);
u = K \ f;
y = C * u;
end
