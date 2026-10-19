% assembly of the K matrix 
function K = assemble_stiffness(theta, cfg)

n  = cfg.n;
ke = exp(theta(:)) / cfg.l;

main_diag = [ke(1:n-1) + ke(2:n); ke(n)];
off_diag  = -ke(2:n);

K = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);
end
