function lambda = local_spectrum(theta, cfg, C, f, L, sigma_obs)
J = adjoint_jacobian(theta, cfg, C, f);
lambda = svd((J * L) / sigma_obs).^2;
end
