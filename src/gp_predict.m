function Y = gp_predict(gp, Xte)
K = exp(-sq_dist(Xte, gp.Xtr) / (2*gp.ell^2));
Y = gp.m0 + gp.sf2 .* (K * gp.alpha)';
end
