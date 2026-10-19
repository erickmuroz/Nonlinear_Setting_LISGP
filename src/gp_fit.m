
function gp = gp_fit(Xtr, Ytr, jitter)
nt = size(Xtr, 2);
d  = size(Ytr, 1);

D2  = sq_dist(Xtr, Xtr);
off = D2(~eye(nt));
ell0 = sqrt(median(off(off > 0)) / 2);

best = Inf;
ell  = ell0;
for e = ell0 * logspace(-1.2, 1.5, 30)
    K = exp(-D2/(2*e^2)) + jitter*eye(nt);
    [R, p] = chol(K);
    if p > 0, continue; end
    logdet = 2*sum(log(diag(R)));
    nlml = 0;
    for k = 1:d
        yc = Ytr(k,:)' - mean(Ytr(k,:));
        a  = R \ (R' \ yc);
        s2 = max((yc'*a)/nt, realmin);
        nlml = nlml + 0.5*nt*log(s2) + 0.5*logdet;
    end
    if nlml < best
        best = nlml;
        ell  = e;
    end
end

Ktr   = exp(-D2/(2*ell^2));
alpha = zeros(nt, d);
sf2   = zeros(d, 1);
m0    = zeros(d, 1);
for k = 1:d
    yk     = Ytr(k,:)';
    m0(k)  = mean(yk);
    sf2(k) = max(var(yk), eps);
    R = chol(sf2(k)*Ktr + jitter*sf2(k)*eye(nt), 'lower');
    alpha(:,k) = R' \ (R \ (yk - m0(k)));
end
gp = struct('Xtr', Xtr, 'ell', ell, 'sf2', sf2, 'm0', m0, 'alpha', alpha);
end
