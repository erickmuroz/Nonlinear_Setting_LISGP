function d = weighted_subspace_distance(Pa, ga, Pb, gb, tau)
ka = max(1, sum(ga >= tau));
kb = max(1, sum(gb >= tau));
A  = Pa(:,1:ka) .* sqrt(ga(1:ka)' / sum(ga(1:ka)));
B  = Pb(:,1:kb) .* sqrt(gb(1:kb)' / sum(gb(1:kb)));
d  = sqrt(max(0, 1 - norm(A'*B, 'fro')^2 / (norm(A'*A, 'fro') * norm(B'*B, 'fro'))));
end
