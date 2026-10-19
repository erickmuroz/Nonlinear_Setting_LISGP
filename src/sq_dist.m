function D2 = sq_dist(A, B)
D2 = sum(A.^2, 1)' + sum(B.^2, 1) - 2*(A'*B);
D2 = max(D2, 0);
end
