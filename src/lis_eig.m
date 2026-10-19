function [Psi, gamma] = lis_eig(S)
[Psi, D] = eig((S + S') / 2);
[gamma, ord] = sort(real(diag(D)), 'descend');
Psi = Psi(:, ord);
end
