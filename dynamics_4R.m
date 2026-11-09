function [M, C, G] = dynamics_4R(theta, dtheta, params)
% Termini dinamici del 4R planare: M(q) ddq + C(q,dq) dq + G(q) = tau
% M e G si ottengono dai Jacobiani dei baricentri, C dai simboli di Christoffel.
% Il risultato e' confrontato con derive_dynamics_symbolic.m

n = 4;
theta = theta(:);
dtheta = dtheta(:);
L = params.L(:);
m = params.m(:);
lc = params.lc(:);
I = params.I(:);
g = params.g;

ph = cumsum(theta);
c = cos(ph);
s = sin(ph);

M = zeros(n);
G = zeros(n, 1);
dM = zeros(n, n, n);   % dM(:,:,k) = dM/dq_k

for i = 1:n
    a = [L(1:i-1); lc(i)];          % distanze lungo i link fino al baricentro i
    Jv = zeros(2, n);
    dJ = zeros(2, n, n);            % dJ(:,j,k) = d(Jv(:,j))/dq_k
    for j = 1:i
        for l = j:i
            Jv(:, j) = Jv(:, j) + a(l) * [-s(l); c(l)];
        end
        for k = 1:n
            for l = max(j, k):i
                dJ(:, j, k) = dJ(:, j, k) - a(l) * [c(l); s(l)];
            end
        end
    end
    Jw = [ones(1, i), zeros(1, n - i)];

    M = M + m(i) * (Jv' * Jv) + I(i) * (Jw' * Jw);
    G = G + m(i) * g * Jv(2, :)';
    for k = 1:n
        dJk = dJ(:, :, k);
        dM(:, :, k) = dM(:, :, k) + m(i) * (dJk' * Jv + Jv' * dJk);
    end
end

C = zeros(n);
for k = 1:n
    for j = 1:n
        for r = 1:n
            C(k, j) = C(k, j) + 0.5 * (dM(k, j, r) + dM(k, r, j) - dM(j, r, k)) * dtheta(r);
        end
    end
end
end
