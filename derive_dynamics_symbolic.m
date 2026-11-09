%% DERIVE_DYNAMICS_SYMBOLIC  Derivazione simbolica Eulero-Lagrange del 4R planare
% Richiede Symbolic Math Toolbox. Ricava M(q), C(q,dq), G(q) e li confronta
% numericamente con dynamics_4R.
clear; clc;

n = 4;
q  = sym('q',  [n 1], 'real');
dq = sym('dq', [n 1], 'real');
L  = sym('L',  [n 1], 'positive');
lc = sym('lc', [n 1], 'positive');
m  = sym('m',  [n 1], 'positive');
I  = sym('I',  [n 1], 'positive');
syms g positive

%% Posizione dei baricentri
ph = cumsum(q);
pc = sym(zeros(2, n));
pj = [0; 0];                       % origine del giunto i
for i = 1:n
    u = [cos(ph(i)); sin(ph(i))];
    pc(:, i) = pj + lc(i) * u;
    pj = pj + L(i) * u;
end

%% Energia cinetica T = 1/2 dq' M dq  ->  M(q)
M = sym(zeros(n));
for i = 1:n
    Jv = jacobian(pc(:, i), q);
    Jw = [ones(1, i), zeros(1, n - i)];
    M = M + m(i) * (Jv.' * Jv) + I(i) * (Jw.' * Jw);
end
M = simplify(M);

%% Energia potenziale V = sum m_i g h_i  ->  G(q)
V = sum(m .* g .* pc(2, :).');
G = simplify(gradient(V, q));

%% Coriolis / centrifughi (simboli di Christoffel)
C = sym(zeros(n));
for k = 1:n
    for j = 1:n
        for r = 1:n
            C(k, j) = C(k, j) + 1/2 * (diff(M(k, j), q(r)) + diff(M(k, r), q(j)) ...
                      - diff(M(j, r), q(k))) * dq(r);
        end
    end
end
C = simplify(C);

disp('M(q) ='); disp(M);
disp('G(q) ='); disp(G);
disp('C(q,dq) ='); disp(C);

%% Verifica: Mdot - 2C antisimmetrica
Mdot = sym(zeros(n));
for k = 1:n
    Mdot = Mdot + diff(M, q(k)) * dq(k);
end
N = Mdot - 2*C;
rnd = @() rand(n, 1);
vars = [q; dq; L; lc; m; I; g];
vals = [rnd(); rnd(); 0.5*ones(n,1)+rnd()*0.1; 0.2*ones(n,1); 1+rnd(); 0.05*ones(n,1); 9.81];
Nn = double(subs(N, vars, vals));
fprintf('||N + N''|| (antisimmetria) = %.2e\n', norm(Nn + Nn'));

%% Confronto numerico con dynamics_4R
params.L  = [0.5 0.4 0.3 0.2];
params.m  = [2 1.5 1 0.5];
params.lc = params.L / 2;
params.I  = params.m .* params.L.^2 / 12;
params.g  = 9.81;
qn  = [0.3; -0.7; 0.5; 0.9];
dqn = [0.4; 0.2; -0.6; 0.8];

vars = [q; dq; L; lc; m; I; g];
vals = [qn; dqn; params.L(:); params.lc(:); params.m(:); params.I(:); params.g];
Ms = double(subs(M, vars, vals));
Cs = double(subs(C, vars, vals));
Gs = double(subs(G, vars, vals));
[Mn, Cn, Gn] = dynamics_4R(qn, dqn, params);

fprintf('Errore M: %.2e\n', norm(Ms - Mn));
fprintf('Errore C: %.2e\n', norm(Cs - Cn));
fprintf('Errore G: %.2e\n', norm(Gs - Gn));

save('dynamics_symbolic.mat', 'M', 'C', 'G');
