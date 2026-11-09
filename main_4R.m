%% PROGETTO: Robot manipolatore planare 4R con gravita'
% Script principale: FK, IK, dinamica, controllore a coppia calcolata,
% simulazione con ode45, grafici e animazione.
clear; clc; close all;

%% Parametri del robot
params.L  = [0.5 0.4 0.3 0.2];          % lunghezze link [m]
params.m  = [2.0 1.5 1.0 0.5];          % masse [kg]
params.lc = params.L / 2;               % baricentri a meta' link
params.I  = params.m .* params.L.^2 / 12;   % inerzie baricentriche [kg m^2]
params.g  = 9.81;                       % gravita' lungo -y [m/s^2]
L = params.L;

%% 1) Cinematica diretta
theta_test = deg2rad([30 -40 25 15]);
[T, x, y, phi] = fk_4R(theta_test, L);
disp('T = A1*A2*A3*A4 ='); disp(T);
fprintf('FK: x = %.4f, y = %.4f, phi = %.4f rad\n\n', x, y, phi);
assert(abs(phi - sum(theta_test)) < 1e-12, 'phi deve essere la somma dei giunti');
assert(norm(T(1:3, 1:3)' * T(1:3, 1:3) - eye(3)) < 1e-12, 'La parte rotazionale di T non e'' ortogonale');

figure('Name', 'Cinematica diretta');
plot_robot_4R(gca, theta_test, L);
title('Robot 4R e terne DH (configurazione di test)');

%% 2) Cinematica inversa: FK(IK(p)) = p, gomito basso e alto
figure('Name', 'Cinematica inversa'); hold on;
err_ik = [];
for el = [1 -1]
    opts = struct('elbow', el, 'theta4', theta_test(4));
    th = ik_4R(x, y, phi, L, opts);
    [~, xc, yc, phic] = fk_4R(th, L);
    err = norm([x - xc, y - yc, phi - phic]);
    assert(err < 1e-10, 'FK(IK(p)) diverso da p (gomito %d)', el);
    err_ik(end+1) = err;
    fprintf('IK elbow = %+d: theta = [%s] deg, errore FK(IK(p)) = %.2e\n', ...
        el, num2str(rad2deg(th'), '%8.2f'), err);
    subplot(1, 2, (el == -1) + 1);
    plot_robot_4R(gca, th, L);
    title(sprintf('Gomito %s', ternary(el == 1, 'basso', 'alto')));
end
fprintf('\n');

% Un punto fuori dallo spazio di lavoro deve generare un errore esplicito
try
    ik_4R(sum(L) + 1, 0, 0, L);
    error('Test irraggiungibilita'' fallito: nessun errore sollevato');
catch ME
    assert(strcmp(ME.identifier, 'ik_4R:Unreachable'), 'Errore inatteso: %s', ME.message);
    disp('Test punto irraggiungibile superato.');
end

%% 3) Test della dinamica: conservazione dell'energia (moto libero senza gravita' attiva sui giunti)
[Mt, Ct, Gt] = dynamics_4R(deg2rad([20 -30 40 -10])', [0.3; -0.2; 0.5; 0.1], params);
assert(norm(Mt - Mt') < 1e-12 && all(eig(Mt) > 0), 'M deve essere simmetrica definita positiva');
assert(all(isfinite([Mt(:); Ct(:); Gt(:)])), 'M, C, G contengono NaN o Inf');

x0 = [deg2rad([20 -30 40 -10])'; 0; 0; 0; 0];
free = @(~, s) [s(5:8); ...
    local_acc(s(1:4), s(5:8), zeros(4, 1), params)];
[tt, ss] = ode45(free, [0 3], x0, odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
E = zeros(numel(tt), 1);
for k = 1:numel(tt)
    Mk = dynamics_4R(ss(k, 1:4)', ss(k, 5:8)', params);
    E(k) = 0.5 * ss(k, 5:8) * Mk * ss(k, 5:8)' + potential_energy(ss(k, 1:4)', params);
end
fprintf('Test energia: variazione massima = %.2e J (E0 = %.4f J)\n\n', max(abs(E - E(1))), E(1));
assert(max(abs(E - E(1))) < 1e-6, 'L''energia non si conserva: modello dinamico errato');
figure('Name', 'Test energia');
plot(tt, E - E(1)); grid on; xlabel('t [s]'); ylabel('E - E_0 [J]');
title('Conservazione dell''energia (moto libero)');

%% 4) Traiettoria di riferimento: cerchio nello spazio operativo
Tend = 8;
tg = 0:0.01:Tend;
xc0 = 0.7; yc0 = 0.4; R = 0.15; w = 2*pi/Tend;
phi_d = 0.5;
elbow = -1;                      % gomito alto
opts = struct('elbow', elbow, 'theta4', 0.3);

qref = zeros(4, numel(tg));
for k = 1:numel(tg)
    qref(:, k) = ik_4R(xc0 + R*cos(w*tg(k)), yc0 + R*sin(w*tg(k)), phi_d, L, opts);
end
qref = unwrap(qref, [], 2);

pp   = spline(tg, qref);
ppd  = ppder(pp);
ppdd = ppder(ppd);
qd_f   = @(t) ppval(pp, t);
dqd_f  = @(t) ppval(ppd, t);
ddqd_f = @(t) ppval(ppdd, t);

%% 5) Controllore a coppia calcolata (inverse dynamics + PD)
%   tau = M(q) (ddq_d + Kd de + Kp e) + C(q,dq) dq + G(q)
%   => errore: de'' + Kd de' + Kp e = 0 (asintoticamente stabile)
Kp = diag([100 100 100 100]);
Kd = 2 * sqrt(Kp);               % smorzamento critico
ctrl = @(t, q, dq) control_law(t, q, dq, qd_f, dqd_f, ddqd_f, Kp, Kd, params);

q0  = qd_f(0) + deg2rad([8; -8; 6; -6]);   % errore iniziale
dq0 = dqd_f(0);
sim_f = @(t, s) [s(5:8); local_acc(s(1:4), s(5:8), ctrl(t, s(1:4), s(5:8)), params)];
[t, S] = ode45(sim_f, tg, [q0; dq0], odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
q = S(:, 1:4)';
dq = S(:, 5:8)';

qd = qref;
e = qd - q;
tau = zeros(4, numel(t));
pos = zeros(2, numel(t));
pos_d = zeros(2, numel(t));
for k = 1:numel(t)
    tau(:, k) = ctrl(t(k), q(:, k), dq(:, k));
    [~, pos(1, k), pos(2, k)] = fk_4R(q(:, k), L);
    [~, pos_d(1, k), pos_d(2, k)] = fk_4R(qd(:, k), L);
end
fprintf('Controllore: errore finale ||e|| = %.2e rad, errore max dopo 2 s = %.2e rad\n', ...
    norm(e(:, end)), max(vecnorm(e(:, t > 2))));
fprintf('Picchi delle coppie [N m]: %s\n', num2str(max(abs(tau), [], 2)', '%8.2f'));
assert(max(vecnorm(e(:, t > 2))) < 1e-2, 'Il controllore non converge sulla traiettoria');

% Riepilogo da copiare nel README
fprintf('\n===== RIEPILOGO VALIDAZIONE =====\n');
fprintf('Versione MATLAB: %s\n', version);
fprintf('Errore FK(IK(p)) gomito basso / alto: %.3e / %.3e\n', err_ik(1), err_ik(2));
fprintf('Energia, variazione massima in moto libero: %.3e J\n', max(abs(E - E(1))));
fprintf('Errore di giunto iniziale ||e||: %.3e rad\n', norm(e(:, 1)));
fprintf('Errore di giunto massimo dopo 2 s: %.3e rad\n', max(vecnorm(e(:, t > 2))));
fprintf('Errore di giunto finale ||e||: %.3e rad\n', norm(e(:, end)));
fprintf('Errore finale end-effector: %.3e m\n', norm(pos(:, end) - pos_d(:, end)));
fprintf('Picchi delle coppie [N m]: %s\n', num2str(max(abs(tau), [], 2)', '%8.2f'));
fprintf('=================================\n\n');

%% 6) Grafici
lbl = {'1', '2', '3', '4'};
figure('Name', 'Risultati controllo', 'Position', [100 100 1100 700]);
subplot(2, 2, 1);
plot(t, q', 'LineWidth', 1.2); hold on; plot(t, qd', '--k');
grid on; xlabel('t [s]'); ylabel('q [rad]'); title('q(t) (tratteggio = riferimento)');
legend(strcat('q_', lbl), 'Location', 'best');

subplot(2, 2, 2);
plot(t, e', 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('e [rad]'); title('Errore e(t) = q_d - q');
legend(strcat('e_', lbl), 'Location', 'best');

subplot(2, 2, 3);
plot(t, tau', 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('\tau [N m]'); title('Coppie \tau(t)');
legend(strcat('\tau_', lbl), 'Location', 'best');

subplot(2, 2, 4);
plot(pos_d(1, :), pos_d(2, :), 'k--', 'LineWidth', 1.5); hold on;
plot(pos(1, :), pos(2, :), 'm', 'LineWidth', 1.2);
axis equal; grid on; xlabel('x [m]'); ylabel('y [m]');
title('Traiettoria end-effector'); legend('desiderata', 'effettiva');


%% 7) Animazione
fig = figure('Name', 'Animazione', 'Position', [150 150 700 650]);
ax = axes(fig);
reach = sum(L) * 1.05;
step = 3;
for k = 1:step:numel(t)
    cla(ax);
    plot(ax, pos_d(1, :), pos_d(2, :), 'k--', 'LineWidth', 1);     % target
    hold(ax, 'on');
    plot(ax, pos_d(1, k), pos_d(2, k), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    plot(ax, pos(1, 1:k), pos(2, 1:k), 'm', 'LineWidth', 1.2);    % scia
    plot_robot_4R(ax, q(:, k), L);
    axis(ax, [-0.3 reach -0.3 reach]);
    title(ax, sprintf('t = %.2f s', t(k)));
    drawnow;
end

%if ~exist('results', 'dir'), mkdir('results'); end

%% 8) Salvataggio delle figure in results/
%nomi = {'Cinematica diretta', 'Cinematica inversa', 'Test energia', ...
%        'Risultati controllo', 'Animazione'};
%file = {'cinematica_diretta', 'cinematica_inversa', 'test_energia', ...
%        'risultati_controllo', 'animazione_fine'};
%for i = 1:numel(nomi)
%    f = findobj('Type', 'figure', 'Name', nomi{i});
%    exportgraphics(f(1), fullfile('results', [file{i} '.png']), 'Resolution', 150);
%end
%fprintf('Figure salvate nella cartella results/\n');

%% Funzioni locali
function tau = control_law(t, q, dq, qd_f, dqd_f, ddqd_f, Kp, Kd, params)
e  = qd_f(t) - q;
de = dqd_f(t) - dq;
[M, C, G] = dynamics_4R(q, dq, params);
tau = M * (ddqd_f(t) + Kd * de + Kp * e) + C * dq + G;
end

function ddq = local_acc(q, dq, tau, params)
[M, C, G] = dynamics_4R(q, dq, params);
ddq = M \ (tau - C * dq - G);
end

function V = potential_energy(q, params)
ph = cumsum(q(:));
V = 0; yj = 0;
for i = 1:4
    yj_c = yj + params.lc(i) * sin(ph(i));
    V = V + params.m(i) * params.g * yj_c;
    yj = yj + params.L(i) * sin(ph(i));
end
end

function ppd = ppder(pp)
[b, c, l, k, d] = unmkpp(pp);
ppd = mkpp(b, repmat(k-1:-1:1, d*l, 1) .* c(:, 1:k-1), d);
end

function r = ternary(c, a, b)
if c, r = a; else, r = b; end
end
