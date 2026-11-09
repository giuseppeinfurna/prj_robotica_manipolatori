%% ========================================================================
% SIMULAZIONE DINAMICA E CONTROLLO DI UN ROBOT PLANARE 4R
% ========================================================================
% Questo script simula il comportamento dinamico di un manipolatore planare
% a 4 giunti rotazionali (4R) con gravita' (g lungo -y), sottoposto a:
%   - controllo nello spazio dei giunti (JOINT)
%   - oppure controllo nello spazio operativo (TASK) su (x, y, phi)
%     dell'end-effector. Il robot e' ridondante (4 giunti, 3 task):
%     il grado di liberta' in eccesso e' smorzato nel null-space.
%
% UTILIZZO:
%   - Impostare la variabile CONTROL_MODE a:
%       * "JOINT" per controllo nello spazio dei giunti
%       * "TASK"  per controllo nello spazio operativo
%   - Eseguire lo script in MATLAB. Verranno generati:
%       * Verifiche numeriche (FK/IK, Jacobiano, dinamica)
%       * Grafico della traiettoria
%       * Grafico dell'errore e delle coppie
%       * Animazione del robot 4R
%
% NOTE:
%   - Il modello e' planare, con tutti i giunti rotazionali nel piano (z uscente).
%   - La cinematica e' implementata tramite matrici DH standard.
%   - La dinamica e' M(q) ddq + C(q,dq) + G(q) = tau, con C vettore dei
%     termini centrifughi/Coriolis (simboli di Christoffel).
%   - Controllo a coppia calcolata (inverse dynamics + PD).
%   - L'integrazione della dinamica e' effettuata con Eulero esplicito.
% ========================================================================

clear; close all; clc;

%% =============================================================
% PARAMETRI ROBOT 4R
%% =============================================================
L = [0.5; 0.4; 0.3; 0.2];       % lunghezze dei link [m]
m = [2.0; 1.5; 1.0; 0.5];       % masse [kg]
l = L/2;                        % baricentri a meta' link [m]
I = (1/12).*m.*L.^2;            % inerzie baricentriche [kg*m^2]
g = 9.81;                       % gravita' [m/s^2]

%% =============================================================
% MODALITA' DI CONTROLLO
%% =============================================================
% Possibili:
%   CONTROL_MODE = "TASK"
%   CONTROL_MODE = "JOINT"
CONTROL_MODE = "TASK";
disp("Modalità di controllo attiva: " + CONTROL_MODE);

% Gain joint-space
Kp_j = 100;
Kd_j = 2*sqrt(Kp_j);            % smorzamento critico

% Gain task-space (x, y, phi)
Kp_t = 100;
Kd_t = 2*sqrt(Kp_t);
Kn   = 5;                       % smorzamento nel null-space
lambda = 1e-3;                  % damping della pseudoinversa

%% =============================================================
% VERIFICHE NUMERICHE
%% =============================================================
theta_test = deg2rad([30; -40; 25; 15]);
dtheta_test = [0.4; 0.2; -0.6; 0.8];

% FK e IK: FK(IK(p)) = p, gomito basso e alto
[~, x_t, y_t, phi_t] = fk_4R_DH(theta_test, L);
for elbow = [1 -1]
    th = ik_4R_planar(x_t, y_t, phi_t, L, elbow, theta_test(4));
    [~, xc, yc, phic] = fk_4R_DH(th, L);
    err_ik = norm([x_t - xc; y_t - yc; phi_t - phic]);
    assert(err_ik < 1e-10, "FK(IK(p)) diverso da p");
    fprintf("IK gomito %+d: errore FK(IK(p)) = %.2e\n", elbow, err_ik);
end
disp("Test cinematica superato.")

% Jacobiano confrontato con differenze finite della FK
J = jacobian_4R(theta_test, L);
J_num = zeros(3,4);
h = 1e-6;
for j = 1:4
    dth = zeros(4,1); dth(j) = h;
    [~, xp, yp, php] = fk_4R_DH(theta_test + dth, L);
    [~, xm, ym, phm] = fk_4R_DH(theta_test - dth, L);
    J_num(:,j) = [xp-xm; yp-ym; php-phm]/(2*h);
end
assert(norm(J - J_num) < 1e-6, "Jacobiano non coerente con la FK");
disp("Test Jacobiano superato.")

% Dimensioni, simmetria e definita positivita' di M
M = M_4R(theta_test, L, m, l, I);
C = C_4R(theta_test, dtheta_test, L, m, l, I);
G = g_4R(theta_test, L, m, l, g);
assert(isequal(size(M), [4, 4]), "M deve essere una matrice 4x4");
assert(isequal(size(C), [4, 1]), "C deve essere un vettore 4x1");
assert(isequal(size(G), [4, 1]), "G deve essere un vettore 4x1");
assert(all(isfinite([M(:); C; G])), "M, C o G contengono NaN o Inf");
assert(norm(M - M') < 1e-12, "La matrice M non è simmetrica");
assert(all(eig((M+M')/2) > 0), "M non è definita positiva");
disp("Test M, C, G superato.")

% Conservazione dell'energia in moto libero (integrazione con ode45)
free = @(~, s) [s(5:8); ...
    M_4R(s(1:4), L, m, l, I) \ (-C_4R(s(1:4), s(5:8), L, m, l, I) - g_4R(s(1:4), L, m, l, g))];
s0 = [deg2rad([20; -30; 40; -10]); 0; 0; 0; 0];
[tt, ss] = ode45(free, [0 3], s0, odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
E = zeros(numel(tt), 1);
for k = 1:numel(tt)
    E(k) = 0.5 * ss(k,5:8) * M_4R(ss(k,1:4)', L, m, l, I) * ss(k,5:8)' ...
         + energia_potenziale(ss(k,1:4)', L, m, l, g);
end
dE = max(abs(E - E(1)));
fprintf("Variazione massima dell'energia in moto libero: %.2e J\n", dE);
assert(dE < 1e-6, "L'energia non si conserva: modello dinamico errato");
disp("Test energia superato.")

%% =============================================================
% TRAIETTORIA DESIDERATA (operativa)
%% =============================================================
dt   = 1e-3;
Tend = 8;
t    = 0:dt:Tend;
N    = numel(t);

xc0 = 0.7; yc0 = 0.4; R = 0.15; w = 2*pi/Tend;
phi_des = 0.5;

xd   = [xc0 + R*cos(w*t);  yc0 + R*sin(w*t);  phi_des*ones(1,N)];
dxd  = [-R*w*sin(w*t);      R*w*cos(w*t);      zeros(1,N)];
ddxd = [-R*w^2*cos(w*t);   -R*w^2*sin(w*t);    zeros(1,N)];

%% =============================================================
% TRAIETTORIA NELLO SPAZIO DEI GIUNTI (solo se in modalità JOINT)
%% =============================================================
elbow = -1;                     % gomito alto
theta4_fix = 0.3;               % parametro libero della ridondanza

qd = zeros(4,N);
for k = 1:N
    qd(:,k) = ik_4R_planar(xd(1,k), xd(2,k), xd(3,k), L, elbow, theta4_fix);
end
if any(isnan(qd(:)))
    error("Traiettoria fuori dallo spazio di lavoro raggiungibile.");
end
qd   = unwrap(qd, [], 2);
dqd  = gradient(qd, dt);
ddqd = gradient(dqd, dt);

%% =============================================================
% STATI INIZIALI
%% =============================================================
q  = zeros(4,N);
dq = zeros(4,N);
q(:,1)  = qd(:,1) + deg2rad([8; -8; 6; -6]);    % errore iniziale
dq(:,1) = dqd(:,1);

tau     = zeros(4,N);
pos     = zeros(2,N);
err_j   = zeros(4,N);
err_t   = zeros(3,N);

%% =============================================================
% SIMULAZIONE DINAMICA
%% =============================================================
for k = 1:N
    th  = q(:,k);
    dth = dq(:,k);

    M = M_4R(th, L, m, l, I);
    C = C_4R(th, dth, L, m, l, I);
    G = g_4R(th, L, m, l, g);

    [~, x, y, phi] = fk_4R_DH(th, L);
    pos(:,k)   = [x; y];
    err_j(:,k) = qd(:,k) - th;
    err_t(:,k) = xd(:,k) - [x; y; phi];
    err_t(3,k) = atan2(sin(err_t(3,k)), cos(err_t(3,k)));

    % =============================================================
    % CONTROLLO NELLO SPAZIO DEI GIUNTI
    % =============================================================
    if CONTROL_MODE == "JOINT"
        a = ddqd(:,k) + Kd_j*(dqd(:,k) - dth) + Kp_j*err_j(:,k);

    % =============================================================
    % CONTROLLO NELLO SPAZIO OPERATIVO
    % =============================================================
    else
        J  = jacobian_4R(th, L);
        dJ = (jacobian_4R(th + 1e-6*dth, L) - J) / 1e-6;
        de = dxd(:,k) - J*dth;

        Jp = J' / (J*J' + lambda^2*eye(3));
        a  = Jp*(ddxd(:,k) + Kd_t*de + Kp_t*err_t(:,k) - dJ*dth) ...
           - (eye(4) - Jp*J) * (Kn*dth);
    end

    tau(:,k) = M*a + C + G;

    %% ---- DINAMICA ----
    if k < N
        ddth = M \ (tau(:,k) - C - G);
        dq(:,k+1) = dth + ddth*dt;
        q(:,k+1)  = th  + dth*dt;
    end
end

fprintf("Errore finale: ||e_giunti|| = %.2e rad, ||e_task|| = %.2e\n", ...
    norm(err_j(:,end)), norm(err_t(:,end)));

%% =============================================================
% PLOT TRAIETTORIA
%% =============================================================
figure('Name', 'Traiettoria end-effector');
plot(xd(1,:), xd(2,:), 'k--', 'LineWidth', 1.5); hold on;
plot(pos(1,:), pos(2,:), 'm', 'LineWidth', 1.2);
axis equal; grid on;
xlabel('x [m]'); ylabel('y [m]');
title("Traiettoria end-effector (" + CONTROL_MODE + ")");
legend('desiderata', 'effettiva');

%% =============================================================
% PLOT ERRORE
%% =============================================================
figure('Name', 'Errore e coppie');
subplot(2,1,1);
if CONTROL_MODE == "JOINT"
    plot(t, err_j', 'LineWidth', 1.2);
    ylabel('e = q_d - q [rad]');
    legend('e_1', 'e_2', 'e_3', 'e_4');
else
    plot(t, err_t', 'LineWidth', 1.2);
    ylabel('e = x_d - x');
    legend('e_x [m]', 'e_y [m]', 'e_\phi [rad]');
end
grid on; xlabel('t [s]'); title('Errore');

subplot(2,1,2);
plot(t, tau', 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('\tau [N m]'); title('Coppie');
legend('\tau_1', '\tau_2', '\tau_3', '\tau_4');

%% =============================================================
% ANIMAZIONE
%% =============================================================
figure('Name', 'Animazione', 'Position', [150 150 700 650]);
reach = sum(L) * 1.05;
for k = 1:40:N
    cla; hold on; grid on; axis equal;
    axis([-0.3 reach -0.3 reach]);
    plot(xd(1,:), xd(2,:), 'k--');
    plot(xd(1,k), xd(2,k), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    plot(pos(1,1:k), pos(2,1:k), 'm', 'LineWidth', 1.2);
    plot_4R(q(:,k), L);
    title(sprintf('t = %.2f s', t(k)));
    drawnow;
end

%% =========================================================================
% FUNZIONI DI DINAMICA e CINEMATICA
% =========================================================================

%% Jacobiani dei baricentri e loro derivate parziali
function [Jv, dJv] = jac_com_4R(theta, L, l)
    % Jv(:,j,i)    = d p_ci / d q_j
    % dJv(:,j,k,i) = d^2 p_ci / (d q_j d q_k)
    n = 4;
    L = L(:); l = l(:);
    ph = cumsum(theta(:));
    Jv  = zeros(2, n, n);
    dJv = zeros(2, n, n, n);
    for i = 1:n
        a = [L(1:i-1); l(i)];          % distanze lungo i link fino al baricentro i
        for j = 1:i
            for r = j:i
                Jv(:,j,i) = Jv(:,j,i) + a(r)*[-sin(ph(r)); cos(ph(r))];
            end
            for k = 1:n
                for r = max(j,k):i
                    dJv(:,j,k,i) = dJv(:,j,k,i) - a(r)*[cos(ph(r)); sin(ph(r))];
                end
            end
        end
    end
end

%% Matrice di inerzia M e sue derivate dM(:,:,r) = dM/dq_r
function [M, dM] = inertia_4R(theta, L, m, l, I)
    n = 4;
    [Jv, dJv] = jac_com_4R(theta, L, l);
    M  = zeros(n);
    dM = zeros(n, n, n);
    for i = 1:n
        Jvi = Jv(:,:,i);
        Jw  = [ones(1,i), zeros(1,n-i)];
        M = M + m(i)*(Jvi'*Jvi) + I(i)*(Jw'*Jw);
        for r = 1:n
            dJ = dJv(:,:,r,i);
            dM(:,:,r) = dM(:,:,r) + m(i)*(dJ'*Jvi + Jvi'*dJ);
        end
    end
end

%% Dinamica - Matrice di inerzia M
function M = M_4R(theta, L, m, l, I)
    [M, ~] = inertia_4R(theta, L, m, l, I);
end

%% Dinamica - Vettore centrifugo/Coriolis C(q,dq)
function C = C_4R(theta, dtheta, L, m, l, I)
    n = 4;
    [~, dM] = inertia_4R(theta, L, m, l, I);
    dtheta = dtheta(:);
    C = zeros(n, 1);
    for k = 1:n
        for j = 1:n
            for r = 1:n
                C(k) = C(k) + 0.5*(dM(k,j,r) + dM(k,r,j) - dM(j,r,k))*dtheta(j)*dtheta(r);
            end
        end
    end
end

%% Dinamica - Vettore di gravita' G(q)
function G = g_4R(theta, L, m, l, g)
    n = 4;
    Jv = jac_com_4R(theta, L, l);
    G = zeros(n, 1);
    for i = 1:n
        G = G + m(i)*g*Jv(2,:,i)';
    end
end

%% Energia potenziale (gravita' lungo -y)
function V = energia_potenziale(theta, L, m, l, g)
    ph = cumsum(theta(:));
    V = 0; yj = 0;
    for i = 1:4
        V  = V + m(i)*g*(yj + l(i)*sin(ph(i)));
        yj = yj + L(i)*sin(ph(i));
    end
end

%% IK analitica 4R planare (theta4 fissato per la ridondanza)
function theta = ik_4R_planar(x, y, phi, L, elbow, theta4)
    L1 = L(1); L2 = L(2); L3 = L(3); L4 = L(4);

    % polso e fine del link 2
    wx = x - L4*cos(phi);
    wy = y - L4*sin(phi);
    phi3 = phi - theta4;
    x2 = wx - L3*cos(phi3);
    y2 = wy - L3*sin(phi3);

    D = (x2^2 + y2^2 - L1^2 - L2^2)/(2*L1*L2);
    if abs(D) > 1
        theta = NaN(4,1);
        return;
    end

    t2 = atan2(elbow*sqrt(1-D^2), D);
    t1 = atan2(y2, x2) - atan2(L2*sin(t2), L1 + L2*cos(t2));
    t3 = phi3 - (t1 + t2);
    t4 = phi - (t1 + t2 + t3);

    theta = [t1; t2; t3; t4];
end

%% FK DH 4R
function [T, x, y, phi, Tall] = fk_4R_DH(theta, L)
    Tall = zeros(4, 4, 5);
    Tall(:,:,1) = eye(4);
    T = eye(4);
    for i = 1:4
        T = T * dh_matrix(L(i), 0, 0, theta(i));
        Tall(:,:,i+1) = T;
    end
    x = T(1,4);
    y = T(2,4);
    phi = atan2(T(2,1), T(1,1));
end

%% Jacobiano 4R (righe: x, y, phi)
function J = jacobian_4R(theta, L)
    ph = cumsum(theta(:));
    J = zeros(3,4);
    for j = 1:4
        J(1,j) = -sum(L(j:4).*sin(ph(j:4)));
        J(2,j) =  sum(L(j:4).*cos(ph(j:4)));
    end
    J(3,:) = 1;
end

%% Matrice DH standard
function A = dh_matrix(a, alpha, d, theta)
    A = [ cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta);
          sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
          0,           sin(alpha),            cos(alpha),             d;
          0,           0,                     0,                      1 ];
end

%% Plot del robot
function plot_4R(theta, L)
    [~, ~, ~, ~, Tall] = fk_4R_DH(theta, L);
    P = squeeze(Tall(1:2,4,:));     % 2x5, origini O0..O4

    plot(P(1,:), P(2,:), 'ko-', 'LineWidth', 2, 'MarkerSize', 5);
    plot(P(1,end), P(2,end), 'ms', 'MarkerFaceColor', 'm', 'MarkerSize', 7);
end
