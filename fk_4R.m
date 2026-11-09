function [T, x, y, phi, Tall] = fk_4R(theta, L)
% FK_4R  Cinematica diretta del manipolatore planare 4R (DH standard).
%
% Tabella DH (alpha_i = 0, d_i = 0 per tutti i giunti):
%   i | theta_i | d_i | a_i | alpha_i
%   1 | q1      | 0   | L1  | 0
%   2 | q2      | 0   | L2  | 0
%   3 | q3      | 0   | L3  | 0
%   4 | q4      | 0   | L4  | 0
%
% Uscite:
%   T    : matrice 4x4 base -> end-effector, T = A1*A2*A3*A4
%   x, y : posizione dell'end-effector
%   phi  : orientamento theta1+theta2+theta3+theta4
%   Tall : 4x4x5, terne T0..T4 (T0 = identita') per la visualizzazione

theta = theta(:);
L = L(:);
Tall = zeros(4, 4, 5);
Tall(:, :, 1) = eye(4);
T = eye(4);
for i = 1:4
    T = T * dh_A(theta(i), 0, L(i), 0);
    Tall(:, :, i + 1) = T;
end
x = T(1, 4);
y = T(2, 4);
phi = sum(theta);
end

function A = dh_A(th, d, a, al)
% Matrice di trasformazione omogenea DH standard.
A = [cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th);
     sin(th),  cos(th)*cos(al), -cos(th)*sin(al), a*sin(th);
     0,        sin(al),          cos(al),         d;
     0,        0,                0,               1];
end
