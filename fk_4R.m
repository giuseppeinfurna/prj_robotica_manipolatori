function [T, x, y, phi, Tall] = fk_4R(theta, L)
% Cinematica diretta del 4R planare con DH standard (d = 0, alpha = 0, a_i = L_i)
% T = A1*A2*A3*A4, phi = somma dei giunti, Tall = terne T0..T4 (4x4x5) per il disegno

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
