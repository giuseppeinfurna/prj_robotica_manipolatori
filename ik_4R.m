function theta = ik_4R(x, y, phi, L, opts)
% IK_4R  Cinematica inversa analitica del manipolatore planare 4R.
%
%   theta = ik_4R(x, y, phi, L)
%   theta = ik_4R(x, y, phi, L, opts)
%
% Con 3 coordinate di task (x, y, phi) e 4 giunti il robot e' ridondante di
% grado 1: va fissato un parametro libero. Qui si fissa theta4:
%   opts.theta4 (default 0.3 rad) : valore desiderato di theta4
%   opts.elbow  (default +1)      : +1 gomito basso, -1 gomito alto
%
% Procedura:
%   1) polso (punto a monte del link 4):  W = [x - L4 cos(phi); y - L4 sin(phi)]
%   2) orientamento del link 3:           phi3 = phi - theta4
%   3) fine link 2:                       P2 = W - L3 [cos(phi3); sin(phi3)]
%   4) problema 2R per theta1, theta2 (legge dei coseni, 2 soluzioni)
%   5) theta3 = phi3 - theta1 - theta2,   theta4 = phi - (theta1+theta2+theta3)
%
% Se il punto non e' raggiungibile solleva l'errore 'ik_4R:Unreachable'.

if nargin < 5, opts = struct(); end
if ~isfield(opts, 'theta4'), opts.theta4 = 0.3; end
if ~isfield(opts, 'elbow'),  opts.elbow = 1;   end

L1 = L(1); L2 = L(2); L3 = L(3); L4 = L(4);

xw = x - L4*cos(phi);
yw = y - L4*sin(phi);

phi3 = phi - opts.theta4;
x2 = xw - L3*cos(phi3);
y2 = yw - L3*sin(phi3);

c2 = (x2^2 + y2^2 - L1^2 - L2^2) / (2*L1*L2);
if abs(c2) > 1
    error('ik_4R:Unreachable', ...
        'Punto (%.3f, %.3f, phi = %.3f) fuori dallo spazio di lavoro raggiungibile.', x, y, phi);
end

s2 = opts.elbow * sqrt(1 - c2^2);
th2 = atan2(s2, c2);
th1 = atan2(y2, x2) - atan2(L2*s2, L1 + L2*c2);
th3 = phi3 - th1 - th2;
th4 = phi - (th1 + th2 + th3);

theta = [th1; th2; th3; th4];
end
