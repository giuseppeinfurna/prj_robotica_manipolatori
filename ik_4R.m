function theta = ik_4R(x, y, phi, L, opts)
% Cinematica inversa analitica del 4R planare.
% Con (x, y, phi) e 4 giunti il robot e' ridondante, quindi si fissa theta4:
%   opts.theta4 (default 0.3 rad), opts.elbow (+1 gomito basso, -1 gomito alto)
% Si parte dal polso, si toglie il link 3 e si risolve il 2R con la legge dei coseni.
% Se il punto e' fuori dallo spazio di lavoro l'errore e' 'ik_4R:Unreachable'.

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
        'Punto (%.3f, %.3f, phi = %.3f) fuori dallo spazio di lavoro.', x, y, phi);
end

s2 = opts.elbow * sqrt(1 - c2^2);
th2 = atan2(s2, c2);
th1 = atan2(y2, x2) - atan2(L2*s2, L1 + L2*c2);
th3 = phi3 - th1 - th2;
th4 = phi - (th1 + th2 + th3);

theta = [th1; th2; th3; th4];
end
