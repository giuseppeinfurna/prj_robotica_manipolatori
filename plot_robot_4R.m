function plot_robot_4R(ax, theta, L, frameScale)
% Disegna il robot 4R con le terne DH (x in rosso, y in verde), terna 0 alla base

if nargin < 4, frameScale = 0.08; end
[~, ~, ~, ~, Tall] = fk_4R(theta, L);

hold(ax, 'on');
P = squeeze(Tall(1:2, 4, :));   % 2x5 posizioni delle origini O0..O4
plot(ax, P(1, :), P(2, :), '-', 'Color', [0.1 0.3 0.7], 'LineWidth', 4);
plot(ax, P(1, :), P(2, :), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
plot(ax, P(1, end), P(2, end), 's', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'MarkerSize', 9);

for i = 1:5
    o = Tall(1:2, 4, i);
    xa = Tall(1:2, 1, i) * frameScale;
    ya = Tall(1:2, 2, i) * frameScale;
    quiver(ax, o(1), o(2), xa(1), xa(2), 0, 'r', 'LineWidth', 1.2, 'MaxHeadSize', 0.6);
    quiver(ax, o(1), o(2), ya(1), ya(2), 0, 'g', 'LineWidth', 1.2, 'MaxHeadSize', 0.6);
    plot(ax, o(1), o(2), '.', 'Color', 'b', 'MarkerSize', 12);
    text(ax, o(1) + 0.02, o(2) - 0.03, sprintf('O_%d', i - 1), 'FontSize', 8);
end
axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'x [m]');
ylabel(ax, 'y [m]');
end
