function plot_spectra(lam_pri, lam_pos, gamma_pri, gamma_pos, cfg, out_base)
w = 15.4;
h = 5.8;
FS = 7;

blue   = [0.00 0.45 0.74];
vermil = [0.85 0.33 0.10];
black  = [0.08 0.08 0.08];
red    = [0.75 0.10 0.10];
lgrey  = [0.80 0.80 0.80];
lverm  = [0.95 0.78 0.65];

fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 w h], ...
             'Color', 'w', 'Renderer', 'painters');

pos = [0.080 0.22 0.245 0.62;
       0.405 0.22 0.245 0.62;
       0.730 0.22 0.245 0.62];
ax = gobjects(1, 3);
for k = 1:3
    ax(k) = axes('Parent', fig, 'Units', 'normalized', 'Position', pos(k,:));
end

m = cfg.m;
ylo = 1e-3;
yhi = 1e4;

semilogy(ax(1), 1:m, lam_pri, '-', 'Color', lgrey, 'LineWidth', 0.4); hold(ax(1), 'on');
semilogy(ax(1), 1:m, mean(lam_pri, 2), 'o-', 'Color', black, 'LineWidth', 1.2, ...
         'MarkerSize', 3.5, 'MarkerFaceColor', black);
semilogy(ax(1), [1 m], [cfg.tau_loc cfg.tau_loc], '--', 'Color', red, 'LineWidth', 0.9);
title(ax(1), '(a) Prior samples, local', 'FontWeight', 'normal', 'FontSize', FS);
ylabel(ax(1), '$\lambda_i$', 'Interpreter', 'latex', 'FontSize', FS + 1);

semilogy(ax(2), 1:m, lam_pos, '-', 'Color', lverm, 'LineWidth', 0.4); hold(ax(2), 'on');
semilogy(ax(2), 1:m, mean(lam_pos, 2), 'o-', 'Color', vermil, 'LineWidth', 1.2, ...
         'MarkerSize', 3.5, 'MarkerFaceColor', vermil);
semilogy(ax(2), [1 m], [cfg.tau_loc cfg.tau_loc], '--', 'Color', red, 'LineWidth', 0.9);
title(ax(2), '(b) Posterior samples, local', 'FontWeight', 'normal', 'FontSize', FS);

n_plot = 20;
h1 = semilogy(ax(3), 1:n_plot, gamma_pri(1:n_plot), 'o-', 'Color', black, 'LineWidth', 1.2, ...
              'MarkerSize', 3.5, 'MarkerFaceColor', black); hold(ax(3), 'on');
h2 = semilogy(ax(3), 1:n_plot, gamma_pos(1:n_plot), 's-', 'Color', vermil, 'LineWidth', 1.2, ...
              'MarkerSize', 3.5, 'MarkerFaceColor', vermil);
h3 = semilogy(ax(3), [1 n_plot], [cfg.tau_g cfg.tau_g], '--', 'Color', red, 'LineWidth', 0.9);
h4 = semilogy(ax(3), [m m], [ylo yhi], ':', 'Color', blue, 'LineWidth', 1.2);
title(ax(3), '(c) Global', 'FontWeight', 'normal', 'FontSize', FS);
ylabel(ax(3), '$\gamma_j$', 'Interpreter', 'latex', 'FontSize', FS + 1);

for k = 1:3
    set(ax(k), 'FontSize', FS, 'Box', 'off', 'LineWidth', 0.6, 'TickDir', 'out', ...
        'YGrid', 'on', 'XGrid', 'on', 'GridAlpha', 0.12, 'YLim', [ylo yhi], ...
        'YTick', 10.^(-2:2:4), 'Layer', 'top');
    xlabel(ax(k), 'index', 'FontSize', FS);
end
set(ax(1), 'XLim', [1 m], 'XTick', 1:m);
set(ax(2), 'XLim', [1 m], 'XTick', 1:m, 'YTickLabel', []);
set(ax(3), 'XLim', [1 n_plot], 'XTick', [1 5 10 15 20]);

lg = legend(ax(3), [h1 h2 h3 h4], {'prior', 'posterior', '$\tau_g$', '$j=m$'}, ...
            'Interpreter', 'latex', 'FontSize', FS - 0.5, 'Location', 'northeast', 'Box', 'off');
lg.ItemTokenSize = [14 8];

set(fig, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 w h], 'PaperSize', [w h]);
drawnow;
print(fig, [out_base '.pdf'], '-dpdf', '-painters');
exportgraphics(fig, [out_base '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
close(fig);
end
