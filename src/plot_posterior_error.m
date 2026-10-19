function plot_posterior_error(E, cfg, r, out_base)
w = 15.4;
h = 6.0;
FS = 7;
x = cfg.ntr_list;

vermil = [0.85 0.33 0.10];
black  = [0.08 0.08 0.08];

fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 w h], ...
             'Color', 'w', 'Renderer', 'painters');
pos = [0.085 0.20 0.395 0.68;
       0.585 0.20 0.395 0.68];

A   = {E.mean_alone, E.cov_alone};
B   = {E.mean_lisgp, E.cov_lisgp};
fa  = [E.mean_full, E.cov_full];
fl  = [E.mean_lis,  E.cov_lis];
ttl = {'(a) Posterior mean', '(b) Posterior covariance'};
yl  = {'$e_\mu$', '$e_\Gamma$'};

for k = 1:2
    a = axes('Parent', fig, 'Units', 'normalized', 'Position', pos(k,:));
    hold(a, 'on');
    set(a, 'XScale', 'log', 'YScale', 'log');
    lo = min(A{k}, [], 1);  hi = max(A{k}, [], 1);
    fill(a, [x fliplr(x)], [lo fliplr(hi)], black, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    lo = min(B{k}, [], 1);  hi = max(B{k}, [], 1);
    fill(a, [x fliplr(x)], [lo fliplr(hi)], vermil, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    h1 = plot(a, x, median(A{k}, 1), 'o-', 'Color', black, 'LineWidth', 1.2, ...
              'MarkerSize', 3.5, 'MarkerFaceColor', black);
    h2 = plot(a, x, median(B{k}, 1), 's-', 'Color', vermil, 'LineWidth', 1.2, ...
              'MarkerSize', 3.5, 'MarkerFaceColor', vermil);
    h3 = plot(a, [x(1) x(end)], [fa(k) fa(k)], ':', 'Color', black, 'LineWidth', 1.0);
    h4 = plot(a, [x(1) x(end)], [fl(k) fl(k)], ':', 'Color', vermil, 'LineWidth', 1.0);
    allv = [A{k}(:); B{k}(:); fa(k); fl(k)];
    set(a, 'FontSize', FS, 'Box', 'off', 'LineWidth', 0.6, 'TickDir', 'out', ...
        'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.12, 'Layer', 'top', ...
        'XLim', [x(1)*0.9 x(end)*1.1], 'YLim', [min(allv)*0.8 max(allv)*1.25], ...
        'XTick', x(1:2:end), 'XTickLabel', arrayfun(@num2str, x(1:2:end), 'UniformOutput', false));
    title(a, ttl{k}, 'FontWeight', 'normal', 'FontSize', FS);
    xlabel(a, '$n_\mathrm{tr}$', 'Interpreter', 'latex', 'FontSize', FS + 1);
    ylabel(a, yl{k}, 'Interpreter', 'latex', 'FontSize', FS + 1);
    if k == 1
        lg = legend(a, [h1 h2 h3 h4], {'GP alone', sprintf('LIS+GP, $r=%d$', r), ...
                    'full space, exact model', sprintf('LIS alone, $r=%d$', r)}, ...
                    'Interpreter', 'latex', 'FontSize', FS - 0.5, 'Location', 'northeast', 'Box', 'off');
        lg.ItemTokenSize = [14 8];
    end
end

set(fig, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 w h], 'PaperSize', [w h]);
drawnow;
print(fig, [out_base '.pdf'], '-dpdf', '-painters');
exportgraphics(fig, [out_base '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
close(fig);
end
