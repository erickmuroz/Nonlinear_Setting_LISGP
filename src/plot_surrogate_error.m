function plot_surrogate_error(res, cfg, out_base)
w = 15.4;
h = 6.6;
FS = 7;
x = cfg.ntr_list;

cols = [0.35 0.62 0.88; 0.85 0.33 0.10; 0.00 0.45 0.74];
black = [0.08 0.08 0.08];
grey  = [0.50 0.50 0.50];

fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 w h], ...
             'Color', 'w', 'Renderer', 'painters');
ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.085 0.19 0.895 0.76]);
set(ax, 'XScale', 'log', 'YScale', 'log');
hold(ax, 'on');

lo = min(res.alone, [], 1);  hi = max(res.alone, [], 1);
fill(ax, [x fliplr(x)], [lo fliplr(hi)], black, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
hl = gobjects(1, 1 + numel(cfg.r_list));
hl(1) = plot(ax, x, median(res.alone, 1), 'o-', 'Color', black, 'LineWidth', 1.2, ...
             'MarkerSize', 3.5, 'MarkerFaceColor', black);

names = {'GP alone'};
for a = 1:numel(cfg.r_list)
    E  = res.lis(:,:,a);
    lo = min(E, [], 1);  hi = max(E, [], 1);
    fill(ax, [x fliplr(x)], [lo fliplr(hi)], cols(a,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    hl(a+1) = plot(ax, x, median(E, 1), 's-', 'Color', cols(a,:), 'LineWidth', 1.2, ...
                   'MarkerSize', 3.5, 'MarkerFaceColor', cols(a,:));
    plot(ax, [x(1) x(end)], [res.floor(a) res.floor(a)], ':', 'Color', cols(a,:), 'LineWidth', 1.0);
    names{end+1} = sprintf('LIS+GP, $r=%d$', cfg.r_list(a));
end

plot(ax, [x(1) x(end)], [1 1], '--', 'Color', grey, 'LineWidth', 0.8);
text(ax, x(1)*1.05, 1.12, '$\sigma_\mathrm{obs}$', 'Interpreter', 'latex', 'Color', grey, 'FontSize', FS);

a  = find(cfg.r_list == res.r_main, 1);
n1 = budget_at(x, median(res.lis(:,:,a), 1), 1);
n0 = budget_at(x, median(res.alone, 1), 1);
if ~isnan(n0) && ~isnan(n1)
    plot(ax, [n1 n0], [1 1], '-', 'Color', cols(a,:), 'LineWidth', 1.6);
    plot(ax, [n1 n1], [0.93 1.07], '-', 'Color', cols(a,:), 'LineWidth', 1.0);
    plot(ax, [n0 n0], [0.93 1.07], '-', 'Color', cols(a,:), 'LineWidth', 1.0);
    text(ax, sqrt(n1*n0), 1.09, sprintf('$\\times %.1f$', n0/n1), 'Interpreter', 'latex', ...
         'Color', cols(a,:), 'FontSize', FS, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

allv = [res.alone(:); res.lis(:); res.floor(:); 1];
ylim_ = [min(allv)*0.9 max(allv)*1.15];
set(ax, 'FontSize', FS, 'Box', 'off', 'LineWidth', 0.6, 'TickDir', 'out', ...
    'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.12, 'Layer', 'top', ...
    'XLim', [x(1)*0.9 x(end)*1.1], 'YLim', ylim_, 'YTick', [0.5 1 2 4 8], 'YTickLabel', {'0.5','1','2','4','8'}, 'XTick', x, 'XTickLabel', arrayfun(@num2str, x, 'UniformOutput', false));
xlabel(ax, '$n_\mathrm{tr}$', 'Interpreter', 'latex', 'FontSize', FS + 1);
ylabel(ax, '$\mathrm{RMSE}(\hat G)\,/\,\sigma_\mathrm{obs}$', 'Interpreter', 'latex', 'FontSize', FS + 1);

lg = legend(ax, hl, names, 'Interpreter', 'latex', 'FontSize', FS - 0.5, ...
            'Location', 'northeast', 'Box', 'off');
lg.ItemTokenSize = [14 8];

set(fig, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 w h], 'PaperSize', [w h]);
drawnow;
print(fig, [out_base '.pdf'], '-dpdf', '-painters');
exportgraphics(fig, [out_base '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
close(fig);
end
