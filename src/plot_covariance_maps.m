function plot_covariance_maps(G, nt, out_cov, out_err)
FS = 7;
w  = 15.4;

M = {G.pr, G.ref, G.lis, G.lisgp, G.alone};
T = {'(a) Prior', '(b) Reference', '(c) LIS alone', ...
     sprintf('(d) LIS+GP, $n_\\mathrm{tr}=%d$', nt), sprintf('(e) GP alone, $n_\\mathrm{tr}=%d$', nt)};
cmax = max(G.pr(:));
h = 4.3;
fig = new_fig(w, h);
x0 = 0.055; dx = 0.178; wd = 0.155;
for k = 1:5
    ax = axes('Parent', fig, 'Position', [x0 + (k-1)*dx, 0.17, wd, 0.66]);
    style(ax, M{k}, [0 cmax], T{k}, FS, k == 1);
end
cb = colorbar(ax, 'Position', [x0 + 4*dx + wd + 0.012, 0.17, 0.011, 0.66]);
cb.FontSize = FS - 0.5;
cb.Label.String = 'covariance';
export(fig, w, h, out_cov);

U   = abs(G.pr - G.ref);
nrm = max(U(:));
E   = {abs(G.lis - G.ref)/nrm, abs(G.lisgp - G.ref)/nrm, abs(G.alone - G.ref)/nrm};
T   = {'(a) LIS alone', sprintf('(b) LIS+GP, $n_\\mathrm{tr}=%d$', nt), ...
       sprintf('(c) GP alone, $n_\\mathrm{tr}=%d$', nt)};
emax = max(cellfun(@(e) max(e(:)), E));
h = 5.6;
fig = new_fig(w, h);
x0 = 0.07; dx = 0.29; wd = 0.25;
for k = 1:3
    ax = axes('Parent', fig, 'Position', [x0 + (k-1)*dx, 0.15, wd, 0.72]);
    style(ax, E{k}, [0 emax], T{k}, FS, k == 1);
end
cb = colorbar(ax, 'Position', [x0 + 2*dx + wd + 0.015, 0.15, 0.013, 0.72]);
cb.FontSize = FS - 0.5;
cb.Label.String = '$|\hat\Gamma-\Gamma_\mathrm{ref}|\,/\,\max|\Gamma_\mathrm{pr}-\Gamma_\mathrm{ref}|$';
cb.Label.Interpreter = 'latex';
export(fig, w, h, out_err);
end

function fig = new_fig(w, h)
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 w h], ...
             'Color', 'w', 'Renderer', 'painters');
end

function style(ax, A, clim, ttl, FS, first)
imagesc(ax, A, clim);
axis(ax, 'square');
colormap(ax, parula);
title(ax, ttl, 'Interpreter', 'latex', 'FontSize', FS);
set(ax, 'FontSize', FS, 'TickDir', 'out', 'XTick', [1 25 50], 'YTick', [1 25 50], 'LineWidth', 0.5);
xlabel(ax, 'node $k$', 'Interpreter', 'latex', 'FontSize', FS + 1);
if first
    ylabel(ax, 'node $j$', 'Interpreter', 'latex', 'FontSize', FS + 1);
else
    set(ax, 'YTickLabel', []);
end
end

function export(fig, w, h, out_base)
set(fig, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 w h], 'PaperSize', [w h]);
drawnow;
print(fig, [out_base '.pdf'], '-dpdf', '-painters');
exportgraphics(fig, [out_base '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
close(fig);
end
