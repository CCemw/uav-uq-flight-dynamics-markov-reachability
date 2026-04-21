% PLOT_STABILITY_RESULTS
% Produces:
%   figs/Fig07_eigenvalue_map.pdf       nominal eigenvalue map in the
%                                       complex plane, key speeds highlighted
%
%   figs/Fig08_Re_lambda_t2.pdf         max Re(lambda) and t_double vs airspeed
%
%   figs/Fig09a_eig_clouds_low.pdf      MC eigenvalue clouds at V=5 and V=10
%
%   figs/Fig09b_eig_clouds_high.pdf     MC eigenvalue clouds at V=15 and V=24
%
%   console table                       MC max Re(lambda), M_u, M_q at all
%                                       target speeds
%
% Reads:  results/stability_results_QAV250.mat

clear; close all; clc;
if ~isfolder('figs'); mkdir('figs'); end   % make the figs folder if it doesn't exist yet

% Style ---------------------------------------------------------------------
set(0, 'DefaultAxesFontName',   'Times New Roman');
set(0, 'DefaultTextFontName',   'Times New Roman');
set(0, 'DefaultAxesFontSize',   13);
set(0, 'DefaultTextFontSize',   13);
set(0, 'DefaultLineLineWidth',  2.0);
set(0, 'DefaultAxesLineWidth',  1.2);
set(0, 'DefaultAxesBox',        'on');
set(0, 'DefaultAxesTickDir',    'in');
set(0, 'DefaultAxesTickLength', [0.015 0.015]);
set(0, 'DefaultFigureColor',    'w');

col_nom    = [0.00 0.00 0.00];
col_mc     = [0.00 0.45 0.70];
col_grey   = [0.50 0.50 0.50];
col_red    = [0.80 0.20 0.20];
col_green  = [0.00 0.62 0.45];
col_amber  = [0.90 0.60 0.00];

% Load ----------------------------------------------------------------------
load('stability_results_QAV250.mat');
V = V_sweep(:);

maxRe_nom = max(eig_real_nom, [], 2);
t2_nom    = log(2) ./ max(maxRe_nom, 1e-10);
t2_nom(maxRe_nom <= 0) = Inf;

% Fig 7: nominal eigenvalue map ---------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 13 9]);

key_speeds  = [5, 10, 15, 24];
key_markers = {'s', 'd', '^', 'v'};
key_colors  = {col_green, col_mc, col_amber, col_red};

for vi = 1:nV
    if ~nom_conv(vi), continue; end
    plot(eig_real_nom(vi,:), eig_imag_nom(vi,:), '.', ...
        'Color', [0.75 0.75 0.75], 'MarkerSize', 6, 'HandleVisibility', 'off');
    hold on;
end
for ki = 1:numel(key_speeds)
    [~, vi] = min(abs(V_sweep - key_speeds(ki)));
    if ~nom_conv(vi), continue; end
    plot(eig_real_nom(vi,:), eig_imag_nom(vi,:), key_markers{ki}, ...
        'Color', key_colors{ki}, 'MarkerFaceColor', key_colors{ki}, ...
        'MarkerSize', 7, 'LineWidth', 1.0, ...
        'DisplayName', sprintf('V = %d m/s', key_speeds(ki)));
end
xline(0, '-', 'Color', 'k', 'LineWidth', 0.6, 'HandleVisibility', 'off');
yline(0, '-', 'Color', 'k', 'LineWidth', 0.6, 'HandleVisibility', 'off');
hold off;
xlabel('Re(\lambda) [1/s]'); ylabel('Im(\lambda) [rad/s]');
grid on;
legend('Location', 'southwest', 'FontSize', 8, 'Box', 'off');

drawnow;
exportgraphics(gcf, 'figs/Fig07_eigenvalue_map.pdf', 'ContentType', 'vector');

% Fig 8: Re(lambda) and t_double vs airspeed --------------------------------
valid = nom_conv & (1:nV)' > 1;

figure('Units', 'centimeters', 'Position', [2 2 15.5 7.5]);
sub1 = subplot(1,2,1);
plot(V(valid), maxRe_nom(valid), '-o', 'Color', col_nom, ...
    'MarkerSize', 4, 'MarkerFaceColor', col_nom, 'LineWidth', 1.8);
hold on;
yline(0, '--', 'Color', col_grey, 'LineWidth', 0.8);
hold off;
xlabel('Airspeed [m/s]'); ylabel('max Re(\lambda) [1/s]');
xlim([0 V(end)]); grid on;

sub2 = subplot(1,2,2);
plot(V(valid), t2_nom(valid), '-o', 'Color', col_nom, ...
    'MarkerSize', 4, 'MarkerFaceColor', col_nom, 'LineWidth', 1.8);
xlabel('Airspeed [m/s]'); ylabel('t_{double} [s]');
xlim([0 V(end)]); ylim([0 1.2]); grid on;

w = 0.35; gap = 0.12; left1 = 0.10; left2 = left1 + w + gap;
p1 = get(sub1, 'Position'); p2 = get(sub2, 'Position');
set(sub1, 'Position', [left1, p1(2), w, p1(4)]);
set(sub2, 'Position', [left2, p2(2), w, p2(4)]);

drawnow;
exportgraphics(gcf, 'figs/Fig08_Re_lambda_t2.pdf', 'ContentType', 'vector');

% Fig 9a / 9b: MC eigenvalue clouds -----------------------------------------
cloud_idx = [1, 2, 3, 6];                       % into V_mc_targets
cloud_V   = V_mc_targets(cloud_idx);
cloud_col = [col_green; col_mc; col_amber; col_red];

plot_eig_clouds(cloud_idx(1:2), cloud_V(1:2), cloud_col(1:2,:), ...
                V, V_mc_targets, mc_eig_conv, mc_eig_real, mc_eig_imag, ...
                eig_real_nom, eig_imag_nom, nom_conv, N_MC, ...
                'figs/Fig09a_eig_clouds_low.pdf');

plot_eig_clouds(cloud_idx(3:4), cloud_V(3:4), cloud_col(3:4,:), ...
                V, V_mc_targets, mc_eig_conv, mc_eig_real, mc_eig_imag, ...
                eig_real_nom, eig_imag_nom, nom_conv, N_MC, ...
                'figs/Fig09b_eig_clouds_high.pdf');

% Console table -------------------------------------------------------------
fprintf('\nMonte Carlo stability summary at target airspeeds\n');
fprintf('%-6s %-9s %-22s %-22s %-22s\n', ...
        'V', 'n/N', 'max Re(lambda)', 'M_u', 'M_q');
fprintf('%s\n', repmat('-', 1, 84));
for ti = 1:numel(V_mc_targets)
    conv_mask = logical(mc_eig_conv(ti,:));
    nc = sum(conv_mask);
    re = mc_eig_maxRe(ti, conv_mask);
    mu = mc_deriv_Mu(ti, conv_mask);
    mq = mc_deriv_Mq(ti, conv_mask);
    fprintf('%-6d %3d/%-5d %+0.3f +/- %-10.3f %+0.3f +/- %-10.3f %+0.3f +/- %-10.3f\n', ...
            V_mc_targets(ti), nc, N_MC, ...
            mean(re), std(re), mean(mu), std(mu), mean(mq), std(mq));
end

fprintf('\nSaved: figs/Fig07_eigenvalue_map.pdf, Fig08_Re_lambda_t2.pdf, Fig09a/b_eig_clouds.pdf\n');


%  LOCAL FUNCTIONS ===========================================================================

function plot_eig_clouds(idx_pair, V_pair, col_pair, V, V_mc_targets, ...
                          mc_eig_conv, mc_eig_real, mc_eig_imag, ...
                          eig_real_nom, eig_imag_nom, nom_conv, N_MC, out_path)
% Two-panel MC eigenvalue cloud, one panel per target speed.
    figure('Units', 'centimeters', 'Position', [2 2 15.5 7.5]);
    for pp = 1:2
        ti = idx_pair(pp);
        subplot(1, 2, pp);

        conv_mask = logical(mc_eig_conv(ti,:));
        n_conv    = sum(conv_mask);
        re        = squeeze(mc_eig_real(ti, conv_mask, :));
        im        = squeeze(mc_eig_imag(ti, conv_mask, :));

        plot(re(:), im(:), '.', 'Color', [col_pair(pp,:) 0.15], ...
            'MarkerSize', 3, 'HandleVisibility', 'off');
        hold on;
        plot(NaN, NaN, '.', 'Color', col_pair(pp,:), 'MarkerSize', 10, ...
            'DisplayName', 'MC samples');

        [~, nom_idx] = min(abs(V - V_pair(pp)));
        if nom_conv(nom_idx)
            plot(eig_real_nom(nom_idx,:), eig_imag_nom(nom_idx,:), 'k+', ...
                'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Nominal');
        end
        xline(0, '-', 'Color', 'k', 'LineWidth', 0.6, 'HandleVisibility', 'off');
        yline(0, '-', 'Color', 'k', 'LineWidth', 0.6, 'HandleVisibility', 'off');

        xl = xlim; yl = ylim;
        text(xl(1) + 0.05*(xl(2)-xl(1)), yl(2) - 0.05*(yl(2)-yl(1)), ...
            sprintf('%d/%d', n_conv, N_MC), ...
            'FontSize', 9, 'FontWeight', 'bold', 'Color', col_pair(pp,:), ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

        xlabel('Re(\lambda) [1/s]');
        if pp == 1; ylabel('Im(\lambda) [rad/s]'); end
        legend('Location', 'southwest', 'FontSize', 9, 'Box', 'off');
        grid on; hold off;
    end

    sub1 = subplot(1,2,1); sub2 = subplot(1,2,2);
    w = 0.38; h = 0.72; gap = 0.12; left1 = 0.10; left2 = left1 + w + gap;
    bot = 0.18;
    set(sub1, 'Position', [left1, bot, w, h]);
    set(sub2, 'Position', [left2, bot, w, h]);

    drawnow;
    exportgraphics(gcf, out_path, 'ContentType', 'vector');
end