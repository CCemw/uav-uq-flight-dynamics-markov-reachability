% PLOT_MANOEUVRE_RESPONSE
% Produces:
%   figs/Fig10_manoeuvre_V<V_plot>.pdf   pitch deviation and pitch rate vs
%                                        time, with MC 5th-95th and 25th-75th
%                                        percentile bands and the nominal
%                                        run overlaid
%   console summary                      converged / diverged counts and
%                                        peak pitch deviation statistics
%
% USER SETTINGS at the top: change V_plot to switch between cases.
%
% Reads:  results/manoeuvre_V<V_plot>_nominal_QAV250.mat
%         results/manoeuvre_V<V_plot>_mc_QAV250.mat

clear; close all; clc;
if ~isfolder('figs'); mkdir('figs'); end   % make the figs folder if it doesn't exist yet

% USER SETTINGS -------------------------------------------------------------
V_plot = 15;     % airspeed of the case to plot [m/s] >>> !! MUST match a saved pair !!
% --------------------------------------------------------------------------

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
col_shade  = [0.50 0.72 0.88];
col_shade2 = [0.75 0.86 0.94];
col_grey   = [0.50 0.50 0.50];
col_amber  = [0.85 0.65 0.00];

% Load ----------------------------------------------------------------------
nom_file = sprintf('results/manoeuvre_V%d_nominal_QAV250.mat', V_plot);
mc_file  = sprintf('results/manoeuvre_V%d_mc_QAV250.mat',      V_plot);
nom = load(nom_file);
mc  = load(mc_file);

t_nom = nom.t_vec;
t_mc  = mc.t_vec;

% Pitch deviation from each sample's own trim
dtheta_mc  = mc.all_theta  - mc.mc_trim_theta;
dtheta_nom = nom.all_theta(1,:) - nom.all_theta(1,1);

% Pitch rate (already absolute, no trim subtraction)
q_mc  = mc.all_q;
q_nom = nom.all_q(1,:);

% Bands across all converged samples ---------------------------------------
v_mask = mc.run_conv;

th_p5  = prctile(dtheta_mc(v_mask,:), 5,  1);
th_p25 = prctile(dtheta_mc(v_mask,:), 25, 1);
th_p50 = prctile(dtheta_mc(v_mask,:), 50, 1);
th_p75 = prctile(dtheta_mc(v_mask,:), 75, 1);
th_p95 = prctile(dtheta_mc(v_mask,:), 95, 1);

q_p5  = prctile(q_mc(v_mask,:), 5,  1);
q_p25 = prctile(q_mc(v_mask,:), 25, 1);
q_p50 = prctile(q_mc(v_mask,:), 50, 1);
q_p75 = prctile(q_mc(v_mask,:), 75, 1);
q_p95 = prctile(q_mc(v_mask,:), 95, 1);

n_conv = sum(v_mask);
n_div  = sum(mc.diverged);

% Two-panel figure ----------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 15.5 7.5]);

% (a) Pitch deviation
sub1 = subplot(1,2,1);
hold on;
fill([t_mc fliplr(t_mc)], [th_p5 fliplr(th_p95)], col_shade2, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', '5th-95th pctl.');
fill([t_mc fliplr(t_mc)], [th_p25 fliplr(th_p75)], col_shade, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', '25th-75th pctl.');
plot(t_mc, th_p50, '-', 'Color', col_mc, 'LineWidth', 1.5, ...
    'DisplayName', 'Median');
plot(t_nom, dtheta_nom, '-', 'Color', col_nom, 'LineWidth', 2.0, ...
    'DisplayName', 'Nominal');
yline(0, '--', 'Color', col_grey, 'LineWidth', 0.8, ...
    'DisplayName', 'Trim ref.');
xline(mc.t_pulse, ':', 'Color', col_amber, 'LineWidth', 1.2, ...
    'DisplayName', 'Disturbance');
hold off;
xlabel('Time [s]'); ylabel('\Delta\theta [deg]');
xlim([0 mc.t_sim]); grid on;
legend('Location', 'southwest', 'FontSize', 9, 'Box', 'off');

% (b) Pitch rate
sub2 = subplot(1,2,2);
hold on;
fill([t_mc fliplr(t_mc)], [q_p5 fliplr(q_p95)], col_shade2, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', '5th-95th pctl.');
fill([t_mc fliplr(t_mc)], [q_p25 fliplr(q_p75)], col_shade, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', '25th-75th pctl.');
plot(t_mc, q_p50, '-', 'Color', col_mc, 'LineWidth', 1.5, ...
    'DisplayName', 'Median');
plot(t_nom, q_nom, '-', 'Color', col_nom, 'LineWidth', 2.0, ...
    'DisplayName', 'Nominal');
yline(0, '--', 'Color', col_grey, 'LineWidth', 0.8, ...
    'DisplayName', 'Trim ref.');
xline(mc.t_pulse, ':', 'Color', col_amber, 'LineWidth', 1.2, ...
    'DisplayName', 'Disturbance');
hold off;
xlabel('Time [s]'); ylabel('q [deg/s]');
xlim([0 mc.t_sim]); grid on;
legend('Location', 'southwest', 'FontSize', 9, 'Box', 'off');

w = 0.38; h = 0.72; gap = 0.12; left1 = 0.08; left2 = left1 + w + gap;
bot = 0.18;
set(sub1, 'Position', [left1, bot, w, h]);
set(sub2, 'Position', [left2, bot, w, h]);

drawnow;
out_pdf = sprintf('figs/Fig10_manoeuvre_V%d.pdf', V_plot);
exportgraphics(gcf, out_pdf, 'ContentType', 'vector');

% Console summary -----------------------------------------------------------
mc_peaks = max(abs(dtheta_mc(v_mask,:)), [], 2);

fprintf('\nClosed-loop manoeuvre summary at V = %d m/s\n', V_plot);
fprintf('Pulse: +%d RPM x %.1f s on rear rotors at t = %.2f s\n', ...
        mc.dRPM_pulse, mc.t_pulse_dur, mc.t_pulse);
fprintf('%s\n', repmat('-', 1, 60));
fprintf('Converged runs:   %d / %d\n', n_conv, mc.N_MC);
fprintf('Diverged runs:    %d  (|theta| exceeded %.0f deg)\n', n_div, mc.theta_diverge);
fprintf('Nominal peak |delta theta|:  %.2f deg\n', max(abs(dtheta_nom)));
fprintf('MC peak |delta theta|:  median %.2f, 95th pctl %.2f, max %.2f deg\n', ...
        median(mc_peaks), prctile(mc_peaks, 95), max(mc_peaks));

fprintf('\nSaved: %s\n', out_pdf);