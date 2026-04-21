% PLOT_TRIM_RESULTS
% Produces:
%   figs/Fig04_trim_sweep.pdf     body pitch and rotor speed vs airspeed,
%                                 with MC 5th-95th percentile bands
%   figs/Fig05_power_bands.pdf    total power required vs airspeed,
%                                 with MC 5th-95th percentile band
%   figs/Fig06_Ptrim_Scurve.pdf   P(trim) vs airspeed, with deterministic
%                                 step overlay and 95% / 50% thresholds
%   console table                 V_max at 95%, 50%, and deterministic
%
% Reads:  results/nominal_trim_results_QAV250.mat
%         results/mc_trim_results_QAV250.mat
% Writes: figs/Fig04_trim_sweep.pdf
%         figs/Fig05_power_bands.pdf
%         figs/Fig06_Ptrim_Scurve.pdf

clear; close all; clc;
if ~isfolder('figs'); mkdir('figs'); end % make the figs folder if it doesn't exist yet

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

col_nom        = [0.00 0.00 0.00];
col_mc         = [0.00 0.45 0.70];
col_shade      = [0.50 0.72 0.88];
col_edge       = [0.25 0.50 0.70];
col_grey       = [0.50 0.50 0.50];
col_red        = [0.80 0.20 0.20];
col_shade_rear = [0.88 0.75 0.75];
col_dark_blue  = [0.00 0.35 0.60];

% Load ----------------------------------------------------------------------
mc  = load('mc_trim_results_QAV250.mat');
nom = load('nominal_trim_results_QAV250.mat');

V_sweep   = mc.V_sweep(:);
nV        = numel(V_sweep);
N_MC      = mc.N_MC;
nom_conv  = mc.nom_conv;
nom_theta = mc.nom_theta;
nom_RPM   = mc.nom_RPM;
mc_conv   = mc.mc_conv;
mc_theta  = mc.mc_theta;
mc_RPM    = mc.mc_RPM;
mc_power  = mc.mc_power;
P_trim    = mc.P_trim;

% Interpolate fine-grid nominal power onto MC speed grid
nom_power = NaN(nV, 1);
for vi = 1:nV
    [d, idx] = min(abs(nom.V_sweep(:) - V_sweep(vi)));
    if d < 0.01 && nom.nom_conv(idx) && ~isnan(nom.nom_power(idx))
        nom_power(vi) = nom.nom_power(idx);
    end
end

% Find last airspeed with usable data
V_max_plot = max(V_sweep(nom_conv));
for vi = nV:-1:1
    if sum(mc_conv(vi,:)) >= 5
        V_max_plot = max(V_max_plot, V_sweep(vi));
        break;
    end
end

% MC percentiles ------------------------------------------------------------
mc_th_med   = NaN(nV,1); mc_th_lo   = NaN(nV,1); mc_th_hi   = NaN(nV,1);
mc_RPMf_med = NaN(nV,1); mc_RPMf_lo = NaN(nV,1); mc_RPMf_hi = NaN(nV,1);
mc_RPMr_med = NaN(nV,1); mc_RPMr_lo = NaN(nV,1); mc_RPMr_hi = NaN(nV,1);
mc_P_med    = NaN(nV,1); mc_P_lo    = NaN(nV,1); mc_P_hi    = NaN(nV,1);

for vi = 1:nV
    converged = mc_conv(vi, :);
    if sum(converged) < 5, continue; end

    th = mc_theta(vi, converged); th = th(~isnan(th));
    if ~isempty(th)
        mc_th_med(vi) = median(th);
        mc_th_lo(vi)  = prctile(th, 5);
        mc_th_hi(vi)  = prctile(th, 95);
    end

    rpf = squeeze(mc_RPM(vi, converged, 1)); rpf = rpf(~isnan(rpf));
    rpr = squeeze(mc_RPM(vi, converged, 3)); rpr = rpr(~isnan(rpr));
    if ~isempty(rpf)
        mc_RPMf_med(vi) = median(rpf);
        mc_RPMf_lo(vi)  = prctile(rpf, 5);
        mc_RPMf_hi(vi)  = prctile(rpf, 95);
        mc_RPMr_med(vi) = median(rpr);
        mc_RPMr_lo(vi)  = prctile(rpr, 5);
        mc_RPMr_hi(vi)  = prctile(rpr, 95);
    end

    pwr = mc_power(vi, converged); pwr = pwr(~isnan(pwr) & pwr > 0);
    if ~isempty(pwr)
        mc_P_med(vi) = median(pwr);
        mc_P_lo(vi)  = prctile(pwr, 5);
        mc_P_hi(vi)  = prctile(pwr, 95);
    end
end

valid_nom = nom_conv & ~isnan(nom_theta);
valid_mc  = ~isnan(mc_th_lo) & ~isnan(mc_th_hi);
valid_rpm = ~isnan(mc_RPMf_lo);
valid_pwr = ~isnan(mc_P_lo) & ~isnan(mc_P_hi);

% Fig 4: pitch + RPM with MC bands ------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 15.5 7.5]);
sub1 = subplot(1,2,1);
fill([V_sweep(valid_mc); flipud(V_sweep(valid_mc))], ...
     [mc_th_hi(valid_mc); flipud(mc_th_lo(valid_mc))], ...
     col_shade, 'EdgeColor', col_edge, 'EdgeAlpha', 0.5, ...
     'LineWidth', 0.6, 'FaceAlpha', 0.5, 'DisplayName', '5th-95th pctl.');
hold on;
plot(V_sweep(valid_nom), nom_theta(valid_nom), '--', 'Color', col_nom, ...
    'LineWidth', 1.8, 'DisplayName', 'Nominal');
plot(V_sweep(valid_mc), mc_th_med(valid_mc), '-', 'Color', col_mc, ...
    'LineWidth', 1.5, 'DisplayName', 'Median');
hold off;
xlabel('Airspeed [m/s]'); ylabel('\theta [deg]');
xlim([0 V_max_plot]); grid on;
legend('Location', 'southwest', 'FontSize', 8, 'Box', 'off');

sub2 = subplot(1,2,2);
fill([V_sweep(valid_rpm); flipud(V_sweep(valid_rpm))], ...
     [mc_RPMf_hi(valid_rpm); flipud(mc_RPMf_lo(valid_rpm))], ...
     col_shade, 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');
hold on;
fill([V_sweep(valid_rpm); flipud(V_sweep(valid_rpm))], ...
     [mc_RPMr_hi(valid_rpm); flipud(mc_RPMr_lo(valid_rpm))], ...
     col_shade_rear, 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');
plot(V_sweep(valid_nom), nom_RPM(valid_nom, 1), '--', 'Color', col_nom, ...
    'LineWidth', 1.5, 'DisplayName', 'Front (nom.)');
plot(V_sweep(valid_nom), nom_RPM(valid_nom, 3), '--', 'Color', col_grey, ...
    'LineWidth', 1.5, 'DisplayName', 'Rear (nom.)');
plot(V_sweep(valid_rpm), mc_RPMf_med(valid_rpm), '-', 'Color', col_mc, ...
    'LineWidth', 1.5, 'DisplayName', 'Front (med.)');
plot(V_sweep(valid_rpm), mc_RPMr_med(valid_rpm), '-', 'Color', col_red, ...
    'LineWidth', 1.5, 'DisplayName', 'Rear (med.)');
hold off;
xlabel('Airspeed [m/s]'); ylabel('Rotor speed [RPM]');
xlim([0 V_max_plot]); grid on;
legend('Location', 'best', 'FontSize', 7, 'Box', 'off');

w = 0.35; gap = 0.12; left1 = 0.10; left2 = left1 + w + gap;
p1 = get(sub1, 'Position'); p2 = get(sub2, 'Position');
set(sub1, 'Position', [left1, p1(2), w, p1(4)]);
set(sub2, 'Position', [left2, p2(2), w, p2(4)]);

drawnow;
exportgraphics(gcf, 'figs/Fig04_trim_sweep.pdf', 'ContentType', 'vector');

% Fig 5: power required with MC band ----------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 12 6.5]);
fill([V_sweep(valid_pwr); flipud(V_sweep(valid_pwr))], ...
     [mc_P_hi(valid_pwr); flipud(mc_P_lo(valid_pwr))], ...
     col_shade, 'EdgeColor', col_edge, 'EdgeAlpha', 0.5, ...
     'LineWidth', 0.6, 'FaceAlpha', 0.5, 'DisplayName', '5th-95th pctl.');
hold on;
plot(V_sweep(valid_pwr), mc_P_med(valid_pwr), '-', 'Color', col_dark_blue, ...
    'LineWidth', 2.0, 'DisplayName', 'Median');
valid_nom_pwr = ~isnan(nom_power) & nom_power > 0;
if any(valid_nom_pwr)
    plot(V_sweep(valid_nom_pwr), nom_power(valid_nom_pwr), '--', 'Color', col_nom, ...
        'LineWidth', 1.5, 'DisplayName', 'Nominal');
end
hold off;
xlabel('Airspeed [m/s]'); ylabel('Total power [W]');
xlim([0 V_max_plot]); grid on;
legend('Location', 'northwest', 'FontSize', 9, 'Box', 'off');

drawnow;
exportgraphics(gcf, 'figs/Fig05_power_bands.pdf', 'ContentType', 'vector');

% Fig 6: P(trim) S-curve with deterministic step overlay --------------------
det_limit = max(V_sweep(nom_conv));
Pt = P_trim(:);

figure('Units', 'centimeters', 'Position', [2 2 12 6.5]);
V_det = [0, det_limit, det_limit + 0.25, max(V_sweep)+2];
P_det = [1, 1, 0, 0];
plot(V_det, P_det, ':', 'Color', col_red, 'LineWidth', 1.5, ...
    'DisplayName', 'Deterministic');
hold on;
plot(V_sweep, Pt, '-', 'Color', col_dark_blue, 'LineWidth', 2.5, ...
    'DisplayName', 'With uncertainty');

thresholds  = [0.95, 0.50];
thresh_labs = {'95%', '50%'};
thresh_sty  = {'--', '-.'};
V_cross     = NaN(1, 2);
for t = 1:numel(thresholds)
    yl = yline(thresholds(t), thresh_sty{t}, thresh_labs{t}, ...
        'Color', [0.45 0.45 0.45], 'LineWidth', 0.7, 'FontSize', 10, ...
        'Alpha', 0.6, ...
        'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom');
    yl.HandleVisibility = 'off';
    idx = find(Pt < thresholds(t), 1, 'first');
    if ~isempty(idx) && idx > 1
        V_cross(t) = interp1(Pt(idx-1:idx), V_sweep(idx-1:idx), thresholds(t));
        xl = xline(V_cross(t), ':', 'Color', [0.80 0.80 0.80], 'LineWidth', 0.4);
        xl.HandleVisibility = 'off';
    end
end
hold off;
xlabel('Airspeed [m/s]'); ylabel('P(trim)');
xlim([0 max(V_sweep)]); ylim([-0.05 1.05]); grid on;
legend('Location', 'southwest', 'FontSize', 9, 'Box', 'off');

drawnow;
exportgraphics(gcf, 'figs/Fig06_Ptrim_Scurve.pdf', 'ContentType', 'vector');

% Console table -------------------------------------------------------------
fprintf('\nMaximum airspeed at selected confidence levels\n');
fprintf('%-15s %-15s %-15s\n', 'Confidence', 'P(trim)', 'V_max [m/s]');
fprintf('%s\n', repmat('-', 1, 47));
for t = 1:numel(thresholds)
    if ~isnan(V_cross(t))
        fprintf('%-15s %-15.3f %-15.2f\n', thresh_labs{t}, thresholds(t), V_cross(t));
    else
        fprintf('%-15s %-15.3f %-15s\n', thresh_labs{t}, thresholds(t), 'N/A');
    end
end
fprintf('%-15s %-15.3f %-15.2f\n', 'Deterministic', 1.000, det_limit);

fprintf('\nSaved: figs/Fig04_trim_sweep.pdf, Fig05_power_bands.pdf, Fig06_Ptrim_Scurve.pdf\n');