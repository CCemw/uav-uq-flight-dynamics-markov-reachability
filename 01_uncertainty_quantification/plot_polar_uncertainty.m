% PLOT_POLAR_UNCERTAINTY
% Produces:
%   figs/Fig04_airfoil_geometry.pdf    The three measured blade-section
%                                      airfoil profiles plotted on top of
%                                      each other.
%   figs/Fig05_Cl_Cd_polars.pdf        The lift and drag coefficients of
%                                      all three airfoils against angle of
%                                      attack at M = 0.30, shown side by
%                                      side.
%   figs/Fig05b_Cl_Cd_sigma.pdf        The same polars reduced to a mean
%                                      curve with a pointwise standard
%                                      deviation band derived from the
%                                      range of the three airfoils divided
%                                      by the d2 constant.
%
% Reads:  ../00_inputs/polars/afl{13,14,15}.dat      (geometry coordinates)
%         ../00_inputs/polars/afl{13,14,15}.c81.txt  (CL and CD tables)

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

col_nom      = [0.00 0.00 0.00];
col_v1       = [0.00 0.45 0.70];
col_v2       = [0.90 0.60 0.00];
col_shade    = [0.50 0.72 0.88];
col_edge     = [0.25 0.50 0.70];

labels   = {'Variant 1', 'Nominal', 'Variant 2'};
colours  = [col_v1; col_nom; col_v2];
styles   = {':', '-', ':'};
lwidths  = [1.6, 2.0, 1.6];



% Configuration============================================================
geom_files = {'../00_inputs/polars/afl13.dat', ...
              '../00_inputs/polars/afl14.dat', ...
              '../00_inputs/polars/afl15.dat'};
c81_files  = {'../00_inputs/polars/afl13.c81.txt', ...
              '../00_inputs/polars/afl14.c81.txt', ...
              '../00_inputs/polars/afl15.c81.txt'};
M_plot = 0.30;     % Mach number at which to show the polars
d2_n3  = 1.693;    % Tippett constant for n=3 samples
% =========================================================================



% Fig 4: airfoil geometry ---------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 15.5 5.5]);
hold on;
for k = 1:3
    raw = readmatrix(geom_files{k}, 'FileType', 'text');
    xy  = raw(2:end, :);
    plot(xy(:,1), xy(:,2), styles{k}, ...
         'Color', colours(k,:), 'LineWidth', lwidths(k), ...
         'DisplayName', labels{k});
end
hold off;
axis equal;
xlim([-0.02 1.02]);
xlabel('x/c'); ylabel('y/c');
ax = gca;
pos = get(ax, 'Position');
set(ax, 'Position', [pos(1)-0.03, pos(2), pos(3)+0.03, pos(4)]);
legend('Location', 'northeast', 'FontSize', 10, 'Box', 'off');
drawnow;
exportgraphics(gcf, 'figs/Fig04_airfoil_geometry.pdf', 'ContentType', 'vector');

for k = 1:3
    A = readmatrix(c81_files{k});
    AF(k).M_CL    = A(1,   2:11);
    AF(k).alphaCL = A(2:100, 1);
    AF(k).CL      = A(2:100, 2:11);
    AF(k).M_CD    = A(101, 2:11);
    AF(k).alphaCD = A(102:200, 1);
    AF(k).CD      = A(102:200, 2:11);
end
alphaCL = AF(1).alphaCL;
alphaCD = AF(1).alphaCD;
[~, iM_CL] = min(abs(AF(1).M_CL - M_plot));
[~, iM_CD] = min(abs(AF(1).M_CD - M_plot));

% Fig 5: Cl and Cd polars for the three airfoils ---------------------------
figure('Units', 'centimeters', 'Position', [2 2 15.5 8]);

sub1 = subplot(1,2,1);
hold on;
for k = 1:3
    plot(AF(k).alphaCL, AF(k).CL(:, iM_CL), styles{k}, ...
         'Color', colours(k,:), 'LineWidth', lwidths(k), ...
         'DisplayName', labels{k});
end
hold off;
xlabel('\alpha [deg]'); ylabel('C_l');
xlim([-25 25]); grid on;
set(gca, 'PlotBoxAspectRatio', [1 1 1]);

sub2 = subplot(1,2,2);
hold on;
for k = 1:3
    plot(AF(k).alphaCD, AF(k).CD(:, iM_CD), styles{k}, ...
         'Color', colours(k,:), 'LineWidth', lwidths(k), ...
         'HandleVisibility', 'off');
end
hold off;
xlabel('\alpha [deg]'); ylabel('C_d');
xlim([-25 25]); grid on;
set(gca, 'PlotBoxAspectRatio', [1 1 1]);

w = 0.33; gap = 0.12; left1 = 0.10; left2 = left1 + w + gap;
p1 = get(sub1, 'Position'); p2 = get(sub2, 'Position');
set(sub1, 'Position', [left1, p1(2), w, p1(4)]);
set(sub2, 'Position', [left2, p2(2), w, p2(4)]);
legend(sub1, 'Location', 'northwest', 'FontSize', 10, 'Box', 'off');

drawnow;
exportgraphics(gcf, 'figs/Fig05_Cl_Cd_polars.pdf', 'ContentType', 'vector');

% Fig 5b: mean polars with pointwise sigma bands ----------------------------
CLs = cat(2, AF(1).CL(:,iM_CL), AF(2).CL(:,iM_CL), AF(3).CL(:,iM_CL));
CDs = cat(2, AF(1).CD(:,iM_CD), AF(2).CD(:,iM_CD), AF(3).CD(:,iM_CD));

CL_mean  = mean(CLs, 2);
CD_mean  = mean(CDs, 2);
sigma_CL = (max(CLs,[],2) - min(CLs,[],2)) / d2_n3;
sigma_CD = (max(CDs,[],2) - min(CDs,[],2)) / d2_n3;

figure('Units', 'centimeters', 'Position', [2 2 15.5 8]);

sub1 = subplot(1,2,1);
h_fill = fill([alphaCL; flipud(alphaCL)], ...
              [CL_mean + sigma_CL; flipud(CL_mean - sigma_CL)], ...
              col_shade, 'EdgeColor', col_edge, 'EdgeAlpha', 0.5, ...
              'LineWidth', 0.6, 'FaceAlpha', 0.6, ...
              'DisplayName', 'Mean \pm \sigma');
hold on;
h_mean = plot(alphaCL, CL_mean, '-', 'Color', col_nom, ...
              'LineWidth', 2.0, 'DisplayName', 'Mean');
hold off;
xlabel('\alpha [deg]'); ylabel('C_l');
xlim([-25 25]); grid on;
set(gca, 'PlotBoxAspectRatio', [1 1 1]);

sub2 = subplot(1,2,2);
fill([alphaCD; flipud(alphaCD)], ...
     [CD_mean + sigma_CD; flipud(CD_mean - sigma_CD)], ...
     col_shade, 'EdgeColor', col_edge, 'EdgeAlpha', 0.5, ...
     'LineWidth', 0.6, 'FaceAlpha', 0.6, ...
     'HandleVisibility', 'off');
hold on;
plot(alphaCD, CD_mean, '-', 'Color', col_nom, 'LineWidth', 2.0, ...
     'HandleVisibility', 'off');
hold off;
xlabel('\alpha [deg]'); ylabel('C_d');
xlim([-25 25]); grid on;
set(gca, 'PlotBoxAspectRatio', [1 1 1]);

p1 = get(sub1, 'Position'); p2 = get(sub2, 'Position');
set(sub1, 'Position', [left1, p1(2), w, p1(4)]);
set(sub2, 'Position', [left2, p2(2), w, p2(4)]);
legend(sub1, [h_fill, h_mean], {'Mean \pm \sigma', 'Nominal'}, ...
       'Location', 'northwest', 'FontSize', 10, 'Box', 'off');

drawnow;
exportgraphics(gcf, 'figs/Fig05b_Cl_Cd_sigma.pdf', 'ContentType', 'vector');