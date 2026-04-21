% MISSION_DTMC
% Runs a step-by-step probability model of the UAV flying through a
% mission, phase by phase. At every gust encounter the aircraft either stays
% in normal flight, drops into a reduced-capability state, aborts, or
% crashes, with the probabilities of each outcome taken from the trim
% results and the local gust statistics. The probabilities are carried
% forward through the whole mission so the user can see how likely each
% final outcome is.
%
% Produces:
%   figs/Fig11a_det_post.pdf     How the probabilities of each flight state
%                                evolve through the mission phases when no
%                                uncertainty is assumed. Shown for two cruise
%                                speeds (12 and 17 m/s) side by side.
%
%   figs/Fig11b_unc_pre.pdf      The same evolution when aerodynamic
%                                uncertainty is included, shown BEFORE the
%                                between-phase contingency routing is
%                                applied. This reveals how much reduced-
%                                capability probability builds up during
%                                each phase.
%
%   figs/Fig11c_unc_post.pdf     The same uncertain case, shown AFTER the
%                                contingency routing is applied. This is
%                                what the aircraft actually experiences,
%                                where any residual reduced-capability
%                                probability at the end of a phase is sent
%                                through an abort sub-mission.
%
%   figs/Fig12_det_vs_unc.pdf    Crash rate per flight hour plotted against
%                                cruise speed, with and without aerodynamic
%                                uncertainty. Overlaid with the JARUS SORA
%                                SAIL thresholds for comparison.
%
%   console: terminal outcomes   Final probabilities of ending the mission
%                                in each state at the two reference cruise
%                                speeds, plus the equivalent crash rate per
%                                flight hour.
%
%   console: V_max compliance    The maximum cruise speed at which the
%                                aircraft meets each SAIL level, compared
%                                between the deterministic case and the
%                                case with aerodynamic uncertainty.
%
% Reads:  ../02_trim_monte_carlo/results/mc_trim_results_QAV250.mat

clear; close all; clc;
if ~isfolder('figs'); mkdir('figs'); end   % make the figs folder if it doesn't exist yet


% USER SETTINGS
% ===============================================================================

% Mission profile: {phase name, V_nom [m/s], altitude [m], duration [s]}
mission = {
    'Takeoff', 3,  10, 30;
    'Cruise',  12, 60, 600;
    'Descent', 12, 40, 120;
    'Loiter',  5,  30, 180;
    'Landing', 3,  10, 30
};

% Abort sub-phase: entered from DF at the end of any non-final phase
abort_phase = struct('V_nom', 5, 'h', 30, 'duration', 60);

% Contingency routing parameters
eta     = 0.80;  % successful-landing probability from Abort state
eta_DF  = 0.40;  % successful-landing probability from Degraded Flight
s       = 0.20;  % severity-routing parameter (fraction of trim-loss events
                 % that are routed straight to the terminal branch rather than to DF)

% Threshold separating classification of Normal Flight (NF) from "Degraded Flight"
P_deg = 0.95;    % P(trim) must be >= P_deg for the encounter to count as NF

% Atmospheric parameters (low-altitude Dryden with urban correction)
W20_kt         = 15;                         % wind at 20 ft [knots]
h_urban_table  = [10, 30, 50, 80, 100];      % altitude breakpoints [m]
f_urban_table  = [2.5, 2.0, 1.65, 1.5, 1.3]; % urban intensity factor f_u(h)
f_Lu           = 0.5;                        % urban correlation-length factor

% Numerical integration of the gust distribution
N_gust_pts = 200;
gust_range = 5;   % integration window in units of sigma_u

% Desired cruise speeds to analyse
V_pair   = [12, 17];    % paired cases for state-evolution figures
V_sweep  = 8:1:22;      % sweep for the lambda_crash vs V_cruise plot

% SAIL thresholds (JARUS SORA v2.5, /flight-hour)
sail_targets = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
sail_names   = {'I','II','III','IV','V','VI'};

% USER INPUT END
% ====================================================================================

mc = load('../02_trim_monte_carlo/results/mc_trim_results_QAV250.mat');
fprintf('Loaded trim MC: %d speeds, %d samples.\n', numel(mc.V_sweep), mc.N_MC);


atm.V_grid       = mc.V_sweep(:)';
atm.P_trim       = max(0, min(1, mc.P_trim(:)'));
atm.W20_ms       = W20_kt * 0.5144;
atm.h_urban_table = h_urban_table;
atm.f_urban_table = f_urban_table;
atm.f_Lu         = f_Lu;
atm.P_deg        = P_deg;
atm.N_gust_pts   = N_gust_pts;
atm.gust_range   = gust_range;

% Deterministic reference: P(trim) = 1 up to the deterministic trim limit,
% then a hard step to zero. Trim limit taken from the nominal trim sweep.
det_limit = max(mc.V_sweep(mc.nom_conv));
atm_det            = atm;
atm_det.V_grid     = 0:0.25:35;
atm_det.P_trim     = ones(size(atm_det.V_grid));
atm_det.P_trim(atm_det.V_grid >= det_limit + 0.25) = 0;

% Mission time in flight-hours
T_mission_hrs = sum(cell2mat(mission(:,4))) / 3600;


% STATE EVOLUTION
n_phases   = size(mission, 1);
n_cases    = 2;                           % {deterministic, with uncertainty}
atm_cases  = {atm_det, atm};

history_post = cell(n_cases, numel(V_pair));
history_pre  = cell(1,       numel(V_pair));

for ci = 1:n_cases
    for vi = 1:numel(V_pair)
        [pi_post, pi_pre] = evolve_mission(V_pair(vi), eta, eta_DF, s, ...
                                           mission, abort_phase, atm_cases{ci});
        history_post{ci, vi} = pi_post;
        if ci == 2  % uncertain case: also store the pre-routing history
            history_pre{vi} = pi_pre;
        end
    end
end

plot_state_evolution(history_post(1,:), V_pair, mission, ...
                     'figs/Fig11a_det_post.pdf');
plot_state_evolution(history_pre,       V_pair, mission, ...
                     'figs/Fig11b_unc_pre.pdf');
plot_state_evolution(history_post(2,:), V_pair, mission, ...
                     'figs/Fig11c_unc_post.pdf');



% LAMBDA_CRASH VS CRUISE SPEED (FIG 12)
lambda_det = zeros(numel(V_sweep), 1);
lambda_unc = zeros(numel(V_sweep), 1);
for vi = 1:numel(V_sweep)
    pi_d = run_mission(V_sweep(vi), eta, eta_DF, s, mission, abort_phase, atm_det);
    pi_u = run_mission(V_sweep(vi), eta, eta_DF, s, mission, abort_phase, atm);
    lambda_det(vi) = max(pi_d(4) / T_mission_hrs, 1e-15);
    lambda_unc(vi) = pi_u(4) / T_mission_hrs;
end

plot_lambda_vs_V(V_sweep, lambda_det, lambda_unc, ...
                 sail_targets, sail_names, V_pair, ...
                 'figs/Fig12_det_vs_unc.pdf');

fprintf('\nTerminal outcomes (with eta=%.2f, eta_DF=%.2f, s=%.2f)\n', ...
        eta, eta_DF, s);
fprintf('%-8s  %12s  %12s  %12s  %14s\n', ...
        'V [m/s]', 'P(NF)', 'P(Abort)', 'P(Crash)', 'lambda/FH');
fprintf('%s\n', repmat('-', 1, 64));
for V = V_pair
    pi_f = run_mission(V, eta, eta_DF, s, mission, abort_phase, atm);
    fprintf('%-8d  %12.4e  %12.4e  %12.4e  %14.4e\n', ...
            V, pi_f(1), pi_f(3), pi_f(4), pi_f(4)/T_mission_hrs);
end

fprintf('\nMaximum cruise speed compliant with each SAIL level\n');
fprintf('%-10s  %-14s  %-14s\n', 'SAIL', 'Deterministic', 'With unc.');
fprintf('%s\n', repmat('-', 1, 44));
for si = 1:numel(sail_targets)
    v_det = V_sweep(lambda_det < sail_targets(si));
    v_unc = V_sweep(lambda_unc < sail_targets(si));
    det_str = ternary(~isempty(v_det), sprintf('%d m/s', max(v_det)), '---');
    unc_str = ternary(~isempty(v_unc), sprintf('%d m/s', max(v_unc)), '---');
    fprintf('%-10s  %-14s  %-14s\n', ...
            ['SAIL ' sail_names{si}], det_str, unc_str);
end



% LOCAL FUNCTIONS
% =================================================================
function [history_post, history_pre] = evolve_mission(V_cruise, eta, eta_DF, s, ...
                                                      mission, abort_phase, atm)
% Evolve the DTMC through the mission, logging pre- and post-routing state
% probabilities at the end of each phase.
    n_phases       = size(mission, 1);
    history_post   = zeros(n_phases + 1, 4);
    history_pre    = zeros(n_phases + 1, 4);
    history_post(1,:) = [1, 0, 0, 0];
    history_pre(1,:)  = [1, 0, 0, 0];
    pi_state          = [1, 0, 0, 0];

    for ph = 1:n_phases
        V_nom    = mission{ph, 2};
        if ph == 2                    % "Cruise" phase uses the commanded V_cruise
            V_nom = V_cruise;
        end
        h        = mission{ph, 3};
        duration = mission{ph, 4};

        [pN, pDF_pre, pL, N_enc] = compute_phase_probabilities(V_nom, h, duration, atm);
        T_k = build_transition_matrix(pN, pDF_pre, pL, eta, eta_DF, s);

        for n = 1:N_enc
            pi_state = pi_state * T_k;
        end
        history_pre(ph+1, :) = pi_state;

        % Between-phase routing: any residual DF population enters the abort
        % sub-phase, which terminates in either a successful abort landing or
        % a crash. Landing is assumed to fully absorb DF 
        if ph < n_phases && pi_state(2) > 1e-15
            pi_DF = pi_state(2);
            [pNa, pDFa, pLa, Na] = compute_phase_probabilities(...
                abort_phase.V_nom, abort_phase.h, abort_phase.duration, atm);
            T_a = build_transition_matrix(pNa, pDFa, pLa, eta, eta_DF, s);
            pa  = [0, 1, 0, 0];
            for n = 1:Na
                pa = pa * T_a;
            end
            p_ok = pa(1) + pa(2) + pa(3);
            pi_state = [pi_state(1), ...
                        0, ...
                        pi_state(3) + pi_DF * p_ok, ...
                        pi_state(4) + pi_DF * pa(4)];
        end
        history_post(ph+1, :) = pi_state;
    end
end


function pi_final = run_mission(V_cruise, eta, eta_DF, s, mission, abort_phase, atm)
% Evolve the DTMC and return only the final state distribution.
    [hist_post, ~] = evolve_mission(V_cruise, eta, eta_DF, s, ...
                                     mission, abort_phase, atm);
    pi_final = hist_post(end, :);
end


function T = build_transition_matrix(pN, pDF_pre, pL, eta, eta_DF, s)
% 4x4 row-stochastic transition matrix, rows/cols = {NF, DF, Abort, Crash}.
    T = [pN,  pDF_pre + (1-s)*pL,  s*eta*pL,       s*(1-eta)*pL;
         pN,  pDF_pre,             eta_DF*pL,      (1-eta_DF)*pL;
         0,   0,                   1,              0;
         0,   0,                   0,              1];
end


function [p_N, p_DF, p_L, N_enc] = compute_phase_probabilities(V_nom, h, duration, atm)
% Probabilities that a single gust encounter at this phase produces an
% NF, DF, or trim-loss outcome, obtained by convolving the P(trim) envelope
% with a Gaussian gust distribution at local Dryden parameters. N_enc is the
% number of independent encounters assumed to occur during the phase.
    [sigma_u, L_u, N_enc] = dryden_urban_params(V_nom, h, duration, atm);

    dv  = linspace(-atm.gust_range*sigma_u, atm.gust_range*sigma_u, atm.N_gust_pts);
    ddv = dv(2) - dv(1);
    p_N = 0; p_DF = 0; p_L = 0;

    for i = 1:atm.N_gust_pts
        V_tot = V_nom + dv(i);
        Pt    = lookup_Ptrim(V_tot, atm);
        w     = normpdf(dv(i), 0, sigma_u) * ddv;
        if Pt >= atm.P_deg
            p_N  = p_N  + Pt * w;
        else
            p_DF = p_DF + Pt * w;
        end
        p_L  = p_L + (1 - Pt) * w;
    end
end


function [sigma_u, L_u, N_enc] = dryden_urban_params(V_nom, h, duration, atm)
% Low-altitude Dryden standard deviation and correlation length with urban
% intensity scaling, and resulting encounter count over the phase duration.
    h_ft   = h * 3.281;
    sigma_d = (0.1 * atm.W20_ms) / (0.177 + 0.000823*h_ft)^0.4;
    L_d_m   = (h_ft / (0.177 + 0.000823*h_ft)^1.2) * 0.3048;

    f_u = interp1(atm.h_urban_table, atm.f_urban_table, h, 'linear', 'extrap');
    f_u = max(1.0, f_u);

    sigma_u = sigma_d * f_u;
    L_u     = L_d_m * atm.f_Lu;

    if V_nom > 0.5
        N_enc = max(1, ceil(duration / (L_u / V_nom)));
    else
        N_enc = max(1, ceil(duration / 5));   % near-hover case
    end
end


function Pt = lookup_Ptrim(V, atm)
% Linearly interpolate the P(trim) S-curve at V (outside grid is clamped..)
    if V < 0
        Pt = 1;
    elseif V > atm.V_grid(end)
        Pt = 0;
    else
        Pt = interp1(atm.V_grid, atm.P_trim, V, 'linear', 1);
        Pt = max(0, min(1, Pt));
    end
end


function plot_state_evolution(hist_list, V_pair, mission, out_file)
% Two-panel state-evolution plot (V_pair(1) and V_pair(2)).
    col_nf    = [0.00 0.45 0.70];
    col_df    = [0.90 0.60 0.00];
    col_abort = [0.00 0.62 0.45];
    col_crash = [0.80 0.20 0.20];
    state_colours = [col_nf; col_df; col_abort; col_crash];
    state_names   = {'Normal Flight', 'Degraded Flight', 'Abort', 'Crash'};
    state_styles  = {'-o', '-s', '-^', '-v'};

    n_phases = size(mission, 1);
    x_pts    = 0:n_phases;
    labels   = [{'Start'}, mission(:,1)'];
    y_lim    = [-0.02 1.05];

    figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 15.5 7.5]);
    for vi = 1:numel(V_pair)
        subplot(1, 2, vi);
        h = hist_list{vi};
        for st = 1:4
            plot(x_pts, h(:,st), state_styles{st}, ...
                 'Color', state_colours(st,:), 'LineWidth', 1.6, ...
                 'MarkerSize', 4, 'MarkerFaceColor', state_colours(st,:), ...
                 'DisplayName', state_names{st});
            hold on;
        end
        for pd = 0:n_phases
            line([pd+0.5 pd+0.5], y_lim, 'Color', [0.3 0.3 0.3], ...
                 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        hold off;
        set(gca, 'XTick', x_pts, 'XTickLabel', labels, 'XTickLabelRotation', 30);
        ylabel('State probability');
        xlim([-0.5 n_phases+0.5]); ylim(y_lim);
        set(gca, 'YGrid', 'on', 'XGrid', 'off');
        if vi == 1
            legend('Location', 'east', 'FontSize', 8, 'Box', 'off');
        end
    end

    sub1 = subplot(1,2,1); sub2 = subplot(1,2,2);
    w = 0.38; gap = 0.08; left1 = 0.08; left2 = left1 + w + gap;
    p1 = get(sub1, 'Position'); p2 = get(sub2, 'Position');
    set(sub1, 'Position', [left1, p1(2)+0.05, w, p1(4)-0.05]);
    set(sub2, 'Position', [left2, p2(2)+0.05, w, p2(4)-0.05]);

    drawnow;
    exportgraphics(gcf, out_file, 'ContentType', 'vector');
end


function plot_lambda_vs_V(V_sweep, lambda_det, lambda_unc, ...
                          sail_targets, sail_names, V_pair, out_file)
% lambda_crash vs cruise speed with SAIL thresholds.
    col_mc  = [0.00 0.45 0.70];
    col_red = [0.80 0.20 0.20];
    col_sail = [0.35 0.35 0.35];

    figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 13 8]);
    semilogy(V_sweep, lambda_unc, '-o', 'Color', col_red, 'LineWidth', 2.0, ...
             'MarkerSize', 4, 'MarkerFaceColor', col_red, ...
             'DisplayName', 'With manufacturing uncertainty');
    hold on;
    semilogy(V_sweep, lambda_det, '--s', 'Color', col_mc, 'LineWidth', 1.8, ...
             'MarkerSize', 4, 'MarkerFaceColor', col_mc, ...
             'DisplayName', 'Deterministic');
    for si = 1:numel(sail_targets)
        yline(sail_targets(si), '--', ['SAIL ' sail_names{si}], ...
              'Color', col_sail, 'LineWidth', 0.9, 'FontSize', 8, ...
              'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
    end
    for v = V_pair
        xline(v, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5, ...
              'HandleVisibility', 'off');
    end
    hold off;
    xlabel('Cruise speed [m/s]');
    ylabel('\lambda_{crash} [per flight hour]');
    ylim([1e-8 1e1]);
    xlim([V_sweep(1) V_sweep(end)]);
    set(gca, 'YTick', 10.^(-8:1), 'YMinorTick', 'off', ...
             'YGrid', 'off', 'YMinorGrid', 'off', 'XGrid', 'off');
    legend('Location', 'northwest', 'FontSize', 9, 'Box', 'off');

    drawnow;
    exportgraphics(gcf, out_file, 'ContentType', 'vector');
end


function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end