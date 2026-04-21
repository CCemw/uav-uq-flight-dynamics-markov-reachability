% CLOSED_LOOP_MANOEUVRE_MC
% Closed-loop pitch response to a symmetric rear-rotor RPM pulse, with
% optional (on/off) Monte Carlo aerodynamic uncertainty propagated through the
% perturbed trim and the time-domain response.
%
% Trims internally (ie does not load trim file). The cascaded angle/rate PID
% gains are exposed at the top alongside V_trim, the pulse magnitude and
% pulse duration so the user can change them accordingly.
%
% USER SETTINGS at the top: change V_trim, dRPM_pulse, t_pulse_dur, RUN_MODE,
% and the controller gains. Re-run for each case. Filename is built from
% V_trim and RUN_MODE.
%
% Reads:  ../00_inputs/platform/platform_parameters.m
%         ../00_inputs/polars/afl{13,14,15}.c81.txt
% Writes: results/manoeuvre_V<V_trim>_<RUN_MODE>_QAV250.mat  
% >>> used for plotting

clear; clc; close all;

% USER SETTINGS -------------------------------------------------------------
V_trim       = 15;          % trim airspeed [m/s]
dRPM_pulse   = 300;         % rear-rotor pulse magnitude [RPM]
t_pulse_dur  = 0.3;         % pulse duration [s]
RUN_MODE     = 'mc';        % 'nominal' or 'mc'

% Cascaded angle/rate PID gains for pitch control
Kp_angle     = 7.0;         % outer loop: pitch angle [-]
Kp_rate      = 24.0;        % inner loop: pitch rate proportional [-]
Ki_rate      = 3.0;         % inner loop: pitch rate integral [-]
Kd_rate      = 0.8;         % inner loop: pitch rate derivative [-]
f_D_cutoff   = 50;          % derivative low-pass cutoff [Hz]
% --------------------------------------------------------------------------

% Configuration -------------------------------------------------------------
N_MC          = 500;
RNG_SEED      = 42;
t_pulse       = 0.5;        % pulse start time [s]
t_sim         = 2.5;        % simulation length [s]
dt            = 0.001;      % integration step [s]
sigma_twist   = deg2rad(1.0);
d2_n3         = 1.693;
CD_floor      = 0.005;

% Actuator limits (only used here for the CL manoeuvre)
RPM_min       = 2000;
RPM_max       = 20000;

% Divergence abort (user-set): if |theta| exceeds this the run is logged
% as "diverged" and the time history is held constant from that point on
theta_diverge = 60;          % [deg]

% Trim-solver iteration bounds (inf)
RPM_bounds = [-50000, 50000];

airfoil_files = {'afl13.c81.txt', 'afl14.c81.txt', 'afl15.c81.txt'};
nominal_idx   = 2;

out_file = sprintf('results/manoeuvre_V%d_%s_QAV250.mat', V_trim, RUN_MODE);
airfoil_files = {'afl13.c81.txt', 'afl14.c81.txt', 'afl15.c81.txt'};
nominal_idx   = 2;

out_file = sprintf('results/manoeuvre_V%d_%s_QAV250.mat', V_trim, RUN_MODE);

% Setup ---------------------------------------------------------------------
p = platform_parameters();
p.mu_switch = 0.01;

for k = 1:numel(airfoil_files)
    A = readmatrix(airfoil_files{k});
    AF(k).M       = A(1, 2:11);
    AF(k).alphaCL = A(2:100, 1);
    AF(k).CL      = A(2:100, 2:11);
    AF(k).alphaCD = A(102:200, 1);
    AF(k).CD      = A(102:200, 2:11);
end
nom_CL = AF(nominal_idx).CL;
nom_CD = AF(nominal_idx).CD;

CLs = cat(3, AF(1).CL, AF(2).CL, AF(3).CL);
CDs = cat(3, AF(1).CD, AF(2).CD, AF(3).CD);
sigma_CL_grid = (max(CLs,[],3) - min(CLs,[],3)) / d2_n3;
sigma_CD_grid = (max(CDs,[],3) - min(CDs,[],3)) / d2_n3;

p.CL_alpha_vec = AF(nominal_idx).alphaCL;
p.CL_n_alpha   = length(p.CL_alpha_vec);
p.CL_M_vec     = AF(nominal_idx).M;
p.CL_n_M       = size(nom_CL, 2);
p.CL_alpha_min = p.CL_alpha_vec(1);
p.CL_alpha_max = p.CL_alpha_vec(end);
p.CL_M_min     = p.CL_M_vec(1);
p.CL_M_max     = p.CL_M_vec(end);
p.CL_d_alpha   = diff(p.CL_alpha_vec);
p.CL_d_M       = diff(p.CL_M_vec);

p.CD_alpha_vec = AF(nominal_idx).alphaCD;
p.CD_n_alpha   = length(p.CD_alpha_vec);
p.CD_M_vec     = AF(nominal_idx).M;
p.CD_n_M       = size(nom_CD, 2);
p.CD_alpha_min = p.CD_alpha_vec(1);
p.CD_alpha_max = p.CD_alpha_vec(end);
p.CD_M_min     = p.CD_M_vec(1);
p.CD_M_max     = p.CD_M_vec(end);
p.CD_d_alpha   = diff(p.CD_alpha_vec);
p.CD_d_M       = diff(p.CD_M_vec);

for r = 1:4
    p.rotor_CL{r} = nom_CL;
    p.rotor_CD{r} = nom_CD;
end
p.RPM_hover = find_hover_rpm(p);

fprintf('Hover RPM: %d | Mass: %.3f kg\n', p.RPM_hover, p.m);

% Nominal trim sweep up to V_trim+2 m/s ------------------------------------
V_sweep    = 0:1:max(V_trim+2, 10);
nV         = numel(V_sweep);
nom_theta  = NaN(nV,1);
nom_phi    = NaN(nV,1);
nom_RPM    = NaN(nV,4);
nom_conv   = false(nV,1);

x_prev = [0; 0; p.RPM_hover; p.RPM_hover; p.RPM_hover; p.RPM_hover];

fprintf('(Nominal) Trim Speed Sweep: %.2f to %.2f m/s, %d points\n', ...
        V_sweep(1), V_sweep(end), nV);
for vi = 1:nV
    V = V_sweep(vi);
    guesses = build_guesses(x_prev, vi, nom_conv, nom_theta, nom_phi, ...
                            nom_RPM, p.RPM_hover);
    [x_best, conv] = try_guesses(guesses, V, p, RPM_bounds);
    nom_conv(vi) = conv;
    if conv
        nom_theta(vi) = rad2deg(x_best(1));
        nom_phi(vi)   = rad2deg(x_best(2));
        nom_RPM(vi,:) = x_best(3:6)';
        x_prev = x_best;
    end
end
idx_trim = find(V_sweep == V_trim & nom_conv', 1);
if isempty(idx_trim)
    error('Trim failed at V = %d m/s.', V_trim);
end
theta_trim = deg2rad(nom_theta(idx_trim));
phi_trim   = deg2rad(nom_phi(idx_trim));
RPM_trim   = nom_RPM(idx_trim,:)';
fprintf('Trim at V=%d m/s: theta=%+.2f deg, RPM=[%.0f %.0f %.0f %.0f]\n', ...
        V_trim, rad2deg(theta_trim), RPM_trim');

% CRN draws ----------------------------------------------
use_mc = strcmpi(RUN_MODE, 'mc');
N_runs = 1;
if use_mc
    N_runs = N_MC;
    rng(RNG_SEED, 'twister');
    z_CL_all    = randn(4, N_MC);
    z_CD_all    = randn(4, N_MC);
    z_twist_all = randn(4, N_MC);
end

% Time loop allocations ----------------------------------------------------
t_vec    = 0:dt:t_sim;
N_steps  = numel(t_vec);
all_theta = NaN(N_runs, N_steps);
all_q     = NaN(N_runs, N_steps);
all_V     = NaN(N_runs, N_steps);
all_phi   = NaN(N_runs, N_steps);
run_conv       = false(N_runs, 1);
diverged       = false(N_runs, 1);
peak_overshoot = NaN(N_runs, 1);
settling_time  = NaN(N_runs, 1);
ss_error       = NaN(N_runs, 1);
mc_trim_theta  = NaN(N_runs, 1);
mc_trim_RPM    = NaN(N_runs, 4);
mc_peak_q      = NaN(N_runs, 1);

% Simulate -----------------------------------------------------------------
fprintf('(%s) Closed-loop manoeuvre: V=%d m/s, +%d RPM x %.1fs pulse, %d run(s)\n', ...
        upper(RUN_MODE), V_trim, dRPM_pulse, t_pulse_dur, N_runs);
alpha_D = dt / (dt + 1/(2*pi*f_D_cutoff));

tic_total = tic;
for s = 1:N_runs
    p_run = p;
    if use_mc
        for r = 1:4
            p_run.rotor_CL{r}      = nom_CL + z_CL_all(r,s) * sigma_CL_grid;
            p_run.rotor_CD{r}      = max(CD_floor, ...
                                          nom_CD + z_CD_all(r,s) * sigma_CD_grid);
            twist_offset           = z_twist_all(r,s) * sigma_twist;
            p_run.theta_blade_r{r} = p.theta_blade + twist_offset;
            p_run.theta_mat_r{r}   = repmat(p_run.theta_blade_r{r}, 1, p.M_azi);
        end
        [x_mc, conv_mc] = try_guesses({[theta_trim; phi_trim; RPM_trim]}, ...
                                      V_trim, p_run, RPM_bounds);
        if ~conv_mc, continue; end
        theta_s = x_mc(1); phi_s = x_mc(2); RPM_s = x_mc(3:6);
        mc_trim_theta(s) = rad2deg(theta_s);
        mc_trim_RPM(s,:) = RPM_s';
    else
        theta_s = theta_trim;
        phi_s   = phi_trim;
        RPM_s   = RPM_trim;
    end
    state = [V_trim*cos(theta_s); 0; V_trim*sin(theta_s);
             0; 0; 0;
             phi_s; theta_s];

    int_rate_err = 0;
    prev_q_filt  = 0;
    did_diverge  = false;

    for n = 1:N_steps
        t                 = t_vec(n);
        all_theta(s, n)   = rad2deg(state(8));
        all_q(s, n)       = rad2deg(state(5));
        all_V(s, n)       = norm(state(1:3));
        all_phi(s, n)     = rad2deg(state(7));

        if abs(rad2deg(state(8))) > theta_diverge || any(isnan(state))
            did_diverge = true;
            if n < N_steps
                all_theta(s, n+1:end) = all_theta(s, n);
                all_q(s, n+1:end)     = all_q(s, n);
                all_V(s, n+1:end)     = all_V(s, n);
                all_phi(s, n+1:end)   = all_phi(s, n);
            end
            break;
        end
        if n == N_steps, break; end

        RPM_cmd = RPM_s;
        if t >= t_pulse && t <= t_pulse + t_pulse_dur
            RPM_cmd(3) = RPM_cmd(3) + dRPM_pulse;
            RPM_cmd(4) = RPM_cmd(4) + dRPM_pulse;
        end

        % Cascaded angle -> rate PID with derivative low-pass
        theta_err    = theta_s - state(8);
        q_cmd        = Kp_angle * theta_err;
        q_meas       = state(5);
        rate_err     = q_cmd - q_meas;
        P_out        = Kp_rate * rate_err;
        int_rate_err = max(-50, min(50, int_rate_err + rate_err * dt));
        I_out        = Ki_rate * int_rate_err;
        q_filt       = alpha_D * q_meas + (1 - alpha_D) * prev_q_filt;
        dq_filt      = (q_filt - prev_q_filt) / dt;
        prev_q_filt  = q_filt;
        D_out        = -Kd_rate * dq_filt;
        u_ctrl       = P_out + I_out + D_out;
        dRPM_ctrl    = u_ctrl * 60 / (2*pi);

        RPM_cmd(1)  = RPM_cmd(1) + dRPM_ctrl;
        RPM_cmd(2)  = RPM_cmd(2) + dRPM_ctrl;
        RPM_cmd(3)  = RPM_cmd(3) - dRPM_ctrl;
        RPM_cmd(4)  = RPM_cmd(4) - dRPM_ctrl;
        RPM_cmd     = max(RPM_min, min(RPM_max, RPM_cmd));

        s8 = state(1:8);
        k1 = eom_full(s8,             RPM_cmd, p_run);
        k2 = eom_full(s8 + 0.5*dt*k1, RPM_cmd, p_run);
        k3 = eom_full(s8 + 0.5*dt*k2, RPM_cmd, p_run);
        k4 = eom_full(s8 + dt*k3,     RPM_cmd, p_run);
        state(1:8) = s8 + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    end

    run_conv(s) = true;
    diverged(s) = did_diverge;

    if ~did_diverge
        idx_post   = t_vec >= t_pulse;
        theta_post = all_theta(s, idx_post);
        theta_tgt  = rad2deg(theta_s);
        peak_overshoot(s) = max(abs(theta_post - theta_tgt));
        band       = 0.5;
        settled    = abs(theta_post - theta_tgt) < band;
        if any(settled)
            lo = find(~settled, 1, 'last');
            if isempty(lo)
                settling_time(s) = 0;
            else
                tp = t_vec(idx_post);
                settling_time(s) = tp(lo) - t_pulse;
            end
        end
        idx_ss      = t_vec >= (t_sim - 0.5);
        ss_error(s) = mean(all_theta(s, idx_ss)) - theta_tgt;
    else
        peak_overshoot(s) = Inf;
        settling_time(s)  = Inf;
    end
    mc_peak_q(s) = max(abs(all_q(s,:)));

    if use_mc && mod(s, 25) == 0
        n_div = sum(diverged(1:s));
        n_ok  = sum(run_conv(1:s) & ~diverged(1:s));
        eta   = (toc(tic_total) / s) * (N_runs - s);
        fprintf('  %3d/%d: %d recovered, %d diverged  | ETA %.0fmin\n', ...
                s, N_runs, n_ok, n_div, eta/60);
    end
end
n_total = sum(run_conv);
n_div   = sum(diverged);
n_ok    = n_total - n_div;
fprintf('Done in %.1f min. %d/%d completed: %d recovered, %d diverged.\n', ...
        toc(tic_total)/60, n_total, N_runs, n_ok, n_div);

% Save ---------------------------------------------------------------------
if ~isfolder('results'); mkdir('results'); end
save_vars = {'all_theta', 'all_q', 'all_V', 'all_phi', 't_vec', ...
             'run_conv', 'diverged', ...
             'peak_overshoot', 'settling_time', 'ss_error', ...
             'mc_trim_theta', 'mc_trim_RPM', 'mc_peak_q', ...
             'V_trim', 'RUN_MODE', 'N_MC', 'RNG_SEED', ...
             'dRPM_pulse', 't_pulse_dur', 't_pulse', 't_sim', ...
             'Kp_angle', 'Kp_rate', 'Ki_rate', 'Kd_rate', 'f_D_cutoff', ...
             'theta_diverge'};
if use_mc
    save_vars = [save_vars, {'sigma_CL_grid', 'sigma_CD_grid', 'sigma_twist', ...
                             'z_CL_all', 'z_CD_all', 'z_twist_all'}];
end
save(out_file, save_vars{:});
fprintf('Saved: %s\n', out_file);


%  LOCAL FUNCTIONS ===========================================================================

function dx = eom_full(state, RPM, p)
% Full nonlinear 8-state equations of motion (body axes).
    global LAMBDA_CACHE  %#ok<GVMIS>
    if isempty(LAMBDA_CACHE), LAMBDA_CACHE = 0.08*ones(4,1); end

    u   = state(1); v  = state(2); w  = state(3);
    pb  = state(4); qb = state(5); rb = state(6);
    phi = state(7); theta = state(8);

    Vel     = [u; v; w];
    Omega_b = [pb; qb; rb];
    F_aero  = [0; 0; 0];
    M_aero  = [0; 0; 0];

    for i = 1:4
        r_i   = p.r_cg(:,i);
        v_hub = Vel + cross(Omega_b, r_i);
        Vz    = v_hub(3);
        Vxy   = v_hub(1:2);
        V0    = norm(Vxy);
        Omega = abs(2*pi*RPM(i)/60);
        mu_x  = V0 / (Omega*p.R + 1e-10);

        if mu_x < p.mu_switch
            Vi = solve_vi_hover(RPM(i), i, p);
            [T, Q] = bemt_hover(Vi, RPM(i), i, p);
            Fx = 0; Fy = 0;
        else
            [Vi, lam] = solve_vi_forward(RPM(i), V0, Vz, p.spin(i), i, ...
                                         p, LAMBDA_CACHE(i));
            LAMBDA_CACHE(i) = lam;
            [T, Q, Fx, Fy] = bemt_forward(Vi, RPM(i), V0, Vz, ...
                                          p.spin(i), i, p);
        end

        if V0 > 0.05 && mu_x >= p.mu_switch
            dir  = Vxy/V0;
            F_xy = -Fx*dir + Fy*[-dir(2); dir(1)];
        else
            F_xy = [0; 0];
        end
        F_r    = [F_xy; -T];
        F_aero = F_aero + F_r;
        M_aero = M_aero + cross(r_i, F_r) + [0; 0; -p.spin(i)*Q];
    end

    V_h   = Vel(1:2);
    V_mag = norm(V_h);
    if V_mag > 0.1
        q_dyn  = (V_mag/p.Vref_fuse)^2;
        F_aero = F_aero + [-p.Dref*q_dyn*V_h/V_mag; -p.Lref*q_dyn];
        M_aero = M_aero + [0; p.Mref*q_dyn; 0];
    end

    g_body = [-p.g*sin(theta);
               p.g*cos(theta)*sin(phi);
               p.g*cos(theta)*cos(phi)];

    udot = F_aero(1)/p.m + g_body(1) + rb*v - qb*w;
    vdot = F_aero(2)/p.m + g_body(2) + pb*w - rb*u;
    wdot = F_aero(3)/p.m + g_body(3) + qb*u - pb*v;

    Iomega    = p.I*Omega_b;
    omega_dot = p.Iinv*(M_aero - cross(Omega_b, Iomega));

    phi_dot   = pb + (qb*sin(phi) + rb*cos(phi))*tan(theta);
    theta_dot = qb*cos(phi) - rb*sin(phi);

    dx = [udot; vdot; wdot; omega_dot; phi_dot; theta_dot];
end


function RPM_hover = find_hover_rpm(p)
% Hover RPM (4 x rotor thrust = vehicle weight).
    W        = p.m * p.g;
    RPM_test = 8000:100:16000;
    T_test   = zeros(size(RPM_test));
    for i = 1:numel(RPM_test)
        Vi        = solve_vi_hover(RPM_test(i), 1, p);
        T_test(i) = 4 * bemt_hover(Vi, RPM_test(i), 1, p);
    end
    idx = find(T_test >= W, 1);
    if ~isempty(idx) && idx > 1
        RPM_hover = round(interp1(T_test(idx-1:idx), RPM_test(idx-1:idx), W));
    else
        RPM_hover = 10484;
    end
end


function guesses = build_guesses(x_prev, vi, nom_conv, nom_theta, nom_phi, ...
                                 nom_RPM, RPM_hover)
% Initial guess list: previous trim, recent trims, fan of fallbacks.
    guesses = {x_prev};
    for back = 1:min(3, vi-1)
        vj = vi - back;
        if nom_conv(vj)
            guesses{end+1} = [deg2rad(nom_theta(vj));
                              deg2rad(nom_phi(vj));
                              nom_RPM(vj,:)'];
        end
    end
    for th = [-5, -10, -15, -20, -25, -30, -35, -40]
        guesses{end+1} = [deg2rad(th); 0; ...
                          RPM_hover; RPM_hover; RPM_hover; RPM_hover];
        rf = max(2000,  RPM_hover - 200*abs(th));
        rr = min(21000, RPM_hover + 100*abs(th));
        guesses{end+1} = [deg2rad(th); 0; rf; rf; rr; rr];
    end
end


function [x_best, conv] = try_guesses(guesses, V, p, rpm_bounds)
% Run trim solver from each guess, keep best solution.
    best_res = Inf;
    x_best   = guesses{1};
    conv     = false;
    for ig = 1:numel(guesses)
        [x_try, c, ~, res] = trim_solve(guesses{ig}, V, p, rpm_bounds);
        if c && res < best_res
            x_best = x_try; best_res = res; conv = true; break;
        elseif res < best_res
            x_best = x_try; best_res = res;
        end
    end
end


function [x, converged, n_iter, res] = trim_solve(x0, V, p, rpm_bounds)
% Damped Newton-Raphson trim solver with backtracking line search.
    x         = x0;
    dx        = [1e-5; 1e-5; 0.5; 0.5; 0.5; 0.5];
    alpha     = 0.6;
    converged = false;

    for iter = 1:50
        R   = trim_residuals(x, V, p);
        res = norm(R);
        if res < 1e-8
            converged = true; n_iter = iter; break;
        end
        if iter > 10 && res > 1e6, break; end

        J = zeros(6);
        for j = 1:6
            xp = x; xp(j) = xp(j) + dx(j);
            J(:,j) = (trim_residuals(xp, V, p) - R) / dx(j);
        end
        step  = -J \ R;
        x_new = x + alpha * step;
        x_new(3:6) = max(rpm_bounds(1), min(rpm_bounds(2), x_new(3:6)));
        R_new = trim_residuals(x_new, V, p);

        for bt = 1:3
            if norm(R_new) < res, break; end
            x_new = x + (alpha * 0.5^bt) * step;
            x_new(3:6) = max(rpm_bounds(1), min(rpm_bounds(2), x_new(3:6)));
            R_new = trim_residuals(x_new, V, p);
        end
        x = x_new;
    end
    if ~converged, n_iter = 50; end
    res = norm(trim_residuals(x, V, p));
    if res < 1e-3, converged = true; end
end


function R = trim_residuals(x, V, p)
% Residual vector: body-axis force and moment.
    [F, M, ~, ~] = flight_forces(x, V, p);
    R = [F; M];
end


function [F_total, M_total, T_each, Q_each] = flight_forces(x, V_fwd, p)
% Total body-axis force and moment from rotors, fuselage, gravity at trim.
    theta_t = x(1); phi_t = x(2); rpm = x(3:6);
    Vel     = [V_fwd*cos(theta_t); 0; V_fwd*sin(theta_t)];
    Rates   = [0; 0; 0];
    F_b     = [0; 0; 0];
    M_b     = [0; 0; 0];
    T_each  = zeros(4,1);
    Q_each  = zeros(4,1);
    lam_cache = 0.08 * ones(4,1);

    for i = 1:4
        r_i   = p.r_cg(:,i);
        v_hub = Vel + cross(Rates, r_i);
        Vz    = v_hub(3);
        Vxy   = v_hub(1:2);
        V0    = norm(Vxy);
        Omega = abs(2*pi*rpm(i)/60);
        mu_x  = V0 / (Omega*p.R + 1e-10);

        if mu_x < p.mu_switch
            Vi = solve_vi_hover(rpm(i), i, p);
            [T, Q] = bemt_hover(Vi, rpm(i), i, p);
            Fx = 0; Fy = 0;
        else
            [Vi, lam] = solve_vi_forward(rpm(i), V0, Vz, p.spin(i), i, ...
                                         p, lam_cache(i));
            lam_cache(i) = lam;
            [T, Q, Fx, Fy] = bemt_forward(Vi, rpm(i), V0, Vz, ...
                                          p.spin(i), i, p);
        end
        T_each(i) = T;
        Q_each(i) = Q;

        if V0 > 0.05 && mu_x >= p.mu_switch
            dir  = Vxy/V0;
            F_xy = -Fx*dir + Fy*[-dir(2); dir(1)];
        else
            F_xy = [0; 0];
        end
        F_r = [F_xy; -T];
        F_b = F_b + F_r;
        M_b = M_b + cross(r_i, F_r) + [0; 0; -p.spin(i)*Q];
    end

    V_h   = Vel(1:2);
    V_mag = norm(V_h);
    if V_mag > 0.1
        q   = (V_mag/p.Vref_fuse)^2;
        F_b = F_b + [-p.Dref*q*V_h/V_mag; -p.Lref*q];
        M_b = M_b + [0; p.Mref*q; 0];
    end

    cph = cos(phi_t); sph = sin(phi_t);
    cth = cos(theta_t); sth = sin(theta_t);
    R_be = [cth,     sph*sth, cph*sth;
            0,       cph,    -sph;
            -sth,    sph*cth, cph*cth];
    F_grav  = R_be' * [0; 0; p.m*p.g];
    F_total = F_b + F_grav;
    M_total = M_b;
end


function [T, Q] = bemt_hover(Vi, RPM, rotor_idx, p)
% Hover BEMT (one rotor): integrate sectional loads at zero advance ratio.
    Omega = 2*pi*RPM/60;
    r     = p.r_blade;
    c     = p.c_blade;
    dr    = p.dr_blade;
    if isfield(p, 'theta_blade_r')
        theta = p.theta_blade_r{rotor_idx};
    else
        theta = p.theta_blade;
    end
    N  = length(r);
    dT = zeros(N,1);
    dQ = zeros(N,1);
    for i = 1:N
        Vt = Omega*r(i);
        VR = sqrt(Vi^2 + Vt^2);
        M  = VR/p.a_sound;
        ph = atan2(Vi, Vt);
        alpha_deg = (theta(i)-ph)*180/pi;
        [Cl, Cd] = c81_lookup(alpha_deg, M, rotor_idx, p);
        dL = 0.5*p.rho*VR^2*c(i)*Cl*dr(i);
        dD = 0.5*p.rho*VR^2*c(i)*Cd*dr(i);
        dT(i) = dL*cos(ph) - dD*sin(ph);
        dQ(i) = (dL*sin(ph) + dD*cos(ph))*r(i);
    end
    T = p.Nb*sum(dT);
    Q = p.Nb*sum(dQ);
end


function Vi = solve_vi_hover(RPM, rotor_idx, p)
% Newton-Raphson induced velocity from hover momentum balance.
    A  = pi*p.R^2;
    Vi = 5;
    for it = 1:30
        T  = bemt_hover(Vi, RPM, rotor_idx, p);
        f  = T - 2*p.rho*A*Vi^2;
        if abs(f) < 1e-6, break; end
        T2 = bemt_hover(Vi+0.01, RPM, rotor_idx, p);
        df = (T2 - 2*p.rho*A*(Vi+0.01)^2 - f) / 0.01;
        if abs(df) < 1e-10, break; end
        Vi = Vi - f/df;
        if Vi <= 0, Vi = 0.5; end
    end
end


function [T, Q, Fx, Fy] = bemt_forward(Vi, RPM, V0, Vz, spin, rotor_idx, p)
% Forward-flight BEMT (one rotor): azimuth-averaged sectional loads.
    Omega   = spin*2*pi*RPM/60;
    Vt_mat  = Omega*p.r_blade + V0*p.sin_psi;
    Va      = Vi + Vz;
    phi_mat = atan2(Va, abs(Vt_mat));
    VR2     = Va^2 + Vt_mat.^2;
    if isfield(p, 'theta_mat_r')
        th = p.theta_mat_r{rotor_idx};
    else
        th = p.theta_mat;
    end
    alpha_mat = (th - phi_mat) * (180/pi);
    M_mat     = sqrt(VR2)/p.a_sound;

    [Cl, Cd] = c81_lookup_matrix(alpha_mat, M_mat, rotor_idx, p);
    dL = 0.5*p.rho*VR2.*p.c_mat.*Cl.*p.dr_mat;
    dD = 0.5*p.rho*VR2.*p.c_mat.*Cd.*p.dr_mat;
    dT = dL.*cos(phi_mat) - dD.*sin(phi_mat);
    dH = dL.*sin(phi_mat) + dD.*cos(phi_mat);

    T   = p.Nb*mean(sum(dT,1));
    Q   = p.Nb*mean(sum(dH.*p.r_mat,1));
    dFx = -dH.*sign(Vt_mat).*p.neg_sin_psi_mat;
    dFy = -dH.*sign(Vt_mat).*p.cos_psi_mat;
    Fx  = p.Nb*mean(sum(dFx,1));
    Fy  = p.Nb*mean(sum(dFy,1));
end


function [Vi, lam] = solve_vi_forward(RPM, V0, Vz, spin, rotor_idx, p, lam0)
% Newton-Raphson Glauert inflow for forward flight.
    Omega = abs(2*pi*RPM/60);
    OmR   = Omega*p.R;
    A     = pi*p.R^2;
    mu_x  = V0/OmR;
    mu_z  = Vz/OmR;
    lam   = lam0;
    for it = 1:25
        [T,~,~,~] = bemt_forward(lam*OmR, RPM, V0, Vz, spin, rotor_idx, p);
        CT    = T/(p.rho*A*OmR^2);
        denom = 2*sqrt(mu_x^2 + (mu_z+lam)^2);
        if denom < 1e-5, denom = 1e-5; end
        f = lam - CT/denom;
        if abs(f) < 1e-8, break; end
        [T2,~,~,~] = bemt_forward((lam+1e-4)*OmR, RPM, V0, Vz, ...
                                  spin, rotor_idx, p);
        f2 = (lam+1e-4) - (T2/(p.rho*A*OmR^2)) / ...
                          (2*sqrt(mu_x^2+(mu_z+lam+1e-4)^2));
        df = (f2-f)/1e-4;
        if abs(df) < 1e-12, break; end
        lam = lam - f/df;
    end
    Vi = lam*OmR;
end


function [Cl, Cd] = c81_lookup(a, m, rotor_idx, p)
% Bilinear interpolation in (alpha, M) for one (a, m) point.
    CL_tbl = p.rotor_CL{rotor_idx};
    CD_tbl = p.rotor_CD{rotor_idx};
    a = max(p.CL_alpha_min, min(a, p.CL_alpha_max));
    m = max(p.CL_M_min,     min(m, p.CL_M_max));
    ia = discretize(a, p.CL_alpha_vec);
    ia(ia >= p.CL_n_alpha) = p.CL_n_alpha - 1;
    if isnan(ia), ia = 1; end
    fa = (a - p.CL_alpha_vec(ia)) / p.CL_d_alpha(ia);
    im = discretize(m, p.CL_M_vec);
    im(im >= p.CL_n_M) = p.CL_n_M - 1;
    if isnan(im), im = 1; end
    fm = (m - p.CL_M_vec(im)) / p.CL_d_M(im);
    Cl = (1-fa)*(1-fm)*CL_tbl(ia,im)   + fa*(1-fm)*CL_tbl(ia+1,im) + ...
         (1-fa)*fm    *CL_tbl(ia,im+1) + fa*fm    *CL_tbl(ia+1,im+1);
    Cd = (1-fa)*(1-fm)*CD_tbl(ia,im)   + fa*(1-fm)*CD_tbl(ia+1,im) + ...
         (1-fa)*fm    *CD_tbl(ia,im+1) + fa*fm    *CD_tbl(ia+1,im+1);
end


function [Cl, Cd] = c81_lookup_matrix(a, m, rotor_idx, p)
% Vectorised bilinear C81 lookup over a grid of (alpha, M) pairs.
    CL_tbl = p.rotor_CL{rotor_idx};
    CD_tbl = p.rotor_CD{rotor_idx};
    a = max(p.CL_alpha_min, min(a, p.CL_alpha_max));
    m = max(p.CL_M_min,     min(m, p.CL_M_max));
    ia = discretize(a, p.CL_alpha_vec);
    ia(ia >= p.CL_n_alpha) = p.CL_n_alpha - 1;
    ia(isnan(ia)) = 1;
    fa = (a - p.CL_alpha_vec(ia)) ./ p.CL_d_alpha(ia);
    im = discretize(m, p.CL_M_vec);
    im(im >= p.CL_n_M) = p.CL_n_M - 1;
    im(isnan(im)) = 1;
    fm = (m - p.CL_M_vec(im)) ./ p.CL_d_M(im);
    n  = p.CL_n_alpha;
    i00 = ia + (im-1)*n;
    i10 = ia + 1 + (im-1)*n;
    i01 = ia + im*n;
    i11 = ia + 1 + im*n;
    Cl = (1-fa).*(1-fm).*CL_tbl(i00) + fa.*(1-fm).*CL_tbl(i10) + ...
         (1-fa).*fm    .*CL_tbl(i01) + fa.*fm    .*CL_tbl(i11);
    Cd = (1-fa).*(1-fm).*CD_tbl(i00) + fa.*(1-fm).*CD_tbl(i10) + ...
         (1-fa).*fm    .*CD_tbl(i01) + fa.*fm    .*CD_tbl(i11);
end