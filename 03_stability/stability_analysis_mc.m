% STABILITY_ANALYSIS_MC
% Produces:
%   results/stability_results_QAV250.mat
%     Nominal eigenvalues across the airspeed range
%     Nominal A matrices and stability derivatives (M_u, M_q, M_w, X_u)
%     Monte Carlo eigenvalue clouds at target speeds (CRN, seed=42)
%     Per-sample MC stability derivatives (M_u, M_q)
%
% Reads:  ../00_inputs/platform/platform_parameters.m
%         ../00_inputs/polars/afl{13,14,15}.c81.txt
% Writes: results/stability_results_QAV250.mat

clear; clc; close all;
global LAMBDA_CACHE
LAMBDA_CACHE = 0.08 * ones(4,1);

% Configuration -------------------------------------------------------------
V_sweep       = 0:1:26;                          % nominal eigenvalue sweep
V_mc_targets  = [5, 10, 15, 17, 20, 24];         % MC cloud speeds [m/s] (adjust)
N_MC          = 500;                             % MC samples per target
RNG_SEED      = 42;                              % CRN seed
sigma_twist   = deg2rad(1.0);                    % blade twist sigma [rad]
d2_n3         = 1.693;                           % Tippett constant for n=3
CD_floor      = 0.005;                           % minimum perturbed CD
RPM_bounds    = [-50000, 50000];
airfoil_files = {'afl13.c81.txt', 'afl14.c81.txt', 'afl15.c81.txt'};
nominal_idx   = 2;
out_file      = 'results/stability_results_QAV250.mat';

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

% Pointwise sigma grids from Range / d2 across the three polars
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

% Nominal trim sweep --------------------------------------------------------
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
        fprintf('  V=%5.2f: th=%+6.2f  RPM=[%5.0f %5.0f %5.0f %5.0f]\n', ...
                V, nom_theta(vi), nom_RPM(vi,:));
    else
        fprintf('  V=%5.2f: FAILED\n', V);
    end
end
fprintf('Nominal: %d/%d trimmed.\n', sum(nom_conv), nV);

% Nominal eigenvalue sweep --------------------------------------------------
state_names  = {'u','v','w','p','q','r','phi','theta'};
eig_real_nom = NaN(nV, 8);
eig_imag_nom = NaN(nV, 8);
A_all        = cell(nV, 1);
deriv_Mu     = NaN(nV, 1);
deriv_Mq     = NaN(nV, 1);
deriv_Mw     = NaN(nV, 1);
deriv_Xu     = NaN(nV, 1);

fprintf('(Nominal) Eigenvalue Sweep: %d trimmed speeds\n', sum(nom_conv));
for vi = 1:nV
    if ~nom_conv(vi), continue; end
    V         = V_sweep(vi);
    theta_rad = deg2rad(nom_theta(vi));
    phi_rad   = deg2rad(nom_phi(vi));
    RPM_trim  = nom_RPM(vi,:)';
    state0    = [V*cos(theta_rad); 0; V*sin(theta_rad); 0; 0; 0; phi_rad; theta_rad];

    A = compute_jacobian(state0, RPM_trim, p);
    A_all{vi} = A;

    ev = eig(A);
    [~, idx] = sort(real(ev), 'descend');
    ev = ev(idx);
    eig_real_nom(vi,:) = real(ev)';
    eig_imag_nom(vi,:) = imag(ev)';

    deriv_Mu(vi) = A(5,1);
    deriv_Mq(vi) = A(5,5);
    deriv_Mw(vi) = A(5,3);
    deriv_Xu(vi) = A(1,1);

    max_re = max(real(ev));
    if max_re > 0
        fprintf('  V=%5.2f: max Re(lambda)=%+.4f  t2=%.2fs  M_u=%+.4f  M_q=%+.4f\n', ...
                V, max_re, 0.693/max_re, deriv_Mu(vi), deriv_Mq(vi));
    else
        fprintf('  V=%5.2f: max Re(lambda)=%+.4f  (stable)\n', V, max_re);
    end
end

% CRN draws -----------------------------------------------------------------
fprintf('Pre-generating %d x 12 CRN draws (seed=%d).\n', N_MC, RNG_SEED);
rng(RNG_SEED, 'twister');
z_CL_all    = randn(4, N_MC);
z_CD_all    = randn(4, N_MC);
z_twist_all = randn(4, N_MC);

% Monte Carlo eigenvalue clouds ---------------------------------------------
nT           = numel(V_mc_targets);
mc_eig_real  = NaN(nT, N_MC, 8);
mc_eig_imag  = NaN(nT, N_MC, 8);
mc_eig_maxRe = NaN(nT, N_MC);
mc_eig_conv  = false(nT, N_MC);
mc_eig_theta = NaN(nT, N_MC);
mc_eig_RPM   = NaN(nT, N_MC, 4);
mc_deriv_Mu  = NaN(nT, N_MC);
mc_deriv_Mq  = NaN(nT, N_MC);

fprintf('(Monte Carlo) Eigenvalue Clouds: %d targets x %d samples = %d solves\n', ...
        nT, N_MC, nT*N_MC);
tic_total = tic;
for ti = 1:nT
    V      = V_mc_targets(ti);
    vi_nom = find(V_sweep == V & nom_conv', 1);
    if isempty(vi_nom)
        fprintf('  V=%5.2f: no nominal trim, skipping.\n', V);
        continue;
    end

    theta_nom = deg2rad(nom_theta(vi_nom));
    phi_nom   = deg2rad(nom_phi(vi_nom));
    RPM_nom   = nom_RPM(vi_nom,:)';

    tic_speed = tic;
    n_conv    = 0;

    for s = 1:N_MC
        % Apply CRN: same z values for sample s at every target speed
        p_loc = p;
        for r = 1:4
            p_loc.rotor_CL{r}      = nom_CL + z_CL_all(r,s) * sigma_CL_grid;
            p_loc.rotor_CD{r}      = max(nom_CD + z_CD_all(r,s) * sigma_CD_grid, CD_floor);
            twist_offset           = z_twist_all(r,s) * sigma_twist;
            p_loc.theta_blade_r{r} = p.theta_blade + twist_offset;
            p_loc.theta_mat_r{r}   = p.theta_mat   + twist_offset;
        end
        LAMBDA_CACHE = 0.08*ones(4,1);

        [x_mc, conv_mc] = try_guesses({[theta_nom; phi_nom; RPM_nom]}, ...
                                       V, p_loc, RPM_bounds);
        if ~conv_mc, continue; end

        mc_eig_conv(ti, s)    = true;
        mc_eig_theta(ti, s)   = rad2deg(x_mc(1));
        mc_eig_RPM(ti, s, :)  = x_mc(3:6)';

        state0 = [V*cos(x_mc(1)); 0; V*sin(x_mc(1)); 0; 0; 0; x_mc(2); x_mc(1)];
        A_mc   = compute_jacobian(state0, x_mc(3:6), p_loc);

        ev = eig(A_mc);
        [~, idx] = sort(real(ev), 'descend');
        ev = ev(idx);
        mc_eig_real(ti, s, :) = real(ev)';
        mc_eig_imag(ti, s, :) = imag(ev)';
        mc_eig_maxRe(ti, s)   = real(ev(1));
        mc_deriv_Mu(ti, s)    = A_mc(5,1);
        mc_deriv_Mq(ti, s)    = A_mc(5,5);

        n_conv = n_conv + 1;
    end

    eta = (toc(tic_total) / ti) * (nT - ti);
    fprintf('  V=%5.2f: %3d/%d converged  | %.0fs  ETA %.0fmin\n', ...
            V, n_conv, N_MC, toc(tic_speed), eta/60);
end
fprintf('Done in %.1f min.\n', toc(tic_total)/60);

% Save ----------------------------------------------------------------------
if ~isfolder('results'); mkdir('results'); end
save(out_file, ...
     'V_sweep', 'nV', 'state_names', ...
     'eig_real_nom', 'eig_imag_nom', 'A_all', ...
     'deriv_Mu', 'deriv_Mq', 'deriv_Mw', 'deriv_Xu', ...
     'nom_theta', 'nom_phi', 'nom_RPM', 'nom_conv', ...
     'V_mc_targets', 'N_MC', 'RNG_SEED', ...
     'mc_eig_real', 'mc_eig_imag', 'mc_eig_maxRe', ...
     'mc_eig_conv', 'mc_eig_theta', 'mc_eig_RPM', ...
     'mc_deriv_Mu', 'mc_deriv_Mq', ...
     'sigma_CL_grid', 'sigma_CD_grid', 'sigma_twist', ...
     'z_CL_all', 'z_CD_all', 'z_twist_all');
fprintf('Saved: %s\n', out_file);


%  LOCAL FUNCTIONS ===========================================================================

function dx = eom_full(state, RPM, p)
% Full nonlinear 8-state equations of motion (body axes).
    global LAMBDA_CACHE
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


function A = compute_jacobian(state0, RPM0, p)
% 8x8 system Jacobian via central differences about (state0, RPM0).
    delta = [0.1; 0.1; 0.1; 0.01; 0.01; 0.01; 0.005; 0.005];
    n = 8;
    A = zeros(n, n);
    for j = 1:n
        sp = state0; sp(j) = sp(j) + delta(j);
        sm = state0; sm(j) = sm(j) - delta(j);
        A(:,j) = (eom_full(sp, RPM0, p) - eom_full(sm, RPM0, p)) / (2*delta(j));
    end
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
    global LAMBDA_CACHE
    theta_t = x(1); phi_t = x(2); rpm = x(3:6);
    Vel     = [V_fwd*cos(theta_t); 0; V_fwd*sin(theta_t)];
    Rates   = [0; 0; 0];
    F_b     = [0; 0; 0];
    M_b     = [0; 0; 0];
    T_each  = zeros(4,1);
    Q_each  = zeros(4,1);

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
                                         p, LAMBDA_CACHE(i));
            LAMBDA_CACHE(i) = lam;
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