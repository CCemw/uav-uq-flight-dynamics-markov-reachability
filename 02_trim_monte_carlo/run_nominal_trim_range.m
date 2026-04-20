% RUN_NOMINAL_TRIM_RANGE
% Deterministic trim sweep across a forward-airspeed range.
% Saves trim states and per-rotor loads to results/.
%
% Reads:  00_inputs/platform/platform_parameters.m
%         00_inputs/polars/afl{13,14,15}.c81.txt   (afl14 is the baseline for QAV250)
% Writes: results/nominal_trim_results_QAV250.mat >> used for results and
% plotting

clear; clc; close all;

% Configuration -------------------------------------------------------------
V_sweep       = 0:0.25:30; % adjust
airfoil_files = {'afl13.c81.txt', 'afl14.c81.txt', 'afl15.c81.txt'};
nominal_idx   = 2;
RPM_bounds    = [-50000, 50000]; % currently inf RPM limit
out_file      = 'results/nominal_trim_results_QAV250.mat';

% Setup ---------------------------------------------------------------------
p = platform_parameters();
p.mu_switch = 0.01;

for k = 1:numel(airfoil_files)
    [AF(k).CL, AF(k).CD, AF(k).alpha_CL, AF(k).alpha_CD, AF(k).M] = ...
        load_c81(airfoil_files{k});
end
nom_CL = AF(nominal_idx).CL;
nom_CD = AF(nominal_idx).CD;
for r = 1:4
    p.rotor_CL{r} = nom_CL;
    p.rotor_CD{r} = nom_CD;
end
p = setup_c81_lookup(p, AF(nominal_idx).alpha_CL, AF(nominal_idx).alpha_CD, ...
                     AF(nominal_idx).M);
p.RPM_hover = find_hover_rpm(p);

fprintf('Hover RPM: %d | Mass: %.3f kg\n', p.RPM_hover, p.m);

% Trim sweep ----------------------------------------------------------------
nV = numel(V_sweep);
nom_theta  = NaN(nV,1);
nom_phi    = NaN(nV,1);
nom_RPM    = NaN(nV,4);
nom_conv   = false(nV,1);
nom_thrust = NaN(nV,1);
nom_torque = NaN(nV,1);
nom_power  = NaN(nV,1);
nom_T_each = NaN(nV,4);
nom_Q_each = NaN(nV,4);

x_prev = [0; 0; p.RPM_hover; p.RPM_hover; p.RPM_hover; p.RPM_hover];

fprintf('(Deterministic) Trim Speed Sweep: %.2f to %.2f m/s, step %.2f m/s\n', ...
        V_sweep(1), V_sweep(end), V_sweep(2)-V_sweep(1), nV);
tic;
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

        [~, ~, T_each, Q_each] = flight_forces(x_best, V, p);
        nom_T_each(vi,:) = T_each';
        nom_Q_each(vi,:) = Q_each';
        nom_thrust(vi)   = sum(T_each);
        nom_torque(vi)   = sum(abs(Q_each));
        nom_power(vi)    = sum(abs(x_best(3:6)) .* (2*pi/60) .* abs(Q_each));

        fprintf('  V=%5.2f: th=%+6.2f  RPM=[%5.0f %5.0f %5.0f %5.0f]  T=%5.2fN  P=%6.1fW\n', ...
                V, nom_theta(vi), nom_RPM(vi,:), nom_thrust(vi), nom_power(vi));
    else
        fprintf('  V=%5.2f: FAILED\n', V);
    end
end
fprintf('Done in %.1fs. Trimmed %d/%d. Max trimmed speed: %.2f m/s.\n', ...
        toc, sum(nom_conv), nV, max(V_sweep(nom_conv)));

if ~isfolder('results'); mkdir('results'); end
save(out_file, ...
     'V_sweep', 'nom_theta', 'nom_phi', 'nom_RPM', 'nom_conv', ...
     'nom_thrust', 'nom_torque', 'nom_power', 'nom_T_each', 'nom_Q_each');
fprintf('Saved: %s\n', out_file);


%  LOCAL FUNCTIONS ===========================================================================

function [CL, CD, alpha_CL, alpha_CD, M_vec] = load_c81(filename)
% C81 table read: returns CL/CD as alpha-by-Mach matrices.
    A        = readmatrix(filename);
    M_vec    = A(1, 2:11);
    alpha_CL = A(2:100, 1);
    CL       = A(2:100, 2:11);
    alpha_CD = A(102:200, 1);
    CD       = A(102:200, 2:11);
end


function p = setup_c81_lookup(p, alpha_CL, alpha_CD, M_vec)
% Pre-compute breakpoints for fast bilinear interpolation.
    p.CL_alpha_vec = alpha_CL;
    p.CL_n_alpha   = length(alpha_CL);
    p.CL_M_vec     = M_vec;
    p.CL_n_M       = length(M_vec);
    p.CL_alpha_min = alpha_CL(1);
    p.CL_alpha_max = alpha_CL(end);
    p.CL_M_min     = M_vec(1);
    p.CL_M_max     = M_vec(end);
    p.CL_d_alpha   = diff(alpha_CL);
    p.CL_d_M       = diff(M_vec);

    p.CD_alpha_vec = alpha_CD;
    p.CD_n_alpha   = length(alpha_CD);
    p.CD_M_vec     = M_vec;
    p.CD_n_M       = length(M_vec);
    p.CD_alpha_min = alpha_CD(1);
    p.CD_alpha_max = alpha_CD(end);
    p.CD_M_min     = M_vec(1);
    p.CD_M_max     = M_vec(end);
    p.CD_d_alpha   = diff(alpha_CD);
    p.CD_d_M       = diff(M_vec);
end


function RPM_hover = find_hover_rpm(p)
% Hover RPM (4 x rotor thrust = vehicle weight.)
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
% Initial guess list for the trim solver: previous trim, recent trims, and
% a fan of attitude-RPM combinations spanning likely operating points (for qav250).
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
            x_best = x_try; best_res = res; conv = true;
            break;
        elseif res < best_res
            x_best = x_try; best_res = res;
        end
    end
end


function [x, converged, n_iter, res] = trim_solve(x0, V, p, rpm_bounds)
% DAMPED Newton-Raphson trim solver with backtracking line search.
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
        if iter > 10 && res > 1e6
            break;
        end

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
% Total body-axis forces and moment from rotors, fuselage, and gravity at trim.
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
    theta = p.theta_blade;
    N     = length(r);
    dT    = zeros(N,1);
    dQ    = zeros(N,1);
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
    alpha_mat = (p.theta_mat - phi_mat) * (180/pi);
    M_mat   = sqrt(VR2)/p.a_sound;

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