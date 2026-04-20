function p = platform_parameters()
% PLATFORM_PARAMETERS (specific airframe)
%
% Returns a struct p containing all geometric, inertial, and aerodynamic
% reference data required by the BEMT and 6-DOF flight dynamics models.
%
% To apply the methodology to a different multi-rotor, edit the values in
% this file. 
%
% **No other script in the repository contains platform-specific data**
%
% EXISTING CONTENTS: QAV250 Holybro airframe with APC 5x4.5" three-bladed
% props.

% Environmental constants ----------------------------------------------------
p.g       = 9.81;        % gravitational acceleration [m/s^2]
p.rho     = 1.225;       % air density at sea level [kg/m^3]
p.a_sound = 340.3;       % speed of sound at sea level [m/s]

% Mass -----------------------------------------------------------------------
p.m_base  = 0.795925;    % airframe mass [kg]
p.m_pay   = 0.500;       % payload mass [kg] % set to zero for no added mass
p.m       = p.m_base + p.m_pay;

% Inertia tensor -------------------------------------------------------------
Ixx = 0.004164;  Iyy = 0.003518;  Izz = 0.007414;
Ixy = -4.882876e-6;  Ixz = +1.13585e-6;  Iyz = -4.05147e-7;
p.I    = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];
p.Iinv = inv(p.I);

% Rotor positions in body frame [m] ------------------------------------------
% Columns correspond to rotor indices 1-4 in the X-configuration.
%   Rotor 1: front-left   (CCW)
%   Rotor 2: front-right  (CW)
%   Rotor 3: rear-right   (CCW)
%   Rotor 4: rear-left    (CW)
p.r_cg = [ 0.076771 -0.098705 0.017442;
           0.076771  0.098705 0.017442;
          -0.077501 -0.098705 0.017442;
          -0.077501  0.098705 0.017442]';

p.spin = [-1; +1; +1; -1];   % rotor spin direction (sign of body-z angular velocity)

% Fuselage aerodynamic reference (from wind tunnel at V_ref) -----------------
p.Vref_fuse = 10;        % reference airspeed for fuselage data [m/s]
p.Dref      = 0.565;     % drag at V_ref [N]
p.Lref      = 0.742;     % lift at V_ref [N]
p.Mref      = 0.065;     % pitching moment at V_ref [N.m]

% Rotor geometry -------------------------------------------------------------
p.R  = 0.0635;           % blade radius [m]
p.Nb = 3;                % number of blades

% Radial stations [m]
p.r_blade = [8;12;16;20;24;28;30;32;34;36;38;40;42;44;46;48; ...
             50;52;54;56;57;58;59;60;61;62] / 1000;

% Local chord at each radial station [m]
p.c_blade = [10.934;11.943;12.771;13.394;13.797;13.972;13.972; ...
             13.913;13.794;13.614;13.375;13.079;12.728;12.324; ...
             11.866;11.354;10.788;10.160;9.453;8.685; ...
             8.220;7.661;7.057;6.381;5.325;3.219] / 1000;

% Geometric blade twist at each radial station [enter in deg]
p.theta_blade = deg2rad([28.9;27.8;26.8;26.0;25.1;24.2;23.7;23.1; ...
                         22.5;21.8;21.1;20.2;19.3;18.3;17.3;16.1; ...
                         14.8;13.5;12.1;10.8;10.3;9.7;9.5;9.9; ...
                         11.1;12.1]);

% Annular widths for blade-element integration [m]
Ns = length(p.r_blade);
x  = p.r_blade / p.R;
dx = zeros(Ns,1);
dx(1)        = x(2) - x(1);
dx(Ns)       = x(Ns) - x(Ns-1);
dx(2:Ns-1)   = 0.5 * (x(3:Ns) - x(1:Ns-2));
p.dr_blade   = p.R * dx;

% BEMT solver settings -------------------------------------------------------
p.mu_switch = 0.01;      % advance ratio threshold for hover to forward flight

% Azimuthal discretisation for forward flight --------------------------------
M_azi = 12; %change if desired
psi   = linspace(0, 2*pi, M_azi+1);
psi(end) = [];
p.M_azi             = M_azi;
p.sin_psi           = sin(psi);
p.cos_psi           = cos(psi);
p.c_mat             = repmat(p.c_blade,     1, M_azi);
p.dr_mat            = repmat(p.dr_blade,    1, M_azi);
p.r_mat             = repmat(p.r_blade,     1, M_azi);
p.theta_mat         = repmat(p.theta_blade, 1, M_azi);
p.neg_sin_psi_mat   = repmat(-p.sin_psi, Ns, 1);
p.cos_psi_mat       = repmat( p.cos_psi, Ns, 1);

end