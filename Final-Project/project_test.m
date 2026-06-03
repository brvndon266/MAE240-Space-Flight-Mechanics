%% GEO Satellite Perturbation Analysis using Gauss Variational Equations
% MAE 240 Project
%
% Four simulation cases:
%   Case 1: Keplerian (unperturbed control)
%   Case 2: Keplerian + Lunisolar Third-Body Gravity
%   Case 3: Keplerian + Solar Radiation Pressure (SRP)
%   Case 4: Keplerian + Lunisolar Gravity + SRP
%
% State vector: x = [a, e, i, Omega, omega, f]
% Units: km, s, radians

clear; close all; clc;

%% -----------------------------------------------------------------------
%% CONSTANTS & PARAMETERS
%% -----------------------------------------------------------------------

% Earth
mu_E = 398600.4418;         % [km^3/s^2]

% Moon
mu_M = 4902.800;            % [km^3/s^2]
a_M  = 384400.0;            % Semi-major axis [km]
T_M  = 27.3217 * 86400;     % Sidereal period [s]
n_M  = 2*pi / T_M;          % Mean motion [rad/s]
i_M  = deg2rad(5.145);      % Inclination to ecliptic [rad]

% Sun
mu_S  = 1.32712440018e11;   % [km^3/s^2]
a_S   = 1.495978707e8;      % Earth-Sun distance [km]
T_S   = 365.25 * 86400;     % Period [s]
n_S   = 2*pi / T_S;         % Mean motion [rad/s]

% SRP
% Acceleration magnitude: a_SRP = (P_SR * C_R * A/m)
% P_SR = 4.56e-6 N/m^2 at 1 AU; typical C_R*A/m ~ 0.02 m^2/kg
% => a_SRP ~ 9.1e-8 m/s^2 = 9.1e-11 km/s^2
% Use 1e-10 km/s^2 as a conservative placeholder — adjust A/m as needed
aSRP_mag = 1e-10;           % [km/s^2]

%% -----------------------------------------------------------------------
%% INITIAL GEO ORBIT
%% -----------------------------------------------------------------------

a0   = 42164.0;             % [km]
e0   = 0.001;               % Small but nonzero to avoid singularity
i0   = deg2rad(0.1);        % Near-equatorial [rad]
Om0  = deg2rad(0.0);        % RAAN [rad]
w0   = deg2rad(0.0);        % Argument of perigee [rad]
f0   = deg2rad(0.0);        % True anomaly [rad]

x0 = [a0; e0; i0; Om0; w0; f0];

%% -----------------------------------------------------------------------
%% SIMULATION TIME  (1 year)
%% -----------------------------------------------------------------------

T_GEO  = 2*pi * sqrt(a0^3 / mu_E);   % GEO period [s]
nYear= 25;
nRevs  = 365 * nYear;
tspan  = [0, nRevs * T_GEO];

opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

%% -----------------------------------------------------------------------
%% RUN ALL FOUR CASES
%% -----------------------------------------------------------------------

fprintf('Running Case 1: Keplerian (control)...\n');
[t1, X1] = ode113(@(t,x) GVE_full(t, x, mu_E, ...
    mu_M, a_M, n_M, i_M, ...
    mu_S, a_S, n_S, ...
    aSRP_mag, ...
    false, false, false), tspan, x0, opts);

fprintf('Running Case 2: Keplerian + Lunisolar Gravity...\n');
[t2, X2] = ode113(@(t,x) GVE_full(t, x, mu_E, ...
    mu_M, a_M, n_M, i_M, ...
    mu_S, a_S, n_S, ...
    aSRP_mag, ...
    true, false, false), tspan, x0, opts);

fprintf('Running Case 3: Keplerian + SRP...\n');
[t3, X3] = ode113(@(t,x) GVE_full(t, x, mu_E, ...
    mu_M, a_M, n_M, i_M, ...
    mu_S, a_S, n_S, ...
    aSRP_mag, ...
    false, true, false), tspan, x0, opts);

fprintf('Running Case 4: Keplerian + Lunisolar Gravity + SRP...\n');
[t4, X4] = ode113(@(t,x) GVE_full(t, x, mu_E, ...
    mu_M, a_M, n_M, i_M, ...
    mu_S, a_S, n_S, ...
    aSRP_mag, ...
    true, true, false), tspan, x0, opts);

fprintf('All cases complete.\n');

%% -----------------------------------------------------------------------
%% POST-PROCESS: convert to days, unwrap angles
%% -----------------------------------------------------------------------

cases = {X1, X2, X3, X4};
times = {t1, t2, t3, t4};
labels = {'Case 1: Keplerian', ...
          'Case 2: + Lunisolar Grav', ...
          'Case 3: + SRP', ...
          'Case 4: + Lunisolar + SRP'};
colors = {'k', 'b', 'r', [0 0.5 0]};

%% -----------------------------------------------------------------------
%% FIGURE 1: Orbital Elements vs Time (all cases overlaid per element)
%% -----------------------------------------------------------------------

elemNames = {'a [km]', 'e', 'i [deg]', '\Omega [deg]', '\omega [deg]', 'f [deg]'};
figure('Name','Orbital Elements','Position',[100 100 1400 900]);

for el = 1:6
    subplot(3,2,el); hold on; grid on;
    for c = 1:4
        X  = cases{c};
        td = times{c} / 86400;
        switch el
            case 1
                y = X(:,1);
            case 2
                y = X(:,2);
            case 3
                y = rad2deg(X(:,3));
            case 4
                y = rad2deg(unwrap(X(:,4)));
            case 5
                y = rad2deg(unwrap(X(:,5)));
            case 6
                y = rad2deg(unwrap(X(:,6)));
        end
        plot(td, y, ...
    'Color', colors{c}, ...
    'LineWidth', 1.2, ...
    'DisplayName', labels{c});
    end
    xlabel('Time [days]');
    ylabel(elemNames{el});
    title(elemNames{el});
    if el == 1
        legend('Location','best','FontSize',7);
    end
end
sgtitle('GEO Orbital Elements — All Cases');

%% -----------------------------------------------------------------------
%% FIGURE 2: Delta elements (deviation from Keplerian)
%% -----------------------------------------------------------------------

figure('Name','Element Deviations from Keplerian','Position',[100 100 1400 900]);

devLabels = {'\Delta a [km]', '\Delta e', '\Delta i [deg]', ...
             '\Delta\Omega [deg]', '\Delta\omega [deg]', '\Delta f [deg]'};

for el = 1:6
    subplot(3,2,el); hold on; grid on;
    for c = 2:4   % skip Case 1 (it's zero by definition)
        X  = cases{c};
        X1i = cases{1};
        td = times{c} / 86400;
        % interpolate Keplerian baseline onto this case's time vector
        switch el
            case 1
                y  = X(:,1);
                y1 = interp1(times{1}, X1i(:,1), times{c});
                dy = y - y1;
            case 2
                y  = X(:,2);
                y1 = interp1(times{1}, X1i(:,2), times{c});
                dy = y - y1;
            case 3
                y  = rad2deg(X(:,3));
                y1 = rad2deg(interp1(times{1}, X1i(:,3), times{c}));
                dy = y - y1;
            case 4
                y  = rad2deg(unwrap(X(:,4)));
                y1 = rad2deg(unwrap(interp1(times{1}, X1i(:,4), times{c})));
                dy = y - y1;
            case 5
                y  = rad2deg(unwrap(X(:,5)));
                y1 = rad2deg(unwrap(interp1(times{1}, X1i(:,5), times{c})));
                dy = y - y1;
            case 6
                y  = rad2deg(unwrap(X(:,6)));
                y1 = rad2deg(unwrap(interp1(times{1}, X1i(:,6), times{c})));
                dy = y - y1;
        end
        plot(td, dy, ...
    'Color', colors{c}, ...
    'LineWidth', 1.2, ...
    'DisplayName', labels{c});
    end
    xlabel('Time [days]');
    ylabel(devLabels{el});
    title(devLabels{el});
    if el == 1
        legend('Location','best','FontSize',7);
    end
end
sgtitle('GEO Element Deviations from Keplerian Baseline');

%% -----------------------------------------------------------------------
%% FIGURE 3: 3D Trajectories
%% -----------------------------------------------------------------------

figure('Name','3D Trajectories','Position',[100 100 1200 900]);
for c = 1:4
    X  = cases{c};
    Nt = size(X,1);
    rH = zeros(Nt,3);
    for k = 1:Nt
        [rI,~] = oe2rv(X(k,1),X(k,2),X(k,3),X(k,4),X(k,5),X(k,6),mu_E);
        rH(k,:) = rI';
    end
    subplot(2,2,c);
    plot3(rH(:,1),rH(:,2),rH(:,3),'Color',colors{c},'LineWidth',0.8);
    hold on;
    plot3(0,0,0,'ko','MarkerSize',8,'LineWidth',2);
    plot3(rH(1,1),rH(1,2),rH(1,3),'go','MarkerSize',6);
    plot3(rH(end,1),rH(end,2),rH(end,3),'r^','MarkerSize',6);
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    title(labels{c}); axis equal; grid on; view(35,25);
end
sgtitle('GEO 3D Trajectories');

fprintf('\nDone. Figures generated.\n');

%% -----------------------------------------------------------------------
%% PRINT START AND END STATES FOR ALL CASES
%% -----------------------------------------------------------------------

fprintf('\n========================================\n');
fprintf('       ORBITAL ELEMENT SUMMARY\n');
fprintf('========================================\n');

elemLabels = {'a [km]', 'e', 'i [deg]', 'Omega [deg]', 'omega [deg]', 'f [deg]'};

for c = 1:4
    X  = cases{c};
    td = times{c};
    
    fprintf('\n--- %s ---\n', labels{c});
    fprintf('%-15s %15s %15s %15s\n', 'Element', 'Start', 'End', 'Delta');
    fprintf('%s\n', repmat('-',1,60));
    
    for el = 1:6
        switch el
            case 1
                val_start = X(1,el);
                val_end   = X(end,el);
            case 2
                val_start = X(1,el);
                val_end   = X(end,el);
            otherwise
                val_start = rad2deg(X(1,el));
                val_end   = rad2deg(X(end,el));
        end
        delta = val_end - val_start;
        fprintf('%-15s %15.6f %15.6f %15.6f\n', ...
                elemLabels{el}, val_start, val_end, delta);
    end
end

%% -----------------------------------------------------------------------
%% SPECIFICALLY HIGHLIGHT SRP CASE (Case 3)
%% -----------------------------------------------------------------------

fprintf('\n========================================\n');
fprintf('   SRP PERTURBATION EFFECT (Case 3)\n');
fprintf('========================================\n');

X_kep = cases{1};   % Keplerian baseline
X_srp = cases{3};   % SRP case

% Eccentricity vector components
ex_start_kep = X_kep(1,2)   * cos(X_kep(1,5));
ey_start_kep = X_kep(1,2)   * sin(X_kep(1,5));
ex_end_kep   = X_kep(end,2) * cos(X_kep(end,5));
ey_end_kep   = X_kep(end,2) * sin(X_kep(end,5));

ex_start_srp = X_srp(1,2)   * cos(X_srp(1,5));
ey_start_srp = X_srp(1,2)   * sin(X_srp(1,5));
ex_end_srp   = X_srp(end,2) * cos(X_srp(end,5));
ey_end_srp   = X_srp(end,2) * sin(X_srp(end,5));

fprintf('\nEccentricity vector:\n');
fprintf('%-25s  ex = %+.6e   ey = %+.6e\n', 'Keplerian start:', ex_start_kep, ey_start_kep);
fprintf('%-25s  ex = %+.6e   ey = %+.6e\n', 'Keplerian end:',   ex_end_kep,   ey_end_kep);
fprintf('%-25s  ex = %+.6e   ey = %+.6e\n', 'SRP start:',       ex_start_srp, ey_start_srp);
fprintf('%-25s  ex = %+.6e   ey = %+.6e\n', 'SRP end:',         ex_end_srp,   ey_end_srp);

delta_ex = ex_end_srp - ex_end_kep;
delta_ey = ey_end_srp - ey_end_kep;
fprintf('\nSRP-induced eccentricity vector shift after 1 year:\n');
fprintf('  delta_ex = %+.6e\n', delta_ex);
fprintf('  delta_ey = %+.6e\n', delta_ey);
fprintf('  |delta_e| = %.6e\n', sqrt(delta_ex^2 + delta_ey^2));

% Convert to position displacement
a0 = X_srp(1,1);
delta_r_km = sqrt(delta_ex^2 + delta_ey^2) * a0;
fprintf('  Position displacement ~ %.4f km = %.1f m\n', delta_r_km, delta_r_km*1000);

% Inertial position start vs end
[r_start, ~] = oe2rv(X_srp(1,1),   X_srp(1,2),   X_srp(1,3), ...
                     X_srp(1,4),   X_srp(1,5),   X_srp(1,6),   mu_E);
[r_end,   ~] = oe2rv(X_srp(end,1), X_srp(end,2), X_srp(end,3), ...
                     X_srp(end,4), X_srp(end,5), X_srp(end,6), mu_E);
[r_kep_end,~] = oe2rv(X_kep(end,1), X_kep(end,2), X_kep(end,3), ...
                      X_kep(end,4), X_kep(end,5), X_kep(end,6), mu_E);

fprintf('\nInertial position (ECI):\n');
fprintf('%-25s  [%12.4f, %12.4f, %12.4f] km\n', 'SRP start:',     r_start(1), r_start(2), r_start(3));
fprintf('%-25s  [%12.4f, %12.4f, %12.4f] km\n', 'SRP end:',       r_end(1),   r_end(2),   r_end(3));
fprintf('%-25s  [%12.4f, %12.4f, %12.4f] km\n', 'Keplerian end:',  r_kep_end(1), r_kep_end(2), r_kep_end(3));

diff_r = r_end - r_kep_end;
fprintf('\nPosition difference (SRP end - Keplerian end):\n');
fprintf('  dX = %+.4f km\n', diff_r(1));
fprintf('  dY = %+.4f km\n', diff_r(2));
fprintf('  dZ = %+.4f km\n', diff_r(3));
fprintf('  |dr| = %.4f km = %.1f m\n', norm(diff_r), norm(diff_r)*1000);
%% =======================================================================
%% LOCAL FUNCTIONS
%% =======================================================================

function dx = GVE_full(t, x, mu_E, ...
                        mu_M, a_M, n_M, i_M, ...
                        mu_S, a_S, n_S, ...
                        aSRP_mag, ...
                        use_grav, use_SRP, ~)

%% Unpack state
a  = x(1);
e  = x(2);
i  = x(3);
Om = x(4);
w  = x(5);
f  = x(6);

% Protect against singularities
e = max(e, 1e-8);
i = max(i, deg2rad(0.01));

%% Orbit geometry
p     = a * (1 - e^2);
r     = p / (1 + e*cos(f));
h     = sqrt(mu_E * p);
eta   = sqrt(1 - e^2);
n_sat = sqrt(mu_E / a^3);

%% RTN frame unit vectors (inertial)
[rI, vI] = oe2rv(a, e, i, Om, w, f, mu_E);

er = rI / norm(rI);
eh = cross(rI, vI);  eh = eh / norm(eh);
et = cross(eh, er);

%% -----------------------------------------------------------------------
%% GRAVITY PERTURBATIONS (lunisolar third-body) via RTN projection
%% -----------------------------------------------------------------------

aI_grav = zeros(3,1);

if use_grav
    % Moon position
    theta_M = n_M * t;
    r_M = a_M * [cos(theta_M); ...
                  cos(i_M)*sin(theta_M); ...
                  sin(i_M)*sin(theta_M)];

    % Sun position
    eps_S   = deg2rad(23.45);
    theta_S = n_S * t;
    r_S = a_S * [cos(theta_S); ...
                  cos(eps_S)*sin(theta_S); ...
                  sin(eps_S)*sin(theta_S)];

    % Third-body accelerations
    r_rel_M = rI - r_M;
    a_moon  = -mu_M * (r_rel_M/norm(r_rel_M)^3 + r_M/norm(r_M)^3);

    r_rel_S = rI - r_S;
    a_sun   = -mu_S * (r_rel_S/norm(r_rel_S)^3 + r_S/norm(r_S)^3);

    aI_grav = a_moon + a_sun;
end

% Project gravity perturbation onto RTN
R = dot(aI_grav, er);
T = dot(aI_grav, et);
N = dot(aI_grav, eh);

%% -----------------------------------------------------------------------
%% GVE: gravity contribution
%% -----------------------------------------------------------------------

u = w + f;

da  = (2*a^2/h) * (e*sin(f)*R + (p/r)*T);
de  = (1/h) * (p*sin(f)*R + ((p+r)*cos(f) + r*e)*T);
di  = (r*cos(u)/h) * N;
dOm = (r*sin(u) / (h*sin(i))) * N;
dw  = (1/(h*e)) * (-p*cos(f)*R + (p+r)*sin(f)*T) ...
      - (r*sin(u)*cos(i) / (h*sin(i))) * N;
df  = h/r^2 + (1/(h*e)) * (p*cos(f)*R - (p+r)*sin(f)*T);

%% -----------------------------------------------------------------------
%% SRP PERTURBATION — Cook averaged model (Eq. 41-44 from Guffanti 2017)
%% -----------------------------------------------------------------------

if use_SRP
    eps_S    = deg2rad(23.45);
    theta_S  = n_S * t;
    lambda_s = atan2(cos(eps_S)*sin(theta_S), cos(theta_S));  % Sun ecliptic longitude

    % Half-angle trig
    ce2 = cos(eps_S/2);   se2 = sin(eps_S/2);
    ci  = cos(i);         si  = sin(i);
    ci2 = cos(i/2);       si2 = sin(i/2);

    % Angle combinations (omega + Omega +/- lambda_s, omega - Omega +/- lambda_s)
    opOm_m = w + Om - lambda_s;
    opOm_p = w + Om + lambda_s;
    omOm_p = w - Om + lambda_s;
    omOm_m = w - Om - lambda_s;

    Fsrp_Bsrp = aSRP_mag;   % [km/s^2]

    % --- Rp: radial component at perigee (Eq. 42) ---
    Rp =  Fsrp_Bsrp * ( ...
           ce2^2 * cos(opOm_m) * ci2^2 ...
         + se2^2 * cos(opOm_p) * ci2^2 ...
         + ce2^2 * cos(omOm_p) * si2^2 ...
         + se2^2 * cos(omOm_m) * si2^2 ...
         + 0.5*(cos(w - lambda_s) - cos(w + lambda_s)) * si * sin(eps_S) );

    % --- Sp: transverse component at perigee (Eq. 43) ---
    Sp = -Fsrp_Bsrp * ( ...
           ce2^2 * sin(opOm_m) * ci2^2 ...
         + se2^2 * sin(opOm_p) * ci2^2 ...
         + ce2^2 * sin(omOm_p) * si2^2 ...
         + se2^2 * sin(omOm_m) * si2^2 ...
         - 0.5*(sin(w + lambda_s) - sin(w - lambda_s)) * si * sin(eps_S) );

    % --- W: normal component — built from W*sin(w) and W*cos(w) (Eq. 44) ---
    Wsinw = -(Fsrp_Bsrp/2) * ( ...
              (cos(opOm_m) - cos(omOm_p)) * si * ce2^2 ...
            + (cos(opOm_p) - cos(omOm_m)) * si * se2^2 ...
            + (cos(w + lambda_s) - cos(w - lambda_s)) * ci * sin(eps_S) );

    Wcosw =  (Fsrp_Bsrp/2) * ( ...
              (sin(opOm_m) - sin(omOm_p)) * si * ce2^2 ...
            + (sin(opOm_p) - sin(omOm_m)) * si * se2^2 ...
            + (sin(w + lambda_s) - sin(w - lambda_s)) * ci * sin(eps_S) );

    W = Wsinw * sin(w) + Wcosw * cos(w);

    % --- Averaged element rates (Eq. 41) ---
    % da/dt = 0 for full sunlight

    de_srp  =  (3*eta) / (2*n_sat*a) * Sp;

    di_srp  = -(3*e*W*cos(w)) / (2*n_sat*a*eta);

    dOm_srp = -(3*e*W*sin(w)) / (2*n_sat*a*eta*sin(i));

    dw_srp  = -(3*eta) / (2*n_sat*a*e) * Rp  -  dOm_srp*cos(i);

    dM_srp  =  (9*e*Rp) / (2*n_sat*a) ...
               - eta * (dw_srp + dOm_srp*cos(i));

    % Convert mean anomaly rate to true anomaly rate
    % dM/dt -> df/dt:  df = dM * (1 + e*cos(f))^2 / eta^3
    df_srp  = dM_srp * (1 + e*cos(f))^2 / eta^3;

else
    de_srp  = 0;
    di_srp  = 0;
    dOm_srp = 0;
    dw_srp  = 0;
    df_srp  = 0;
end

%% -----------------------------------------------------------------------
%% TOTAL RATES
%% -----------------------------------------------------------------------

da  = da;
de  = de  + de_srp;
di  = di  + di_srp;
dOm = dOm + dOm_srp;
dw  = dw  + dw_srp;
df  = df  + df_srp;

dx = [da; de; di; dOm; dw; df];

end
% -------------------------------------------------------------------------

function [rI, vI] = oe2rv(a, e, i, Om, w, f, mu)
% OE2RV  Convert classical orbital elements to inertial position/velocity.

p = a * (1 - e^2);

rPQW = [p*cos(f)/(1+e*cos(f));
        p*sin(f)/(1+e*cos(f));
        0];

vPQW = sqrt(mu/p) * [-sin(f);
                       e + cos(f);
                       0];

R3Om = [ cos(Om) -sin(Om) 0;
         sin(Om)  cos(Om) 0;
         0        0       1];

R1i  = [1  0       0;
        0  cos(i) -sin(i);
        0  sin(i)  cos(i)];

R3w  = [ cos(w) -sin(w) 0;
         sin(w)  cos(w) 0;
         0       0      1];

Q  = R3Om * R1i * R3w;

rI = Q * rPQW;
vI = Q * vPQW;

end