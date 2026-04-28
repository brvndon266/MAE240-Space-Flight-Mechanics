%{
FILENAME:  TBP_demo.m

COURSE:    MAE 240 -- Space Flight Mechaniics (Spring 2026)

DESCRIPTION:
    A sectioned (%%) standalone MATLAB demo script for in-class lecture on
    the Two-Body Problem (2BP). The script mirrors the style and pacing of
    NBP_demo.m, but specializes to the fully integrable Kepler problem.

    The lecture arc moves through three equivalent viewpoints:
      (i)   the full nonlinear 2BP in an inertial barycentric frame,
      (ii)  the exact reduction to relative/body-centric coordinates,
      (iii) the Keplerian one-body form in a fixed central field.

    Along the way the script highlights:
      * full inertial-frame dynamics of the two masses
      * exact reduction to relative motion r_ddot = -mu r / r^3
      * the body-centric specialization mu ~= G m1 when m2 << m1
      * 3-D orbital geometry and orbit-plane triads
      * the relationship among true, eccentric, and mean anomalies
      * conservation of specific angular momentum, specific energy,
        and eccentricity vector
      * round-trip conversion between Cartesian state and orbital elements
      * coordinate transformations among equatorial, ecliptic, and orbital
        (perifocal) frames
      * Earth-centered and Moon-centered test cases

    Available test cases:
      * 'circumterrestrial_cislunar'
      * 'circumterrestrial_translunar'
      * 'circumlunar_selenocentric'

    The circumterrestrial cislunar/translunar labels are used in the simple,
    geometric sense of geocentric radial scale. In this ideal 2BP demo there
    are no third-body perturbations, so the resulting motion is exactly
    Keplerian even at large distance.

INSTRUCTOR NOTES:
    - Run section-by-section (Ctrl+Enter) in lecture, or run the whole file.
    - All figures are guarded by doPlot.
    - Set HK.pauseAtSections = true and HK.useKeyboard = true to pause and
      inspect the workspace at K>> after each section (type dbcont to resume).
    - The script is standalone: all helper functions are included below.
%}

%% Section 0 -- Global controls (read first)
close all
clearvars
clc

HK.enabled         = true;
HK.doClc           = false;
HK.doClear         = false;
HK.doClose         = false;
HK.pauseAtSections = true;
HK.useKeyboard     = true;
HK.clcAfterSection = true;
HK.clearVarsKeep   = {'HK','fmtTitle','doPlot'};

doPlot = true;

format compact;

fmtTitle = defaultFmtTitle();
fmtTitle('2BP DEMO -- Two-Body Problem: inertial formulation, body-centric reduction, and Keplerian structure');

disp('Lecture arc:');
disp('  Sec. 1  Full inertial-frame 2BP');
disp('  Sec. 2  Constants and frame conventions');
disp('  Sec. 3  Test-case selection and orbital elements');
disp('  Sec. 4  Initial states in barycentric and body-centric form');
disp('  Sec. 5  Integrate full 2BP and reduced 1BP');
disp('  Sec. 6  Compare inertial, relative, and restricted body-centric forms');
disp('  Sec. 7  3-D orbit geometry');
disp('  Sec. 8  True / eccentric / mean anomaly geometry');
disp('  Sec. 9  Specific angular momentum vector');
disp('  Sec. 10 Specific energy');
disp('  Sec. 11 Eccentricity vector');
disp('  Sec. 12 Cartesian <-> orbital elements');
disp('  Sec. 13 Equatorial, ecliptic, and orbital-plane transformations');

hkPause(HK);

%% Section 1 -- The two-body problem: inertial and relative formulations
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 1 -- The two-body problem: inertial equations and exact relative reduction');

%{
  PHYSICS REVIEW
  ==============
  Let R1 and R2 denote the inertial position vectors of masses m1 and m2.
  Define the relative vector

      r = R2 - R1 .

  The exact inertial equations are

      R1_ddot = +G m2 (R2 - R1) / |R2 - R1|^3
      R2_ddot = -G m1 (R2 - R1) / |R2 - R1|^3 .

  Introduce the barycenter

      R_cm = (m1 R1 + m2 R2) / (m1 + m2) .

  Since the internal forces cancel pairwise, R_cm_ddot = 0, so the barycenter
  moves at constant velocity and may serve as the origin of an inertial frame.

  EXACT RELATIVE REDUCTION
  ========================
  Subtracting the two inertial equations gives

      r_ddot = -mu * r / |r|^3,

  where

      mu = G (m1 + m2) .

  This has exactly the same form as the familiar one-body problem in a fixed
  central field.  Thus the non-restricted 2BP is completely equivalent to a
  Kepler problem in the relative coordinate.

  BODY-CENTRIC SPECIALIZATION
  ===========================
  When m2 << m1, as for an artificial satellite or small spacecraft,

      mu = G (m1 + m2) ~= G m1,

  and the primary can be treated as fixed at the origin of a nonrotating,
  body-centric frame.
%}

disp('Key equations:');
disp('  R1_ddot = +G*m2*(R2-R1)/|R2-R1|^3');
disp('  R2_ddot = -G*m1*(R2-R1)/|R2-R1|^3');
disp('  r       = R2 - R1');
disp('  r_ddot  = -mu*r/|r|^3,     mu = G*(m1 + m2)');
disp('  if m2 << m1, then mu ~= G*m1 and the fixed-center 1BP emerges.');

hkPause(HK);

%% Section 2 -- Constants, units, and frame conventions
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 2 -- Constants, units, and frame conventions');

%{
  UNIT SYSTEM
  ===========
  We work in km-kg-s throughout.

  INERTIAL FRAMES USED IN THIS DEMO
  =================================
  N_eq  : equatorial inertial frame (ECI-like reference frame)
  N_ecl : ecliptic inertial frame, obtained from N_eq by a fixed obliquity
          rotation about the x-axis
  P     : perifocal/orbital frame {e_hat, q_hat, h_hat}

  For circumlunar motion we define a selenocentric inertial frame whose origin
  is at the Moon center and whose axes remain parallel to the chosen inertial
  axes.  Thus the mathematics is identical to the Earth-centered case, with
  Earth constants replaced by lunar ones.
%}

C.G         = 6.67430e-20;                % [km^3 / (kg s^2)]
C.R_E       = 6378.1366;                  % [km]
C.mu_E      = 3.986004354360959e5;        % [km^3 / s^2]
C.M_E       = C.mu_E / C.G;               % [kg]
C.R_M       = 1737.4;                     % [km]
C.mu_M      = 4.902800066163796e3;        % [km^3 / s^2]
C.M_M       = C.mu_M / C.G;               % [kg]
C.a_M       = 383397.7725;                % [km] mean Earth-Moon distance
C.eps_eq2ecl_deg = 23.439291111;          % [deg] mean obliquity
C.eps_eq2ecl     = deg2rad(C.eps_eq2ecl_deg);

fprintf('Earth: mu = %.9e km^3/s^2,  R = %.4f km\n', C.mu_E, C.R_E);
fprintf('Moon : mu = %.9e km^3/s^2,  R = %.4f km\n', C.mu_M, C.R_M);
fprintf('Earth-Moon mean distance a_M = %.4f km\n', C.a_M);
fprintf('Mean obliquity = %.6f deg\n', C.eps_eq2ecl_deg);

fprintf('\nAvailable demo cases:\n');
fprintf('  1. circumterrestrial_cislunar\n');
fprintf('  2. circumterrestrial_translunar\n');
fprintf('  3. circumlunar_selenocentric\n');

hkPause(HK);

%% Section 3 -- Test-case selection and orbital-element specification
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 3 -- Test-case selection and orbital-element specification');

%{
  CASE SELECTION
  ==============
  The demo is written for a single active case at a time. Change CASE.name to
  switch among the three presets below.

  The orbital elements are specified relative to the chosen body-centric frame:

      [a, e, i, Omega, omega, M0]

  with angles in degrees at this setup stage for readability.

  NOTES ON THE CIRCUMTERRESTRIAL LABELS
  =====================================
  The cislunar/translunar terminology here refers only to geocentric radial
  scale. In an exact 2BP those labels do not imply distinct dynamical regimes,
  because lunar and solar perturbations are absent.
%}

CASE.name = 'circumterrestrial_cislunar';
% CASE.name = 'circumterrestrial_translunar';
% CASE.name = 'circumlunar_selenocentric';

CASE = makeDemoCase(CASE.name, C);

fprintf('Selected case: %s\n', CASE.name);
fprintf('  Description: %s\n', CASE.description);
fprintf('  Primary     : %s\n', CASE.primary_name);
fprintf('  Secondary   : %s\n', CASE.secondary_name);
fprintf('  Frame       : %s\n', CASE.frame_name);
fprintf('  a           : %.3f km\n', CASE.koe_deg(1));
fprintf('  e           : %.6f\n', CASE.koe_deg(2));
fprintf('  i           : %.3f deg\n', CASE.koe_deg(3));
fprintf('  Omega       : %.3f deg\n', CASE.koe_deg(4));
fprintf('  omega       : %.3f deg\n', CASE.koe_deg(5));
fprintf('  M0          : %.3f deg\n', CASE.koe_deg(6));

fprintf('\nDerived radial scales for this case:\n');
fprintf('  rp = %.3f km\n', CASE.rp);
fprintf('  ra = %.3f km\n', CASE.ra);
if strcmp(CASE.primary_name, 'Earth')
    fprintf('  rp / a_M = %.4f\n', CASE.rp / C.a_M);
    fprintf('  ra / a_M = %.4f\n', CASE.ra / C.a_M);
end

hkPause(HK);

%% Section 4 -- Initial states in barycentric and body-centric form
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 4 -- Initial states: orbital elements, body-centric state, and barycentric inertial state');

%{
  FROM ORBIT ELEMENTS TO STATE VECTOR
  ===================================
  We first convert the chosen Keplerian elements to the relative Cartesian state
  in the body-centric inertial frame:

      x_rel = [r0; v0] .

  We then embed that exact relative state into the full barycentric inertial 2BP.
  With the barycenter chosen at the origin and at rest,

      R1 = -(m2 / (m1 + m2)) r
      R2 = +(m1 / (m1 + m2)) r
      V1 = -(m2 / (m1 + m2)) v
      V2 = +(m1 / (m1 + m2)) v .

  This is the exact mapping between the reduced 1BP and the full inertial 2BP.
%}

mu_exact      = CASE.mu_exact;
mu_restricted = CASE.mu_restricted;
Mtot          = CASE.m1 + CASE.m2;

[r_rel0, v_rel0, aux0] = koe2rv_local(CASE.koe_deg, mu_exact, 'deg');

R10 = -(CASE.m2 / Mtot) * r_rel0;
R20 = +(CASE.m1 / Mtot) * r_rel0;
V10 = -(CASE.m2 / Mtot) * v_rel0;
V20 = +(CASE.m1 / Mtot) * v_rel0;

x0_full          = [R10; R20; V10; V20];
x0_rel_exact     = [r_rel0; v_rel0];
x0_rel_restricted = [r_rel0; v_rel0];

fprintf('Initial relative state in %s frame:\n', CASE.frame_name);
fprintf('  r0 = [% .6e  % .6e  % .6e] km\n', r_rel0(1), r_rel0(2), r_rel0(3));
fprintf('  v0 = [% .6e  % .6e  % .6e] km/s\n', v_rel0(1), v_rel0(2), v_rel0(3));

fprintf('\nInitial barycentric inertial states:\n');
fprintf('  R1 = [% .6e  % .6e  % .6e] km\n', R10(1), R10(2), R10(3));
fprintf('  R2 = [% .6e  % .6e  % .6e] km\n', R20(1), R20(2), R20(3));
fprintf('  |R1| / |r| = %.6e\n', norm(R10)/norm(r_rel0));
fprintf('  |R2| / |r| = %.6e\n', norm(R20)/norm(r_rel0));

fprintf('\nExact vs restricted gravitational parameter:\n');
fprintf('  mu_exact      = %.12e km^3/s^2\n', mu_exact);
fprintf('  mu_restricted = %.12e km^3/s^2\n', mu_restricted);
fprintf('  relative difference = %.6e\n', abs(mu_exact - mu_restricted)/mu_exact);

hkPause(HK);

%% Section 5 -- Integrate the full 2BP and the reduced 1BP
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 5 -- Integrate the full inertial 2BP and the reduced body-centric equations');

%{
  NUMERICAL PLAN
  ==============
  We integrate three mathematically related systems:

    (1) Full inertial 2BP in barycentric coordinates
    (2) Exact reduced relative equation with mu = G (m1 + m2)
    (3) Restricted body-centric approximation with mu = G m1

  For a spacecraft around a planet or moon, (2) and (3) are nearly identical.
  The full 2BP and the exact relative system should agree to numerical precision.
%}

T_orbit = 2*pi * sqrt(CASE.koe_deg(1)^3 / mu_exact);
T_end   = CASE.nPeriods * T_orbit;
Nt      = CASE.nTimes;

tspan   = linspace(0, T_end, Nt);
opts    = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);

fprintf('Integrating %.2f orbital periods (T = %.6f hr) with %d output times...\n', ...
    CASE.nPeriods, T_orbit/3600, Nt);

tic
[T, X_full] = ode45(@(t,x) eom_2bp_full(t, x, CASE.m1, CASE.m2, C.G), tspan, x0_full, opts);
t_full = toc;

tic
[~, X_rel_exact] = ode45(@(t,x) eom_2bp_relative(t, x, mu_exact), tspan, x0_rel_exact, opts);
t_rel_exact = toc;

tic
[~, X_rel_restricted] = ode45(@(t,x) eom_2bp_relative(t, x, mu_restricted), tspan, x0_rel_restricted, opts);
t_rel_restricted = toc;

fprintf('  Full inertial 2BP       : %.3f s\n', t_full);
fprintf('  Exact reduced relative  : %.3f s\n', t_rel_exact);
fprintf('  Restricted body-centric : %.3f s\n', t_rel_restricted);

R1_hist = X_full(:, 1:3);
R2_hist = X_full(:, 4:6);
V1_hist = X_full(:, 7:9);
V2_hist = X_full(:,10:12);

r_rel_from_full = R2_hist - R1_hist;
v_rel_from_full = V2_hist - V1_hist;

r_rel_exact_hist = X_rel_exact(:,1:3);
v_rel_exact_hist = X_rel_exact(:,4:6);

r_rel_restr_hist = X_rel_restricted(:,1:3);
v_rel_restr_hist = X_rel_restricted(:,4:6);

hkPause(HK);

%% Section 6 -- Compare full inertial, exact relative, and restricted body-centric forms
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 6 -- Comparing the full inertial 2BP, the exact reduction, and the restricted body-centric 1BP');

%{
  EXPECTATION
  ===========
  * Full inertial 2BP  <-> exact relative 1BP : agreement to numerical precision
  * Exact relative 1BP <-> restricted body-centric 1BP : agreement when m2 << m1

  This section verifies both statements directly by differencing trajectories.
%}

pos_err_full_vs_exact = rownorm(r_rel_from_full - r_rel_exact_hist);
vel_err_full_vs_exact = rownorm(v_rel_from_full - v_rel_exact_hist);

pos_err_exact_vs_restr = rownorm(r_rel_exact_hist - r_rel_restr_hist);
vel_err_exact_vs_restr = rownorm(v_rel_exact_hist - v_rel_restr_hist);

fprintf('Full inertial vs exact reduced (should be numerical roundoff level):\n');
fprintf('  max |dr| = %.6e km\n', max(pos_err_full_vs_exact));
fprintf('  max |dv| = %.6e km/s\n', max(vel_err_full_vs_exact));

fprintf('\nExact reduced vs restricted body-centric:\n');
fprintf('  max |dr| = %.6e km\n', max(pos_err_exact_vs_restr));
fprintf('  max |dv| = %.6e km/s\n', max(vel_err_exact_vs_restr));

if doPlot
    figure('Color','w','Name','2BP -- formulation comparison');
    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    semilogy(T/3600, max(pos_err_full_vs_exact, eps), 'k-', 'LineWidth', 1.7); hold on;
    semilogy(T/3600, max(pos_err_exact_vs_restr, eps), 'r--', 'LineWidth', 1.7);
    grid on; box on;
    ylabel('|dr|  [km]');
    title('Relative-position differences among formulations');
    legend({'full inertial vs exact relative', 'exact relative vs restricted body-centric'}, ...
        'Location','best');

    nexttile;
    semilogy(T/3600, max(vel_err_full_vs_exact, eps), 'k-', 'LineWidth', 1.7); hold on;
    semilogy(T/3600, max(vel_err_exact_vs_restr, eps), 'r--', 'LineWidth', 1.7);
    grid on; box on;
    xlabel('Time  [hr]');
    ylabel('|dv|  [km/s]');
    title('Relative-velocity differences among formulations');
end

hkPause(HK);

%% Section 7 -- 3-D orbit geometry in inertial and body-centric frames
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 7 -- 3-D orbit geometry: barycentric inertial picture and body-centric relative picture');

%{
  WHAT TO LOOK FOR
  ================
  Left panel:
    Full inertial barycentric motion of both bodies.
    For a spacecraft case, the primary moves only slightly around the barycenter.

  Middle panel:
    Exact relative orbit in the body-centric frame. This is the standard Keplerian
    picture: the primary is at the origin and the secondary follows a conic section.

  Right panel:
    Projection into the orbital plane, with periapsis direction e_hat and orbit pole h_hat.
%}

h0_vec = cross(r_rel0, v_rel0);
h0_hat = h0_vec / norm(h0_vec);
e0_vec = cross(v_rel0, h0_vec) / mu_exact - r_rel0 / norm(r_rel0);
e0_hat = e0_vec / norm(e0_vec);
q0_hat = cross(h0_hat, e0_hat);

if doPlot
    figure('Color','w','Name','2BP -- orbit geometry', 'Position', [60 60 1350 460]);
    tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile;
    hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on'); axis(ax1,'equal');
    hR1 = plot3(ax1, R1_hist(:,1), R1_hist(:,2), R1_hist(:,3), 'b-', 'LineWidth', 1.8);
    hR2 = plot3(ax1, R2_hist(:,1), R2_hist(:,2), R2_hist(:,3), 'r-', 'LineWidth', 1.8);
    hBc = plot3(ax1, 0, 0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot3(ax1, R1_hist(1,1), R1_hist(1,2), R1_hist(1,3), 'bo', 'MarkerFaceColor', 'b', 'HandleVisibility','off');
    plot3(ax1, R2_hist(1,1), R2_hist(1,2), R2_hist(1,3), 'ro', 'MarkerFaceColor', 'r', 'HandleVisibility','off');
    xlabel(ax1,'x [km]'); ylabel(ax1,'y [km]'); zlabel(ax1,'z [km]');
    title(ax1,'Full inertial 2BP (barycentric frame)');
    legend(ax1,[hR1, hR2, hBc],{CASE.primary_name, CASE.secondary_name, 'barycenter'}, 'Location','best');
    view(ax1,3);

    ax2 = nexttile;
    hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on'); axis(ax2,'equal');
    plot3(ax2, r_rel_exact_hist(:,1), r_rel_exact_hist(:,2), r_rel_exact_hist(:,3), ...
        'r-', 'LineWidth', 1.8);
    plot3(ax2, 0, 0, 0, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    plot3(ax2, r_rel_exact_hist(1,1), r_rel_exact_hist(1,2), r_rel_exact_hist(1,3), ...
        'ro', 'MarkerFaceColor', 'r');
    drawBody(ax2, CASE.primary_radius);
    xlabel(ax2,'x [km]'); ylabel(ax2,'y [km]'); zlabel(ax2,'z [km]');
    title(ax2, sprintf('Exact relative motion in %s', CASE.frame_name));
    view(ax2,3);

    ax3 = nexttile;
    hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on'); axis(ax3,'equal');
    rP_hist = (CASE.Q_P2N.' * r_rel_exact_hist.').';
    plot(ax3, rP_hist(:,1), rP_hist(:,2), 'r-', 'LineWidth', 1.8);
    plot(ax3, 0, 0, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 8);
    quiver(ax3, 0, 0, CASE.rp*0.35, 0, 0, 'k', 'LineWidth', 1.4, 'MaxHeadSize', 0.7);
    quiver(ax3, 0, 0, 0, CASE.rp*0.35, 0, 'k--', 'LineWidth', 1.2, 'MaxHeadSize', 0.7);
    text(0.38*CASE.rp, 0, 'e-hat', 'Parent', ax3);
    text(0, 0.38*CASE.rp, 'q-hat', 'Parent', ax3);
    xlabel(ax3,'p-axis [km]'); ylabel(ax3,'q-axis [km]');
    title(ax3,'Orbital plane / perifocal frame');
end

fprintf('Initial orbit triad in inertial coordinates:\n');
fprintf('  e_hat = [% .6f  % .6f  % .6f]\n', e0_hat(1), e0_hat(2), e0_hat(3));
fprintf('  q_hat = [% .6f  % .6f  % .6f]\n', q0_hat(1), q0_hat(2), q0_hat(3));
fprintf('  h_hat = [% .6f  % .6f  % .6f]\n', h0_hat(1), h0_hat(2), h0_hat(3));

hkPause(HK);

%% Section 8 -- The true, eccentric, and mean anomalies
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 8 -- True, eccentric, and mean anomalies: geometry and time dependence');

%{
  ANOMALY GEOMETRY
  ================
  For an elliptic orbit:

      M = E - e sin(E)                (Kepler''s equation)
      tan(f/2) = sqrt((1+e)/(1-e)) tan(E/2)

  where
      f : true anomaly      (focus-based angle)
      E : eccentric anomaly (center-based auxiliary angle)
      M : mean anomaly      (uniform-in-time angle)

  The mapping E -> M is explicit. The inverse mapping M -> E is transcendental
  and must be solved numerically. The body''s physical position is located by f.
%}

% Evaluate anomalies over one orbit from the exact reduced solution.
N_anom = numel(T);
f_hist = zeros(N_anom,1);
E_hist = zeros(N_anom,1);
M_hist = zeros(N_anom,1);

for k = 1:N_anom
    kep_k    = rv2koe_local(r_rel_exact_hist(k,:).', v_rel_exact_hist(k,:).', mu_exact, 'rad');
    f_hist(k)= kep_k.nu;
    E_hist(k)= kep_k.E;
    M_hist(k)= kep_k.M;
end

E_grid = linspace(0, 2*pi, 600).';
f_grid = 2*atan2(sqrt(1 + CASE.e)*sin(E_grid/2), sqrt(1 - CASE.e)*cos(E_grid/2));
f_grid = wrapTo2Pi(f_grid);
M_grid = wrapTo2Pi(E_grid - CASE.e*sin(E_grid));

rP0 = CASE.Q_P2N.' * r_rel0;
pt_curr = rP0(1:2);                         % current point on ellipse P
E0 = aux0.E;
f0 = aux0.f;
M0 = aux0.M;
a  = CASE.a;
b  = CASE.a*sqrt(1 - CASE.e^2);
xc = -CASE.a*CASE.e;                        % ellipse center C in focus-centered coordinates

% Auxiliary-circle construction for the eccentric anomaly E:
%   ellipse point  P = [a(cosE - e), b sinE]
%   circle point   Q = [xc + a cosE, a sinE]
% The auxiliary circle is centered at the ellipse center C = (xc,0), not at the focus.
pt_aux = [xc + a*cos(E0); a*sin(E0)];       % point Q on auxiliary circle corresponding to E0

if doPlot
    figure('Color','w','Name','2BP -- anomalies', 'Position', [50 50 1200 780]);
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile;
    hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on'); axis(ax1,'equal');
    tt = linspace(0, 2*pi, 500);
    x_ref = xc + CASE.a*cos(tt);            % auxiliary circle centered at ellipse center C
    y_ref =      CASE.a*sin(tt);
    x_ell = CASE.a*(cos(tt) - CASE.e);
    y_ell = b*sin(tt);
    plot(ax1, x_ref, y_ref, 'k--', 'LineWidth', 1.1);
    plot(ax1, x_ell, y_ell, 'r-', 'LineWidth', 1.8);
    plot(ax1, 0, 0, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    plot(ax1, xc, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(ax1, pt_curr(1), pt_curr(2), 'ro', 'MarkerFaceColor', 'r');
    plot(ax1, pt_aux(1), pt_aux(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, ...
        'HandleVisibility','off');
    line(ax1, [0 pt_curr(1)], [0 pt_curr(2)], 'Color', [0.2 0.2 0.8], 'LineWidth', 1.3);
    line(ax1, [xc pt_aux(1)], [0 pt_aux(2)], 'Color', [0.1 0.6 0.1], 'LineWidth', 1.3);
    line(ax1, [pt_aux(1) pt_curr(1)], [pt_aux(2) pt_curr(2)], ...
        'Color', [0.1 0.1 0.1], 'LineStyle', '--', 'LineWidth', 1.1);
    xlabel(ax1,'p [km]'); ylabel(ax1,'q [km]');
    title(ax1,'Ellipse, auxiliary circle, and current point');
    legend(ax1, {'auxiliary circle', 'ellipse', CASE.primary_name, 'ellipse center', ...
        CASE.secondary_name, 'true-anomaly radius', 'eccentric-anomaly radius', ...
        'vertical projection (same x)'}, 'Location','northwest');

    ax2 = nexttile;
    hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
    plot(ax2, E_grid, M_grid, 'k-', 'LineWidth', 1.8);
    plot(ax2, E_grid, E_grid, 'r--', 'LineWidth', 1.2);
    plot(ax2, E0, M0, 'ko', 'MarkerFaceColor','k');
    xlabel(ax2,'E [rad]'); ylabel(ax2,'M [rad]');
    title(ax2,'Kepler equation: M = E - e sin E');
    legend(ax2, {'M(E)', 'M = E', 'initial state'}, 'Location','best');

    ax3 = nexttile;
    hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
    plot(ax3, T/3600, unwrap(M_hist), 'k-', 'LineWidth', 1.6);
    plot(ax3, T/3600, unwrap(E_hist), 'b--', 'LineWidth', 1.6);
    plot(ax3, T/3600, unwrap(f_hist), 'r-.', 'LineWidth', 1.6);
    xlabel(ax3,'Time [hr]'); ylabel(ax3,'Angle [rad]');
    title(ax3,'Anomaly evolution over time');
    legend(ax3, {'M(t)', 'E(t)', 'f(t)'}, 'Location','best');

    ax4 = nexttile;
    hold(ax4,'on'); grid(ax4,'on'); box(ax4,'on');
    plot(ax4, E_grid, f_grid, 'm-', 'LineWidth', 1.8);
    plot(ax4, E0, f0, 'ko', 'MarkerFaceColor','k');
    xlabel(ax4,'E [rad]'); ylabel(ax4,'f [rad]');
    title(ax4,'Geometric map from eccentric anomaly to true anomaly');
end

fprintf('Initial anomaly values:\n');
fprintf('  M0 = %.6f rad  (%.6f deg)\n', M0, rad2deg(M0));
fprintf('  E0 = %.6f rad  (%.6f deg)\n', E0, rad2deg(E0));
fprintf('  f0 = %.6f rad  (%.6f deg)\n', f0, rad2deg(f0));

hkPause(HK);

%% Section 9 -- First integral: specific angular momentum vector
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 9 -- First integral: specific angular momentum vector h = r x v');

%{
  THEORY
  ======
  For the Kepler problem the force is central, so

      h = r x v

  is constant in both magnitude and direction.  Hence all motion is planar and
  the orbit pole h_hat is fixed in inertial space.
%}

h_hist = cross(r_rel_exact_hist, v_rel_exact_hist, 2);
h0     = h_hist(1,:);
dh     = h_hist - h0;
rel_h  = rownorm(dh) / norm(h0);

fprintf('Angular momentum diagnostics:\n');
fprintf('  |h0| = %.12e km^2/s\n', norm(h0));
fprintf('  max relative drift in h = %.6e\n', max(rel_h));

if doPlot
    figure('Color','w','Name','2BP -- angular momentum');
    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    hold on;
    if abs(h0(1)) > eps*norm(h0)
        plot(T/3600, h_hist(:,1)/h0(1), 'r-', 'LineWidth', 1.4);
    else
        plot(T/3600, h_hist(:,1), 'r-', 'LineWidth', 1.4);
    end
    if abs(h0(2)) > eps*norm(h0)
        plot(T/3600, h_hist(:,2)/h0(2), 'g-', 'LineWidth', 1.4);
    else
        plot(T/3600, h_hist(:,2), 'g-', 'LineWidth', 1.4);
    end
    plot(T/3600, h_hist(:,3)/h0(3), 'b-', 'LineWidth', 1.6);
    grid on; box on;
    ylabel('component history');
    title('Specific angular momentum components');
    legend({'h_x', 'h_y', 'h_z'}, 'Location','best');

    nexttile;
    semilogy(T/3600, max(rel_h, eps), 'k-', 'LineWidth', 1.6);
    grid on; box on;
    xlabel('Time [hr]'); ylabel('relative drift');
    title('Relative drift in the angular momentum vector');
end

hkPause(HK);

%% Section 10 -- First integral: specific mechanical energy
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 10 -- First integral: specific mechanical energy');

%{
  THEORY
  ======
  The relative Kepler problem conserves the specific mechanical energy

      E = 0.5 |v|^2 - mu / |r| .

  For an ellipse, E < 0 and a = -mu / (2E).
%}

rmag_hist = rownorm(r_rel_exact_hist);
vmag_hist = rownorm(v_rel_exact_hist);
energy_hist = 0.5*vmag_hist.^2 - mu_exact ./ rmag_hist;
energy0     = energy_hist(1);
dE          = energy_hist - energy0;

fprintf('Energy diagnostics:\n');
fprintf('  E0 = %.12e km^2/s^2\n', energy0);
fprintf('  max |Delta E| = %.6e km^2/s^2\n', max(abs(dE)));
fprintf('  a from energy = %.6f km\n', -mu_exact/(2*energy0));

if doPlot
    figure('Color','w','Name','2BP -- energy');
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(T/3600, 0.5*vmag_hist.^2, 'b-', 'LineWidth', 1.5); hold on;
    plot(T/3600, -mu_exact./rmag_hist, 'r-', 'LineWidth', 1.5);
    grid on; box on;
    ylabel('Specific energy terms');
    title('Specific kinetic and potential energy');
    legend({'0.5 v^2', '-mu/r'}, 'Location','best');

    nexttile;
    plot(T/3600, energy_hist, 'k-', 'LineWidth', 1.6);
    grid on; box on;
    ylabel('E [km^2/s^2]');
    title('Total specific energy');

    nexttile;
    semilogy(T/3600, max(abs(dE), eps), 'k-', 'LineWidth', 1.6);
    grid on; box on;
    xlabel('Time [hr]');
    ylabel('|Delta E|');
    title('Energy conservation error');
end

hkPause(HK);

%% Section 11 -- First integral: eccentricity vector
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 11 -- First integral: eccentricity vector');

%{
  THEORY
  ======
  The Kepler problem admits the additional vector invariant

      e = (v x h)/mu - r_hat ,

  which lies in the orbit plane, points to periapsis, and has magnitude equal
  to the scalar eccentricity. Together, h, e, and E encode the size, shape,
  and inertial orientation of the orbit.
%}

evec_hist = zeros(size(r_rel_exact_hist));
for k = 1:size(r_rel_exact_hist,1)
    rk = r_rel_exact_hist(k,:).';
    vk = v_rel_exact_hist(k,:).';
    hk = cross(rk, vk);
    evec_hist(k,:) = (cross(vk, hk)/mu_exact - rk/norm(rk)).';
end

e0 = evec_hist(1,:);
de = evec_hist - e0;
rel_e = rownorm(de) / max(norm(e0), eps);

fprintf('Eccentricity vector diagnostics:\n');
fprintf('  e0 = [% .8f  % .8f  % .8f]\n', e0(1), e0(2), e0(3));
fprintf('  |e0| = %.12f\n', norm(e0));
fprintf('  max relative drift in e-vector = %.6e\n', max(rel_e));

if doPlot
    figure('Color','w','Name','2BP -- eccentricity vector');
    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(T/3600, evec_hist(:,1), 'r-', 'LineWidth', 1.5); hold on;
    plot(T/3600, evec_hist(:,2), 'g-', 'LineWidth', 1.5);
    plot(T/3600, evec_hist(:,3), 'b-', 'LineWidth', 1.5);
    grid on; box on;
    ylabel('e_i');
    title('Eccentricity vector components');
    legend({'e_x', 'e_y', 'e_z'}, 'Location','best');

    nexttile;
    semilogy(T/3600, max(rel_e, eps), 'k-', 'LineWidth', 1.6);
    grid on; box on;
    xlabel('Time [hr]');
    ylabel('relative drift');
    title('Relative drift in the eccentricity vector');
end

hkPause(HK);

%% Section 12 -- Cartesian state <-> orbital elements
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 12 -- Cartesian state and orbital-element conversion');

%{
  IDEA
  ====
  The Kepler orbit may be described either by instantaneous Cartesian state
  vectors (r,v) or by six orbital elements.  This section verifies the round-trip

      (a,e,i,Omega,omega,M) -> (r,v) -> (a,e,i,Omega,omega,M)

  and also reconstructs the Cartesian state from the recovered elements.
%}

kep0 = rv2koe_local(r_rel0, v_rel0, mu_exact, 'deg');
[r_rt, v_rt] = koe2rv_local([kep0.a; kep0.e; kep0.i_deg; kep0.Omega_deg; kep0.omega_deg; kep0.M_deg], mu_exact, 'deg');

err_r_rt = norm(r_rt - r_rel0);
err_v_rt = norm(v_rt - v_rel0);

fprintf('Recovered orbital elements from initial Cartesian state:\n');
fprintf('  a     = %.9f km\n',  kep0.a);
fprintf('  e     = %.12f\n',    kep0.e);
fprintf('  i     = %.9f deg\n', kep0.i_deg);
fprintf('  Omega = %.9f deg\n', kep0.Omega_deg);
fprintf('  omega = %.9f deg\n', kep0.omega_deg);
fprintf('  M     = %.9f deg\n', kep0.M_deg);

fprintf('\nRound-trip state reconstruction errors:\n');
fprintf('  |dr| = %.6e km\n',   err_r_rt);
fprintf('  |dv| = %.6e km/s\n', err_v_rt);

hkPause(HK);

%% Section 13 -- Equatorial, ecliptic, and orbital-plane transformations
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 13 -- Coordinate transformations: equatorial, ecliptic, and orbital frames');

%{
  FRAME TRANSFORMATIONS
  =====================
  The same orbit can be represented in different inertial coordinates.

  We use:
    N_eq  : equatorial inertial frame
    N_ecl : ecliptic frame, obtained by a fixed obliquity rotation about x
    P     : perifocal/orbital frame, whose axes are
              p_hat = e_hat,
              q_hat = h_hat x e_hat,
              h_hat = orbit pole.

  The transformation from perifocal to equatorial inertial coordinates is the
  standard 3-1-3 Euler rotation built from (Omega, i, omega).
%}

R_ecl_to_eq = rot1_local(-C.eps_eq2ecl);
R_eq_to_ecl = R_ecl_to_eq.';

r_eq_hist  = r_rel_exact_hist;
r_ecl_hist = (R_eq_to_ecl * r_eq_hist.').';
r_p_hist   = (CASE.Q_P2N.' * r_eq_hist.').';

xhat_eq = [1;0;0]; yhat_eq = [0;1;0]; zhat_eq = [0;0;1];
xhat_ecl = R_eq_to_ecl * xhat_eq;
yhat_ecl = R_eq_to_ecl * yhat_eq;
zhat_ecl = R_eq_to_ecl * zhat_eq;

fprintf('Transformation diagnostics:\n');
fprintf('  Obliquity rotation = %.6f deg\n', C.eps_eq2ecl_deg);
fprintf('  z_eq expressed in ecliptic coordinates = [% .6f  % .6f  % .6f]\n', ...
    (R_eq_to_ecl * zhat_eq).');
fprintf('  h_hat in equatorial coordinates = [% .6f  % .6f  % .6f]\n', h0_hat(1), h0_hat(2), h0_hat(3));
fprintf('  h_hat in ecliptic coordinates   = [% .6f  % .6f  % .6f]\n', (R_eq_to_ecl*h0_hat).');

if doPlot
    figure('Color','w','Name','2BP -- frame transformations', 'Position', [60 60 1320 460]);
    tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile;
    hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on'); axis(ax1,'equal');
    plot3(ax1, r_eq_hist(:,1), r_eq_hist(:,2), r_eq_hist(:,3), 'r-', 'LineWidth', 1.8);
    quiver3(ax1,0,0,0, CASE.a*0.25,0,0, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.6);
    quiver3(ax1,0,0,0, 0,CASE.a*0.25,0, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.6);
    quiver3(ax1,0,0,0, 0,0,CASE.a*0.25, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.6);
    title(ax1,'Equatorial inertial coordinates');
    xlabel(ax1,'x_{eq} [km]'); ylabel(ax1,'y_{eq} [km]'); zlabel(ax1,'z_{eq} [km]');

    ax2 = nexttile;
    hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on'); axis(ax2,'equal');
    plot3(ax2, r_ecl_hist(:,1), r_ecl_hist(:,2), r_ecl_hist(:,3), 'r-', 'LineWidth', 1.8);
    quiver3(ax2,0,0,0, CASE.a*0.25,0,0, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.6);
    quiver3(ax2,0,0,0, 0,CASE.a*0.25,0, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.6);
    quiver3(ax2,0,0,0, 0,0,CASE.a*0.25, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.6);
    title(ax2,'Ecliptic inertial coordinates');
    xlabel(ax2,'x_{ecl} [km]'); ylabel(ax2,'y_{ecl} [km]'); zlabel(ax2,'z_{ecl} [km]');

    ax3 = nexttile;
    hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on'); axis(ax3,'equal');
    plot3(ax3, r_p_hist(:,1), r_p_hist(:,2), r_p_hist(:,3), 'r-', 'LineWidth', 1.8);
    quiver3(ax3,0,0,0, CASE.a*0.25,0,0, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.6);
    quiver3(ax3,0,0,0, 0,CASE.a*0.25,0, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.6);
    quiver3(ax3,0,0,0, 0,0,CASE.a*0.25, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.6);
    title(ax3,'Perifocal / orbital coordinates');
    xlabel(ax3,'p [km]'); ylabel(ax3,'q [km]'); zlabel(ax3,'h [km]');
end

disp(' ');
disp('Wrap-up takeaways:');
disp('  * The full inertial 2BP and the exact relative 1BP are mathematically equivalent.');
disp('  * In the body-centric limit m2 << m1, the gravitational parameter reduces to mu ~= G m1.');
disp('  * The Kepler problem preserves specific angular momentum, specific energy, and eccentricity vector.');
disp('  * The orbit is best understood geometrically in the orbital frame but expressed operationally in an inertial frame.');
disp('  * Mean anomaly evolves uniformly in time; eccentric and true anomaly do not.');

hkPause(HK);

%% ========================================================================
%  Local functions
%  ========================================================================

function CASE = makeDemoCase(name, C)
CASE = struct();
CASE.name = name;
CASE.nPeriods = 2.0;
CASE.nTimes   = 2000;
CASE.m2       = 1200;  % default spacecraft mass [kg]

switch lower(name)
    case 'circumterrestrial_cislunar'
        CASE.description   = 'Earth-centered eccentric orbit whose geocentric radial scale lies well inside the Moon''s orbit.';
        CASE.primary_name  = 'Earth';
        CASE.secondary_name= 'Spacecraft';
        CASE.frame_name    = 'ECI-like Earth-centered inertial frame';
        CASE.m1            = C.M_E;
        CASE.primary_radius= C.R_E;
        CASE.a             = 1.60e5;             % km
        CASE.e             = 0.35;
        CASE.i_deg         = 28.0;
        CASE.Omega_deg     = 35.0;
        CASE.omega_deg     = 40.0;
        CASE.M_deg         = 25.0;

    case 'circumterrestrial_translunar'
        CASE.description   = 'Earth-centered high-altitude elliptical orbit extending beyond the mean lunar distance in the idealized 2BP sense.';
        CASE.primary_name  = 'Earth';
        CASE.secondary_name= 'Spacecraft';
        CASE.frame_name    = 'ECI-like Earth-centered inertial frame';
        CASE.m1            = C.M_E;
        CASE.primary_radius= C.R_E;
        CASE.a             = 5.20e5;             % km
        CASE.e             = 0.24;
        CASE.i_deg         = 18.0;
        CASE.Omega_deg     = 20.0;
        CASE.omega_deg     = 55.0;
        CASE.M_deg         = 12.0;

    case 'circumlunar_selenocentric'
        CASE.description   = 'Moon-centered eccentric orbit in a nonrotating selenocentric inertial frame.';
        CASE.primary_name  = 'Moon';
        CASE.secondary_name= 'Orbiter';
        CASE.frame_name    = 'SCI-like selenocentric inertial frame';
        CASE.m1            = C.M_M;
        CASE.primary_radius= C.R_M;
        CASE.a             = 8.50e3;             % km
        CASE.e             = 0.18;
        CASE.i_deg         = 57.0;
        CASE.Omega_deg     = 42.0;
        CASE.omega_deg     = 33.0;
        CASE.M_deg         = 18.0;
        CASE.m2            = 600;                % kg

    otherwise
        error('Unknown CASE.name: %s', name);
end

CASE.mu_exact      = C.G * (CASE.m1 + CASE.m2);
CASE.mu_restricted = C.G * CASE.m1;
CASE.koe_deg       = [CASE.a; CASE.e; CASE.i_deg; CASE.Omega_deg; CASE.omega_deg; CASE.M_deg];
CASE.rp            = CASE.a * (1 - CASE.e);
CASE.ra            = CASE.a * (1 + CASE.e);
CASE.Q_P2N         = rot3_local(deg2rad(CASE.Omega_deg)) * rot1_local(deg2rad(CASE.i_deg)) * rot3_local(deg2rad(CASE.omega_deg));
end

function dxdt = eom_2bp_full(~, x, m1, m2, G)
R1 = x(1:3);
R2 = x(4:6);
V1 = x(7:9);
V2 = x(10:12);

r12 = R2 - R1;
r3  = norm(r12)^3;

A1 = +G*m2 * r12 / r3;
A2 = -G*m1 * r12 / r3;

dxdt = [V1; V2; A1; A2];
end

function dxdt = eom_2bp_relative(~, x, mu)
r = x(1:3);
v = x(4:6);
r3 = norm(r)^3;
a = -mu * r / r3;
dxdt = [v; a];
end

function [r, v, aux] = koe2rv_local(koe, mu, angletype)
% koe = [a; e; i; Omega; omega; M]
a     = koe(1);
e     = koe(2);
i     = koe(3);
Omega = koe(4);
omega = koe(5);
M     = koe(6);

if strcmpi(angletype,'deg')
    i     = deg2rad(i);
    Omega = deg2rad(Omega);
    omega = deg2rad(omega);
    M     = deg2rad(M);
end

if e >= 1
    error('This demo helper currently assumes an elliptic orbit with e < 1.');
end

M = wrapTo2Pi(M);
E = solveKeplerEllipse(M, e);
f = 2*atan2(sqrt(1 + e)*sin(E/2), sqrt(1 - e)*cos(E/2));
f = wrapTo2Pi(f);

p    = a*(1 - e^2);
r_pf = [p*cos(f)/(1 + e*cos(f));
        p*sin(f)/(1 + e*cos(f));
        0];

v_pf = sqrt(mu/p) * [-sin(f);
                      e + cos(f);
                      0];

Q = rot3_local(Omega) * rot1_local(i) * rot3_local(omega);
r = Q * r_pf;
v = Q * v_pf;

aux = struct();
aux.E = E;
aux.f = f;
aux.M = M;
aux.Q = Q;
end

function kep = rv2koe_local(r, v, mu, angletype)
tol  = 1e-12;
Khat = [0;0;1];

rmag = norm(r);
vmag = norm(v);
h    = cross(r, v);
hmag = norm(h);
hhat = h / hmag;

n    = cross(Khat, h);
nmag = norm(n);
if nmag > tol
    nhat = n / nmag;
else
    nhat = [1;0;0];
end

evec = cross(v, h) / mu - r / rmag;
e    = norm(evec);
if e > tol
    ehat = evec / e;
else
    ehat = nhat;
end

energy = 0.5*vmag^2 - mu/rmag;
if abs(e - 1) > 1e-10
    a = -mu / (2*energy);
else
    a = inf;
end

i = acosSafe(hhat(3));

if nmag > tol
    Omega = wrapTo2Pi(atan2(nhat(2), nhat(1)));
else
    Omega = 0;
end

if e > tol && nmag > tol
    omega = wrapTo2Pi(atan2(dot(cross(nhat, ehat), hhat), dot(nhat, ehat)));
elseif e > tol
    omega = wrapTo2Pi(atan2(ehat(2), ehat(1)));
else
    omega = 0;
end

if e > tol
    nu = wrapTo2Pi(atan2(dot(cross(ehat, r/rmag), hhat), dot(ehat, r/rmag)));
elseif nmag > tol
    nu = wrapTo2Pi(atan2(dot(cross(nhat, r/rmag), hhat), dot(nhat, r/rmag)));
else
    nu = wrapTo2Pi(atan2(r(2), r(1)));
end

if e < 1 - 1e-10
    E = 2*atan2(sqrt(1 - e)*sin(nu/2), sqrt(1 + e)*cos(nu/2));
    E = wrapTo2Pi(E);
    M = wrapTo2Pi(E - e*sin(E));
else
    E = NaN;
    M = NaN;
end

kep = struct();
kep.a         = a;
kep.e         = e;
kep.i         = i;
kep.Omega     = Omega;
kep.omega     = omega;
kep.nu        = nu;
kep.E         = E;
kep.M         = M;
kep.h         = h;
kep.hmag      = hmag;
kep.evec      = evec;
kep.energy    = energy;
kep.i_deg     = rad2deg(i);
kep.Omega_deg = rad2deg(Omega);
kep.omega_deg = rad2deg(omega);
kep.nu_deg    = rad2deg(nu);
kep.E_deg     = rad2deg(E);
kep.M_deg     = rad2deg(M);

if strcmpi(angletype,'deg')
    return
elseif strcmpi(angletype,'rad')
    return
else
    error('angletype must be ''deg'' or ''rad''.');
end
end

function E = solveKeplerEllipse(M, e)
if M < pi
    E = M + 0.5*e;
else
    E = M - 0.5*e;
end

for k = 1:50
    f  = E - e*sin(E) - M;
    fp = 1 - e*cos(E);
    dE = -f / fp;
    E  = E + dE;
    if abs(dE) < 1e-13
        break
    end
end
E = wrapTo2Pi(E);
end

function x = wrapTo2Pi(x)
x = mod(x, 2*pi);
if any(x < 0)
    x(x < 0) = x(x < 0) + 2*pi;
end
end

function v = rownorm(A)
v = sqrt(sum(A.^2, 2));
end

function c = acosSafe(x)
c = acos(max(-1, min(1, x)));
end

function R = rot1_local(a)
ca = cos(a); sa = sin(a);
R = [1  0   0;
     0 ca -sa;
     0 sa  ca];
end

function R = rot3_local(a)
ca = cos(a); sa = sin(a);
R = [ ca -sa  0;
      sa  ca  0;
       0   0  1];
end

function drawBody(ax, radius)
[xs, ys, zs] = sphere(40);
surf(ax, radius*xs, radius*ys, radius*zs, ...
    'FaceColor', [0.35 0.55 0.85], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.25, ...
    'HandleVisibility', 'off');
end

function HK = defaultHK()
HK.enabled         = true;
HK.doClc           = false;
HK.doClear         = false;
HK.doClose         = false;
HK.pauseAtSections = false;
HK.useKeyboard     = false;
HK.clcAfterSection = false;
HK.clearVarsKeep   = {'HK','fmtTitle','doPlot'};
end

function hkSection(HK)
if HK.doClose, close all; end
if HK.doClear
    keep = HK.clearVarsKeep;
    clearvars('-except', keep{:});
end
if HK.doClc, clc; end
end

function hkPause(HK)
if ~(isfield(HK,'pauseAtSections') && HK.pauseAtSections), return; end
fprintf('\n');
if isfield(HK,'useKeyboard') && HK.useKeyboard
    disp('--- INSPECT at K>> (dbcont to continue; dbquit to stop) ---');
    try
        w = evalin('caller','whos');
        for kk = 1:numel(w)
            name = w(kk).name;
            if isvarname(name)
                val = evalin('caller', name); %#ok<NASGU>
                eval([name ' = val;']); %#ok<EVLDIR>
            end
        end
        clear w kk name val
    catch
    end
    keyboard
else
    disp('--- press any key to continue ---');
    pause
end
if isfield(HK,'clcAfterSection') && HK.clcAfterSection, clc; end
end

function fmtTitle = defaultFmtTitle()
fmtTitle = @(s) localFmtTitle(s);
end

function localFmtTitle(s)
line = repmat('-', 1, 78);
fprintf('\n%s\n%s\n%s\n', line, s, line);
end
