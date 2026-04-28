%{
FILENAME:  mae240_bodycentric_NBP_demo.m

COURSE:    MAE 240 -- Space Flight Mechanics (Spring 2026)
REFERENCE: mae240_s26_newtonian_annotated.pdf  (lecture notes only)
BASELINE:  simulate_bodycentric.m

DESCRIPTION:
    A sectioned (%%) standalone MATLAB demo script for in-class lecture on
    the body-centric Newtonian N-body problem and the restricted (N+1)^st-
    body formulation for a spacecraft of small mass.

    The demo considers a Sun-centered, body-centric model of the inner Solar
    System with five massive bodies (Sun, Venus, Earth, Mars, and Jupiter),
    together with a small spacecraft initialized on a heliocentric,
    Hohmann-like transfer from Earth toward Jupiter. The script is organized
    for step-by-step discussion, with each section introducing a specific
    dynamical concept, model component, or numerical result.

    Lecture arc:
      (i)   body-centric reduction of the N-body equations
      (ii)  Solar-System-like initial conditions in a Sun-centered frame
      (iii) numerical integration of the body-centric massive-body dynamics
      (iv)  center-of-mass / barycenter interpretation in a non-inertial frame
      (v)   restricted (N+1)^st-body motion for a spacecraft of small mass
      (vi)  trajectory visualization and transfer diagnostics

INSTRUCTOR NOTES:
    - Run section-by-section (Ctrl+Enter) in lecture, or run the whole file.
    - All figures are guarded by doPlot.
    - The script is standalone: all helper functions are included below.
    - The lecture notes are the only formal reference for MAE 240.
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
fmtTitle('BODY-CENTRIC NBP DEMO -- heliocentric massive-body dynamics and restricted spacecraft motion');

disp('Lecture arc:');
disp('  Sec. 1  Body-centric reduction from the MAE 240 notes');
disp('  Sec. 2  Constants, bodies, and unit system');
disp('  Sec. 3  Solar-System-like initial conditions');
disp('  Sec. 4  Integrate the body-centric N-body equations');
disp('  Sec. 5  Massive-body trajectories and barycenter offset');
disp('  Sec. 6  Restricted spacecraft setup in the body-centric frame');
disp('  Sec. 7  Integrate the restricted spacecraft equations');
disp('  Sec. 8  Full trajectory visualization');
disp('  Sec. 9  Transfer diagnostics and lecture takeaways');

hkPause(HK);

%% Section 1 -- Body-centric reduction of the N-body equations
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 1 -- Body-centric form: direct and indirect terms');

%{
  LECTURE CONNECTION
  ==================
  In the MAE 240 notes, after introducing the barycentric reduction, the
  body-centric variables are defined by making all motion relative to a chosen
  body m_N.  For the remaining massive bodies i = 1,...,N-1, the notes write

      r_i_ddot = -G (m_N + m_i) r_i / |r_i|^3
                 - sum_{j=1, j~=i}^{N-1} G m_j
                   [ (r_i-r_j)/|r_i-r_j|^3 + r_j/|r_j|^3 ]

  where the first contribution is the DIRECT term and the bracketed sum is
  the INDIRECT term.  The notes emphasize that this frame is no longer
  inertial.

  For a spacecraft of negligible mass relative to the chosen central body,
  the restricted (N+1)^st-body equation becomes

      r_sc_ddot = -G m_N r_sc / |r_sc|^3
                  - sum_{j=1}^{N-1} G m_j
                    [ (r_sc-r_j)/|r_sc-r_j|^3 + r_j/|r_j|^3 ].

  That is the exact structure implemented later in this demo.
%}

disp('Key body-centric equations used in this demo:');
disp('  Massive bodies:');
disp('    r_i_ddot = -G (m_N+m_i) r_i / |r_i|^3 - sum_{j~=i} G m_j [ (r_i-r_j)/|r_i-r_j|^3 + r_j/|r_j|^3 ]');
disp('  Restricted spacecraft:');
disp('    r_sc_ddot = -G m_N r_sc / |r_sc|^3 - sum_j G m_j [ (r_sc-r_j)/|r_sc-r_j|^3 + r_j/|r_j|^3 ]');
disp('  Important point: the body-centric frame is not inertial.');

hkPause(HK);

%% Section 2 -- Constants, bodies, and unit system
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 2 -- Constants, body list, and km-kg-s units');

C.G       = 6.67430e-20;                % [km^3 / (kg s^2)]
C.M_sun   = 1.989e30;                   % [kg]
C.AU      = 1.496e8;                    % [km]
C.yr      = 365.25 * 24 * 3600;         % [s]
C.mu_sun  = C.G * C.M_sun;              % [km^3 / s^2]

BODY.names = {'Sun','Venus','Earth','Mars','Jupiter'};
BODY.n     = 5;
BODY.r_AU  = [0.00; 0.72; 1.00; 1.52; 5.20];
BODY.M_nom = [C.M_sun; 4.867e24; 5.972e24; 6.417e23; 1.898e27];

% Keep the simulate_bodycentric option, but default to deterministic masses.
USE.useMassVariation      = false;
USE.massVariationFactor   = 0.05;
USE.randomSeed            = 240;
USE.useRandomPhases       = true;

rng(USE.randomSeed);

if USE.useMassVariation
    scale = 1 + (2*rand(BODY.n-1,1) - 1) * USE.massVariationFactor;
    BODY.m = [BODY.M_nom(1); BODY.M_nom(2:end) .* scale];
else
    BODY.m = BODY.M_nom;
end

fprintf('Unit system: km, kg, s\n');
fprintf('G       = %.8e km^3/(kg s^2)\n', C.G);
fprintf('M_sun   = %.8e kg\n', C.M_sun);
fprintf('AU      = %.8e km\n', C.AU);
fprintf('mu_sun  = %.8e km^3/s^2\n', C.mu_sun);

fprintf('\nBody table:\n');
fprintf('  %-8s  %-12s  %-12s\n', 'Body', 'Mass [kg]', 'Radius [AU]');
fprintf('  %s\n', repmat('-',1,40));
for i = 1:BODY.n
    fprintf('  %-8s  %12.4e  %12.3f\n', BODY.names{i}, BODY.m(i), BODY.r_AU(i));
end

hkPause(HK);

%% Section 3 -- Solar-System-like initial conditions in the Sun-centered frame
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 3 -- Solar-System-like body-centric initial conditions');

%{
  BASE SETUP
  ==========
  This section follows simulate_bodycentric.m:
    * the frame origin is fixed at the Sun,
    * the planets start on approximate circular, coplanar heliocentric orbits,
    * the Sun remains at the origin of the chosen body-centric frame,
    * the initial-value problem is then integrated with the body-centric EOM.

  Because the setup is Sun-centered rather than barycentric, the origin does
  not coincide with the true system barycenter.
%}

n = BODY.n;
r0_all = zeros(3*n,1);
v0_all = zeros(3*n,1);

theta0 = zeros(n,1);
for i = 2:n
    rmag = BODY.r_AU(i) * C.AU;
    if USE.useRandomPhases
        theta0(i) = 2*pi*rand;
    else
        theta0(i) = 2*pi*(i-2)/(n-1);
    end

    pos = [rmag*cos(theta0(i)); rmag*sin(theta0(i)); 0];
    vc  = sqrt(C.mu_sun / rmag);
    vel = [-vc*sin(theta0(i)); vc*cos(theta0(i)); 0];

    r0_all(3*i-2:3*i) = pos;
    v0_all(3*i-2:3*i) = vel;
end

% Sun fixed at the body-center origin
r0_all(1:3) = [0; 0; 0];
v0_all(1:3) = [0; 0; 0];

x0 = [r0_all; v0_all];

fprintf('Initial conditions summary:\n');
fprintf('  %-8s  %-11s  %-28s  %-28s\n', 'Body', 'theta [deg]', 'r_0 [km]', 'v_0 [km/s]');
fprintf('  %s\n', repmat('-',1,88));
for i = 1:n
    ri = r0_all(3*i-2:3*i).';
    vi = v0_all(3*i-2:3*i).';
    fprintf('  %-8s  %11.3f  [%+.3e %+.3e %+.3e]  [%+.3e %+.3e %+.3e]\n', ...
        BODY.names{i}, rad2deg(theta0(i)), ri(1), ri(2), ri(3), vi(1), vi(2), vi(3));
end

hkPause(HK);

%% Section 4 -- Integrate the body-centric N-body equations for the massive bodies
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 4 -- Integrating the body-centric N-body equations');

SIM.tEndYears = 12;
SIM.nTimes    = 5000;
SIM.tEnd      = SIM.tEndYears * C.yr;
SIM.tspan     = linspace(0, SIM.tEnd, SIM.nTimes);
SIM.options   = odeset('AbsTol',1e-13,'RelTol',1e-11);

fprintf('Integrating %d massive bodies for %.1f years (%d output times)...\n', ...
    n, SIM.tEndYears, SIM.nTimes);

tic
[T, X] = ode45(@(t,x) eom_NBP_bodycentric(t, x, BODY.m, C.G), SIM.tspan, x0, SIM.options);
t_body = toc;

fprintf('Massive-body integration complete in %.3f s\n', t_body);
fprintf('State history size: %d x %d\n', size(X,1), size(X,2));

r_all = X(:, 1:3*n);
v_all = X(:, 3*n+1:end);

hkPause(HK);

%% Section 5 -- Massive-body trajectories and barycenter interpretation
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 5 -- Massive-body trajectories and barycenter offset from the Sun');

%{
  LECTURE CONNECTION
  ==================
  The notes distinguish the barycentric frame from the body-centric frame.
  Here the origin is fixed at the Sun, so the center of mass generally does
  not sit at the origin.  This is a useful visual reminder that the chosen
  frame is not inertial.
%}

r_cm = centerOfMassHistory(r_all, BODY.m);
r_cm_norm = vecnorm(r_cm, 2, 2);

fprintf('Maximum distance of the system barycenter from the Sun-centered origin: %.4e km (%.4e AU)\n', ...
    max(r_cm_norm), max(r_cm_norm)/C.AU);

if doPlot
    colors = [0.95 0.85 0.10;
              0.60 0.80 0.95;
              0.20 0.70 0.30;
              0.85 0.40 0.20;
              0.70 0.50 0.85];

    figure('Color','w','Name','MAE 240 -- Body-centric massive-body trajectories');
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile; hold on; grid on; box on; axis equal;
    h = gobjects(n+1,1);
    for i = 1:n
        ri = r_all(:,3*i-2:3*i);
        h(i) = plot3(ri(:,1)/C.AU, ri(:,2)/C.AU, ri(:,3)/C.AU, '-', ...
            'Color', colors(i,:), 'LineWidth', 1.6, 'DisplayName', BODY.names{i});
        plot3(ri(1,1)/C.AU, ri(1,2)/C.AU, ri(1,3)/C.AU, 'o', ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'HandleVisibility','off');
    end
    h(end) = plot3(r_cm(:,1)/C.AU, r_cm(:,2)/C.AU, r_cm(:,3)/C.AU, 'k--', ...
        'LineWidth', 1.6, 'DisplayName', 'Center of Mass');
    xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
    title('Massive-body trajectories in the Sun-centered frame');
    legend(h,'Location','best');
    view(3);

    nexttile;
    plot(T/C.yr, r_cm(:,1)/C.AU, 'r-', 'LineWidth', 1.3); hold on;
    plot(T/C.yr, r_cm(:,2)/C.AU, 'b-', 'LineWidth', 1.3);
    plot(T/C.yr, r_cm(:,3)/C.AU, 'k-', 'LineWidth', 1.3);
    grid on; box on;
    xlabel('Time [yr]');
    ylabel('r_c [AU]');
    legend('x','y','z','Location','best');
    title('Barycenter location relative to the Sun-centered origin');

    title(tl,'Body-centric massive-body dynamics','FontSize',13);
end

hkPause(HK);

%% Section 6 -- Restricted spacecraft setup: Earth departure and heliocentric transfer guess
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 6 -- Restricted spacecraft setup in the body-centric frame');

%{
  LECTURE CONNECTION
  ==================
  The notes next introduce the restricted (N+1)^st-body problem for a
  particle of small mass.  In this demo, that particle is a spacecraft.
  We use the massive-body ephemerides from the previous section and then
  propagate the spacecraft in the heliocentric body-centric frame.

  For a classroom-friendly transfer example, we design the nominal departure
  state so that the Sun-only Hohmann ellipse has apoapsis exactly at
  Jupiter's orbital radius.  The spacecraft is then propagated over the
  full multi-revolution time window, but the speed tuning is performed on
  the *first outbound apoapsis* near the Hohmann half-period.  In that way,
  the complete arc remains visible while the initial transfer still arrives
  at Jupiter's orbital radius without overshooting it at first apoapsis.
%}

SC.iEarth          = 3;
SC.iJupiter        = 5;
SC.offsetAU        = 1e-3;
SC.maxTuneIter     = 20;
SC.speedScaleLo    = 0.98;
SC.speedScaleHi    = 1.02;
SC.apoWindowLoFrac = 0.35;
SC.apoWindowHiFrac = 1.75;

r0_earth = r0_all(3*SC.iEarth-2:3*SC.iEarth);
v0_earth = v0_all(3*SC.iEarth-2:3*SC.iEarth);
r0_jup   = r0_all(3*SC.iJupiter-2:3*SC.iJupiter);

r1 = norm(r0_earth);
r2 = norm(r0_jup);

u_r_earth = r0_earth / norm(r0_earth);
u_t_earth = v0_earth / norm(v0_earth);

r0_sc   = r0_earth + SC.offsetAU * C.AU * u_r_earth;
r1_sc   = norm(r0_sc);
r2_goal = r2;
a_trans = 0.5 * (r1_sc + r2_goal);
T_trans = pi * sqrt(a_trans^3 / C.mu_sun);
v_dep   = sqrt(C.mu_sun * (2/r1_sc - 1/a_trans));
SC.speedScale = 1.0;

SC.tEnd = T(end);
v0_sc = SC.speedScale * v_dep * u_t_earth;
x0_sc = [r0_sc; v0_sc];

r_interp = makeInterpolants(T, r_all, n);

[a_dir0, a_pert0, a_by_body0] = spacecraftAccelerationBreakdown(0, x0_sc, r_interp, BODY.m, C.G);

fprintf('Earth departure radius             : %.4e km (%.4f AU)\n', r1, r1/C.AU);
fprintf('Spacecraft initial radius          : %.4e km (%.4f AU)\n', r1_sc, r1_sc/C.AU);
fprintf('Jupiter orbital radius             : %.4e km (%.4f AU)\n', r2, r2/C.AU);
fprintf('Design apoapsis radius             : %.4e km (%.4f AU)\n', r2_goal, r2_goal/C.AU);
fprintf('Hohmann transfer a                 : %.4e km (%.4f AU)\n', a_trans, a_trans/C.AU);
fprintf('Hohmann half-period estimate       : %.4f years\n', T_trans/C.yr);
fprintf('Spacecraft propagation window      : %.4f years (multi-revolution arc)\n', SC.tEnd/C.yr);
fprintf('Nominal departure speed            : %.4f km/s\n', v_dep);
fprintf('Direct solar acceleration at departure    : %.4e km/s^2\n', norm(a_dir0));
fprintf('Total disturbing acceleration at departure: %.4e km/s^2\n', norm(a_pert0));

fprintf('\nPlanet-by-planet disturbing acceleration norms at departure:\n');
for j = 2:n
    fprintf('  %-8s  %.4e km/s^2\n', BODY.names{j}, norm(a_by_body0(:,j)));
end

hkPause(HK);

%% Section 7 -- Integrate the restricted spacecraft equations
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 7 -- Integrating the restricted spacecraft equations');

fprintf('Integrating restricted spacecraft over the full %.4f-year multi-revolution window...\n', SC.tEnd/C.yr);

% Tune the spacecraft speed so the *first outbound apoapsis* approaches
% Jupiter's orbital radius from below, without exceeding it.
scale_lo = SC.speedScaleLo;
scale_hi = SC.speedScaleHi;
bestScale = scale_lo;
bestGap   = inf;
tuneIter  = 0;
tic

% Build an actual bracket on first-apoapsis radius.  If needed, expand the
% upper bound until the first outbound apoapsis reaches beyond Jupiter's orbit.
SC.speedScale = scale_lo;
v0_sc = SC.speedScale * v_dep * u_t_earth;
x0_sc = [r0_sc; v0_sc];
[T_sc_tmp, X_sc_tmp] = ode45(@(t,x) eom_sc_bodycentric(t, x, r_interp, BODY.m, C.G), [T(1) SC.tEnd], x0_sc, SIM.options);
[firstApo_lo, ~, ~] = firstOutboundApoapsis(T_sc_tmp, X_sc_tmp(:,1:3), X_sc_tmp(:,4:6), T_trans, SC.apoWindowLoFrac, SC.apoWindowHiFrac);

SC.speedScale = scale_hi;
v0_sc = SC.speedScale * v_dep * u_t_earth;
x0_sc = [r0_sc; v0_sc];
[T_sc_tmp, X_sc_tmp] = ode45(@(t,x) eom_sc_bodycentric(t, x, r_interp, BODY.m, C.G), [T(1) SC.tEnd], x0_sc, SIM.options);
[firstApo_hi, ~, ~] = firstOutboundApoapsis(T_sc_tmp, X_sc_tmp(:,1:3), X_sc_tmp(:,4:6), T_trans, SC.apoWindowLoFrac, SC.apoWindowHiFrac);

expandCount = 0;
while firstApo_hi < r2 && expandCount < 8
    scale_lo = scale_hi;
    firstApo_lo = firstApo_hi;
    scale_hi = scale_hi + 0.01;

    SC.speedScale = scale_hi;
    v0_sc = SC.speedScale * v_dep * u_t_earth;
    x0_sc = [r0_sc; v0_sc];
    [T_sc_tmp, X_sc_tmp] = ode45(@(t,x) eom_sc_bodycentric(t, x, r_interp, BODY.m, C.G), [T(1) SC.tEnd], x0_sc, SIM.options);
    [firstApo_hi, ~, ~] = firstOutboundApoapsis(T_sc_tmp, X_sc_tmp(:,1:3), X_sc_tmp(:,4:6), T_trans, SC.apoWindowLoFrac, SC.apoWindowHiFrac);
    expandCount = expandCount + 1;
end

if firstApo_lo > r2
    warning('Lower speed bracket already overshoots Jupiter orbit at first apoapsis.');
end
if firstApo_hi < r2
    warning('Unable to bracket Jupiter-orbit arrival at first apoapsis; using the closest below-target case found.');
end

for k = 1:SC.maxTuneIter
    tuneIter = k;
    SC.speedScale = 0.5 * (scale_lo + scale_hi);
    v0_sc = SC.speedScale * v_dep * u_t_earth;
    x0_sc = [r0_sc; v0_sc];

    [T_sc, X_sc] = ode45(@(t,x) eom_sc_bodycentric(t, x, r_interp, BODY.m, C.G), [T(1) SC.tEnd], x0_sc, SIM.options);

    r_sc_all  = X_sc(:,1:3);
    v_sc_all  = X_sc(:,4:6);
    [firstApo_trial, idxApo_trial, tApo_trial] = firstOutboundApoapsis(T_sc, r_sc_all, v_sc_all, T_trans, SC.apoWindowLoFrac, SC.apoWindowHiFrac);

    if firstApo_trial <= r2
        gap = r2 - firstApo_trial;
        if gap < bestGap
            bestGap     = gap;
            bestScale   = SC.speedScale;
            bestApo     = firstApo_trial;
            bestIdxApo  = idxApo_trial;
            bestTApo    = tApo_trial;
        end
        scale_lo = SC.speedScale;
    else
        scale_hi = SC.speedScale;
    end
end

SC.speedScale = bestScale;
v0_sc = SC.speedScale * v_dep * u_t_earth;
x0_sc = [r0_sc; v0_sc];
[T_sc, X_sc] = ode45(@(t,x) eom_sc_bodycentric(t, x, r_interp, BODY.m, C.G), [T(1) SC.tEnd], x0_sc, SIM.options);

r_sc_all  = X_sc(:,1:3);
v_sc_all  = X_sc(:,4:6);
r_sc_norm = vecnorm(r_sc_all, 2, 2);
[maxRsc, idxMaxR] = max(r_sc_norm);
[firstApoR, idxFirstApo, tFirstApo] = firstOutboundApoapsis(T_sc, r_sc_all, v_sc_all, T_trans, SC.apoWindowLoFrac, SC.apoWindowHiFrac);
if isinf(bestGap)
    warning('No below-target first-apoapsis case was found during tuning; using the lower bound speed scale.');
end

t_sc = toc;

fprintf('Spacecraft integration complete in %.3f s\n', t_sc);
fprintf('Spacecraft state history size: %d x %d\n', size(X_sc,1), size(X_sc,2));
fprintf('Departure speed scale applied: %.8f (after %d tuning iteration(s))\n', SC.speedScale, tuneIter);

% Interpolate Earth/Jupiter trajectories onto the spacecraft time grid.
r_earth_sc = reshape(r_interp{SC.iEarth}(T_sc), [], 3);
r_jup_sc   = reshape(r_interp{SC.iJupiter}(T_sc), [], 3);

d_sc_earth = vecnorm(r_sc_all - r_earth_sc, 2, 2);
d_sc_jup   = vecnorm(r_sc_all - r_jup_sc,   2, 2);
r_sc_norm  = vecnorm(r_sc_all, 2, 2);

[scEnergy, scAngMom] = restrictedInvariants(r_sc_all, v_sc_all, C.mu_sun);
rel_dE_sc = abs(scEnergy - scEnergy(1)) / max(abs(scEnergy(1)), eps);
rel_dH_sc = vecnorm(scAngMom - scAngMom(1,:), 2, 2) / max(norm(scAngMom(1,:)), eps);

fprintf('Max relative drift in Sun-only specific energy proxy: %.4e\n', max(rel_dE_sc));
fprintf('Max relative drift in Sun-only specific angular momentum proxy: %.4e\n', max(rel_dH_sc));

hkPause(HK);

%% Section 8 -- Full trajectory visualization: planets and spacecraft
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 8 -- Trajectory visualization in the heliocentric body-centric frame');

if doPlot
    colors = [0.95 0.85 0.10;
              0.60 0.80 0.95;
              0.20 0.70 0.30;
              0.85 0.40 0.20;
              0.70 0.50 0.85];

    figure('Color','w','Name','MAE 240 -- Body-centric spacecraft transfer');
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile([2 1]); hold on; grid on; box on; axis equal;
    h = gobjects(n+1,1);
    for i = 1:n
        ri = r_all(:,3*i-2:3*i);
        h(i) = plot3(ri(:,1)/C.AU, ri(:,2)/C.AU, ri(:,3)/C.AU, '-', ...
            'Color', colors(i,:), 'LineWidth', 1.4, 'DisplayName', BODY.names{i});
        plot3(ri(1,1)/C.AU, ri(1,2)/C.AU, ri(1,3)/C.AU, 'o', ...
            'Color', colors(i,:), 'LineWidth', 1.0, 'HandleVisibility','off');
    end
    h(end) = plot3(r_sc_all(:,1)/C.AU, r_sc_all(:,2)/C.AU, r_sc_all(:,3)/C.AU, 'k-', ...
        'LineWidth', 1.8, 'DisplayName', 'Spacecraft');
    plot3(r0_sc(1)/C.AU, r0_sc(2)/C.AU, r0_sc(3)/C.AU, 'ks', 'LineWidth', 1.2, 'HandleVisibility','off');
    xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
    title('Body-centric trajectories: massive bodies + spacecraft');
    legend(h, 'Location', 'best');
    view(3);

    nexttile;
    plot(T_sc/C.yr, d_sc_earth/C.AU, 'b-', 'LineWidth', 1.5); hold on;
    plot(T_sc/C.yr, d_sc_jup/C.AU,   'r-', 'LineWidth', 1.5);
    grid on; box on;
    xlabel('Time [yr]');
    ylabel('Distance [AU]');
    legend('SC-Earth','SC-Jupiter','Location','best');
    title('Spacecraft distance to departure and target planets');

    nexttile;
    plot(T_sc/C.yr, r_sc_norm/C.AU, 'k-', 'LineWidth', 1.5); hold on;
    yline(r2/C.AU, 'r--', 'Jupiter orbit', 'LineWidth', 1.1, 'LabelHorizontalAlignment','left');
    plot(tFirstApo/C.yr, firstApoR/C.AU, 'ko', 'MarkerFaceColor','k', 'HandleVisibility','off');
    grid on; box on;
    xlabel('Time [yr]');
    ylabel('|r_{sc}| [AU]');
    title('Spacecraft heliocentric radius');

    title(tl, 'Restricted (N+1)^st-body motion in the Sun-centered body-centric frame', 'FontSize', 13);
end

hkPause(HK);

%% Section 9 -- Transfer diagnostics and lecture takeaways
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 9 -- Transfer diagnostics and lecture takeaways');

[minEarth, idxEarth] = min(d_sc_earth);
[minJup,   idxJup]   = min(d_sc_jup);
[maxRsc,   idxMaxR]  = max(r_sc_norm);
[firstApoR, idxFirstApo, tFirstApo] = firstOutboundApoapsis(T_sc, r_sc_all, v_sc_all, T_trans, SC.apoWindowLoFrac, SC.apoWindowHiFrac);
radiusDeficit = r2 - firstApoR;

fprintf('Closest approach to Earth   : %.4e km (%.4e AU) at t = %.4f yr\n', ...
    minEarth, minEarth/C.AU, T_sc(idxEarth)/C.yr);
fprintf('Closest approach to Jupiter : %.4e km (%.4e AU) at t = %.4f yr\n', ...
    minJup, minJup/C.AU, T_sc(idxJup)/C.yr);
fprintf('First outbound apoapsis radius        : %.4e km (%.4f AU) at t = %.4f yr\n', ...
    firstApoR, firstApoR/C.AU, tFirstApo/C.yr);
fprintf('Maximum spacecraft heliocentric radius: %.4e km (%.4f AU) at t = %.4f yr\n', ...
    maxRsc, maxRsc/C.AU, T_sc(idxMaxR)/C.yr);
fprintf('Hohmann half-period estimate: %.4f yr\n', T_trans/C.yr);
fprintf('Time of minimum SC-Jupiter distance: %.4f yr\n', T_sc(idxJup)/C.yr);
fprintf('First-apoapsis radius deficit relative to Jupiter orbit: %.4e km (%.4e AU)\n', ...
    radiusDeficit, radiusDeficit/C.AU);

fprintf('\nDiscussion cues:\n');
fprintf('  1. The massive bodies are propagated in the Sun-centered body-centric frame.\n');
fprintf('  2. The barycenter therefore moves relative to the origin.\n');
fprintf('  3. The spacecraft sees a direct solar term plus indirect planetary terms.\n');
fprintf('  4. The nominal Sun-only transfer has apoapsis at Jupiter''s orbital radius.\n');
fprintf('  5. The full multi-revolution spacecraft arc is shown, but the speed is tuned so the first outbound apoapsis reaches Jupiter''s orbital radius without exceeding it.\n');

hkPause(HK);

% -------------------------------------------------------------------------
% END OF SCRIPT
% -------------------------------------------------------------------------

function dx = eom_NBP_bodycentric(~, x, m_bodies, G)
% Body-centric Newtonian N-body equations with body 1 chosen as the center.
% The first body is fixed at the origin of the chosen frame.

n = numel(m_bodies);
r = x(1:3*n);
v = x(3*n+1:end);

a = zeros(3*n,1);

% Body 1 remains fixed at the body-centric origin.
a(1:3) = [0;0;0];

for i = 2:n
    ri = r(3*i-2:3*i);
    ri_norm = norm(ri);

    % Direct term from the central body.
    a_i = -G * (m_bodies(1) + m_bodies(i)) * ri / max(ri_norm^3, eps);

    % Indirect terms from the other noncentral bodies.
    for j = 2:n
        if j == i
            continue
        end
        rj = r(3*j-2:3*j);
        rij = ri - rj;
        a_i = a_i - G * m_bodies(j) * ( rij / max(norm(rij)^3, eps) + rj / max(norm(rj)^3, eps) );
    end

    a(3*i-2:3*i) = a_i;
end

dx = [v; a];
end

function dx = eom_sc_bodycentric(t, x, r_interp, m_bodies, G)
% Restricted (N+1)^st-body equation in a Sun-centered body-centric frame.

n = numel(m_bodies);
r_sc = x(1:3);
v_sc = x(4:6);

a_sc = -G * m_bodies(1) * r_sc / max(norm(r_sc)^3, eps);

for j = 2:n
    rj = reshape(r_interp{j}(t), [3,1]);
    r_scj = r_sc - rj;
    a_sc = a_sc - G * m_bodies(j) * ( r_scj / max(norm(r_scj)^3, eps) + rj / max(norm(rj)^3, eps) );
end

dx = [v_sc; a_sc];
end

function r_cm = centerOfMassHistory(r_all, m_bodies)
Nt = size(r_all,1);
n  = numel(m_bodies);
r_cm = zeros(Nt,3);
M = sum(m_bodies);

for k = 1:Nt
    rk = reshape(r_all(k,:), 3, n);
    r_cm(k,:) = (rk * m_bodies / M).';
end
end

function r_interp = makeInterpolants(T, r_all, n)
r_interp = cell(n,1);
for i = 1:n
    ri = r_all(:,3*i-2:3*i);
    r_interp{i} = griddedInterpolant(T, ri, 'spline', 'nearest');
end
end

function [a_dir, a_pert, a_by_body] = spacecraftAccelerationBreakdown(t, x_sc, r_interp, m_bodies, G)
n = numel(m_bodies);
r_sc = x_sc(1:3);
a_dir = -G * m_bodies(1) * r_sc / max(norm(r_sc)^3, eps);
a_pert = zeros(3,1);
a_by_body = zeros(3,n);

for j = 2:n
    rj = reshape(r_interp{j}(t), [3,1]);
    r_scj = r_sc - rj;
    a_j = -G * m_bodies(j) * ( r_scj / max(norm(r_scj)^3, eps) + rj / max(norm(rj)^3, eps) );
    a_by_body(:,j) = a_j;
    a_pert = a_pert + a_j;
end
end

function [eps_hist, h_hist] = restrictedInvariants(r_hist, v_hist, mu)
Nt = size(r_hist,1);
eps_hist = zeros(Nt,1);
h_hist   = zeros(Nt,3);
for k = 1:Nt
    r = r_hist(k,:);
    v = v_hist(k,:);
    eps_hist(k) = 0.5 * dot(v,v) - mu / norm(r);
    h_hist(k,:) = cross(r, v);
end
end


function [rApo, idxApo, tApo] = firstOutboundApoapsis(T, r_hist, v_hist, T_guess, fracLo, fracHi)
r_norm = vecnorm(r_hist, 2, 2);
rdot   = sum(r_hist .* v_hist, 2) ./ max(r_norm, eps);

mask = (T >= fracLo * T_guess) & (T <= fracHi * T_guess);
idxWindow = find(mask);
if numel(idxWindow) < 3
    [rApo, idxApo] = max(r_norm);
    tApo = T(idxApo);
    return
end

idxApo = [];
for kk = idxWindow(1)+1 : idxWindow(end)-1
    if rdot(kk-1) >= 0 && rdot(kk+1) <= 0
        i1 = max(1, kk-2);
        i2 = min(numel(T), kk+2);
        [~, localIdx] = max(r_norm(i1:i2));
        idxApo = i1 + localIdx - 1;
        break
    end
end

if isempty(idxApo)
    [~, localIdx] = max(r_norm(idxWindow));
    idxApo = idxWindow(localIdx);
end

rApo = r_norm(idxApo);
tApo = T(idxApo);
end

function HK = defaultHK()
HK.enabled         = true;
HK.doClc           = false;
HK.doClear         = false;
HK.doClose         = false;
HK.pauseAtSections = true;
HK.useKeyboard     = true;
HK.clcAfterSection = true;
HK.clearVarsKeep   = {'HK','fmtTitle','doPlot'};
end

function hkSection(HK)
if HK.doClc
    clc
end
if HK.doClose
    close all
end
if HK.doClear
    evalin('base', sprintf('clearvars -except %s', strjoin(HK.clearVarsKeep, ' ')));
end
end

function hkPause(HK)
if ~HK.pauseAtSections
    return
end
if HK.useKeyboard
    disp(' ');
    disp('Paused at section boundary. Type dbcont to continue.');
    keyboard
else
    input('Press Enter to continue.','s');
end
if HK.clcAfterSection
    clc
end
end

function fmtTitle = defaultFmtTitle()
fmtTitle = @localTitle;
    function localTitle(txt)
        fprintf('\n%s\n', repmat('=',1,78));
        fprintf('%s\n', txt);
        fprintf('%s\n\n', repmat('=',1,78));
    end
end
