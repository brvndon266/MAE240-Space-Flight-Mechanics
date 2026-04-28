%% MAE 240 HW — Problems 1 and 2
% Problem 1: Barycentric N-body problem
% Problem 2: Bodycentric N-body problem + restricted spacecraft

clear; close all; clc;

%% Constants
G     = 6.67430e-20;          % km^3 kg^-1 s^-2
M_sun = 1.989e30;             % kg
AU    = 1.496e8;              % km
yr    = 365.25*24*3600;       % s

%% Bodies
n = 5;
body_names = {'Sun','Venus','Earth','Mars','Jupiter'};

m_bodies = [
    M_sun
    4.867e24
    5.972e24
    6.417e23
    1.898e27
];

r_AU = [0; 0.723; 1.000; 1.524; 5.203];

%% Initial conditions
r0_all = zeros(3*n,1);
v0_all = zeros(3*n,1);

theta0 = linspace(0,2*pi,n+1);
theta0 = theta0(1:n);

% Sun
r0_all(1:3) = [0;0;0];
v0_all(1:3) = [0;0;0];

% Planets
for i = 2:n
    r_mag = r_AU(i)*AU;
    theta = theta0(i);

    r0_all(3*i-2:3*i) = [
        r_mag*cos(theta)
        r_mag*sin(theta)
        0
    ];

    v_c = sqrt(G*M_sun/r_mag);

    v0_all(3*i-2:3*i) = [
        -v_c*sin(theta)
         v_c*cos(theta)
         0
    ];
end

x0_bary = [r0_all; v0_all];

%% Time span
T_end = 12*yr;
tspan = linspace(0,T_end,5000);
opts = odeset('RelTol',1e-11,'AbsTol',1e-13);

%% ============================================================
%  PROBLEM 1 — BARYCENTRIC N-BODY MODEL
% ============================================================

fprintf('\nRunning Problem 1: barycentric N-body model...\n');

[T_bary,X_bary] = ode45(@(t,x) eom_NBP_bary(t,x,m_bodies,G), ...
                        tspan,x0_bary,opts);

r_bary = X_bary(:,1:3*n);
v_bary = X_bary(:,3*n+1:end);

%% Plot barycentric trajectories
figure('Color','w','Name','Problem 1 — Barycentric N-body');
hold on; grid on; axis equal;

colors = lines(n);

for i = 1:n
    ri = r_bary(:,3*i-2:3*i);

    plot3(ri(:,1)/AU,ri(:,2)/AU,ri(:,3)/AU, ...
        'LineWidth',1.5, ...
        'Color',colors(i,:), ...
        'DisplayName',body_names{i});

    plot3(ri(1,1)/AU,ri(1,2)/AU,ri(1,3)/AU, ...
        'o','Color',colors(i,:), ...
        'HandleVisibility','off');

    plot3(ri(end,1)/AU,ri(end,2)/AU,ri(end,3)/AU, ...
        's','Color',colors(i,:), ...
        'HandleVisibility','off');
end

xlabel('x [AU]');
ylabel('y [AU]');
zlabel('z [AU]');
title('Problem 1: Barycentric N-body trajectories');
legend('Location','best');
view(3);

%% Barycenter check
r_cm_bary = centerOfMassHistory(r_bary,m_bodies);
P_bary = linearMomentumHistory(v_bary,m_bodies);
V_cm0 = P_bary(1,:)/sum(m_bodies);
r_cm_linear = r_cm_bary(1,:) + T_bary(:)*V_cm0;

figure('Color','w','Name','Problem 1 — Barycenter Check');
plot(T_bary/yr,r_cm_bary(:,1)/AU,'LineWidth',1.5); hold on;
plot(T_bary/yr,r_cm_linear(:,1)/AU,'--','LineWidth',1.5);
grid on;
xlabel('Time [yr]');
ylabel('R_{cm,x} [AU]');
legend('Numerical barycenter','Linear prediction');
title('Problem 1: Center of mass check');

%% Energy check
[KE_bary,PE_bary,E_bary] = totalEnergyHistory(r_bary,v_bary,m_bodies,G);
rel_E_err_bary = abs(E_bary - E_bary(1))/abs(E_bary(1));

figure('Color','w','Name','Problem 1 — Energy Error');
semilogy(T_bary/yr,rel_E_err_bary,'LineWidth',1.5);
grid on;
xlabel('Time [yr]');
ylabel('|E - E_0| / |E_0|');
title('Problem 1: Relative energy error');

fprintf('Problem 1 max relative energy error = %.4e\n',max(rel_E_err_bary));

%% ============================================================
%  PROBLEM 2A — BODYCENTRIC N-BODY MODEL
% ============================================================

fprintf('\nRunning Problem 2A: bodycentric N-body model...\n');

% Convert barycentric initial conditions to Sun-centered bodycentric coordinates
r0_body = zeros(3*n,1);
v0_body = zeros(3*n,1);

r_sun0 = r0_all(1:3);
v_sun0 = v0_all(1:3);

for i = 1:n
    r0_body(3*i-2:3*i) = r0_all(3*i-2:3*i) - r_sun0;
    v0_body(3*i-2:3*i) = v0_all(3*i-2:3*i) - v_sun0;
end

x0_body = [r0_body; v0_body];

[T_body,X_body] = ode45(@(t,x) eom_NBP_bodycentric(t,x,m_bodies,G), ...
                        tspan,x0_body,opts);

r_body = X_body(:,1:3*n);
v_body = X_body(:,3*n+1:end);

%% Convert barycentric solution to Sun-centered form for comparison
r_bary_rel = zeros(size(r_bary));

for k = 1:length(T_bary)
    r_sun = r_bary(k,1:3);

    for i = 1:n
        r_bary_rel(k,3*i-2:3*i) = r_bary(k,3*i-2:3*i) - r_sun;
    end
end

%% Plot bodycentric trajectories
figure('Color','w','Name','Problem 2A — Bodycentric N-body');
hold on; grid on; axis equal;

for i = 1:n
    ri = r_body(:,3*i-2:3*i);

    plot3(ri(:,1)/AU,ri(:,2)/AU,ri(:,3)/AU, ...
        'LineWidth',1.5, ...
        'Color',colors(i,:), ...
        'DisplayName',body_names{i});
end

xlabel('x [AU]');
ylabel('y [AU]');
zlabel('z [AU]');
title('Problem 2A: Bodycentric trajectories');
legend('Location','best');
view(3);

%% Compare Problem 1 and Problem 2A
figure('Color','w','Name','Problem 2A — Comparison Error');
hold on; grid on;

for i = 2:n
    err_i = vecnorm(r_body(:,3*i-2:3*i) - r_bary_rel(:,3*i-2:3*i),2,2);
    plot(T_body/yr,err_i,'LineWidth',1.5,'DisplayName',body_names{i});
end

xlabel('Time [yr]');
ylabel('Position difference [km]');
title('Problem 2A: Bodycentric vs barycentric-relative comparison');
legend('Location','best');

%% ============================================================
%  PROBLEM 2B — RESTRICTED SPACECRAFT
% ============================================================

fprintf('\nRunning Problem 2B: restricted spacecraft...\n');

% Use Earth departure and Jupiter transfer-like initial condition
iEarth = 3;
iJupiter = 5;

r0_earth = r_body(1,3*iEarth-2:3*iEarth).';
v0_earth = v_body(1,3*iEarth-2:3*iEarth).';
r0_jup   = r_body(1,3*iJupiter-2:3*iJupiter).';

r1 = norm(r0_earth);
r2 = norm(r0_jup);

u_r = r0_earth/norm(r0_earth);
u_t = v0_earth/norm(v0_earth);

offset = 1e-3*AU;
r0_sc = r0_earth + offset*u_r;

a_trans = 0.5*(norm(r0_sc) + r2);
v_dep = sqrt(G*M_sun*(2/norm(r0_sc) - 1/a_trans));

v0_sc = v_dep*u_t;
x0_sc = [r0_sc; v0_sc];

% Interpolants for the bodycentric planet trajectories
r_interp = makeInterpolants(T_body,r_body,n);

[T_sc,X_sc] = ode45(@(t,x) eom_sc_bodycentric(t,x,r_interp,m_bodies,G), ...
                    tspan,x0_sc,opts);

r_sc = X_sc(:,1:3);
v_sc = X_sc(:,4:6);

%% Plot spacecraft with planets
figure('Color','w','Name','Problem 2B — Restricted Spacecraft');
hold on; grid on; axis equal;

for i = 1:n
    ri = r_body(:,3*i-2:3*i);
    plot3(ri(:,1)/AU,ri(:,2)/AU,ri(:,3)/AU, ...
        'LineWidth',1.2, ...
        'Color',colors(i,:), ...
        'DisplayName',body_names{i});
end

plot3(r_sc(:,1)/AU,r_sc(:,2)/AU,r_sc(:,3)/AU, ...
    'k','LineWidth',2.0,'DisplayName','Spacecraft');

plot3(r_sc(1,1)/AU,r_sc(1,2)/AU,r_sc(1,3)/AU, ...
    'ko','MarkerFaceColor','k','HandleVisibility','off');

plot3(r_sc(end,1)/AU,r_sc(end,2)/AU,r_sc(end,3)/AU, ...
    'ks','MarkerFaceColor','k','HandleVisibility','off');

xlabel('x [AU]');
ylabel('y [AU]');
zlabel('z [AU]');
title('Problem 2B: Restricted spacecraft in bodycentric frame');
legend('Location','best');
view(3);

%% Distance from spacecraft to Earth and Jupiter
r_earth_sc = zeros(length(T_sc),3);
r_jup_sc   = zeros(length(T_sc),3);

for k = 1:length(T_sc)
    r_earth_sc(k,:) = r_interp{iEarth}(T_sc(k)).';
    r_jup_sc(k,:)   = r_interp{iJupiter}(T_sc(k)).';
end

d_sc_earth = vecnorm(r_sc - r_earth_sc,2,2);
d_sc_jup   = vecnorm(r_sc - r_jup_sc,2,2);

figure('Color','w','Name','Problem 2B — Spacecraft Distances');
plot(T_sc/yr,d_sc_earth/AU,'LineWidth',1.5); hold on;
plot(T_sc/yr,d_sc_jup/AU,'LineWidth',1.5);
grid on;
xlabel('Time [yr]');
ylabel('Distance [AU]');
legend('Distance to Earth','Distance to Jupiter');
title('Problem 2B: Spacecraft distance diagnostics');

fprintf('Minimum spacecraft-Earth distance = %.4e km\n',min(d_sc_earth));
fprintf('Minimum spacecraft-Jupiter distance = %.4e km\n',min(d_sc_jup));

%% ============================================================
%  LOCAL FUNCTIONS
% ============================================================

function dxdt = eom_NBP_bary(~,x,m_bodies,G)

n = numel(m_bodies);

r = x(1:3*n);
v = x(3*n+1:end);

a = zeros(3*n,1);

for i = 1:n
    ri = r(3*i-2:3*i);
    ai = zeros(3,1);

    for j = 1:n
        if j == i
            continue
        end

        rj = r(3*j-2:3*j);
        rij = ri - rj;

        ai = ai - G*m_bodies(j)*rij/norm(rij)^3;
    end

    a(3*i-2:3*i) = ai;
end

dxdt = [v; a];

end

function dxdt = eom_NBP_bodycentric(~,x,m_bodies,G)

n = numel(m_bodies);

r = x(1:3*n);
v = x(3*n+1:end);

a = zeros(3*n,1);

% Central body is fixed at the origin in this bodycentric frame
a(1:3) = [0;0;0];

for i = 2:n
    ri = r(3*i-2:3*i);

    % Direct term from central body
    ai = -G*(m_bodies(1) + m_bodies(i))*ri/norm(ri)^3;

    % Indirect terms from other planets
    for j = 2:n
        if j == i
            continue
        end

        rj = r(3*j-2:3*j);
        rij = ri - rj;

        ai = ai - G*m_bodies(j)*(rij/norm(rij)^3 + rj/norm(rj)^3);
    end

    a(3*i-2:3*i) = ai;
end

dxdt = [v; a];

end

function dxdt = eom_sc_bodycentric(t,x,r_interp,m_bodies,G)

n = numel(m_bodies);

r_sc = x(1:3);
v_sc = x(4:6);

% Direct solar term
a_sc = -G*m_bodies(1)*r_sc/norm(r_sc)^3;

% Planet perturbation terms
for j = 2:n
    rj = r_interp{j}(t);
    rj = rj(:);

    r_scj = r_sc - rj;

    a_sc = a_sc - G*m_bodies(j)*(r_scj/norm(r_scj)^3 + rj/norm(rj)^3);
end

dxdt = [v_sc; a_sc];

end

function r_cm = centerOfMassHistory(r_all,m_bodies)

n = numel(m_bodies);
Nt = size(r_all,1);
M = sum(m_bodies);

r_cm = zeros(Nt,3);

for i = 1:n
    ri = r_all(:,3*i-2:3*i);
    r_cm = r_cm + m_bodies(i)*ri;
end

r_cm = r_cm/M;

end

function P = linearMomentumHistory(v_all,m_bodies)

n = numel(m_bodies);
Nt = size(v_all,1);

P = zeros(Nt,3);

for i = 1:n
    vi = v_all(:,3*i-2:3*i);
    P = P + m_bodies(i)*vi;
end

end

function [KE,PE,E] = totalEnergyHistory(r_all,v_all,m_bodies,G)

n = numel(m_bodies);
Nt = size(r_all,1);

KE = zeros(Nt,1);
PE = zeros(Nt,1);

for k = 1:Nt
    for i = 1:n
        vi = v_all(k,3*i-2:3*i);
        KE(k) = KE(k) + 0.5*m_bodies(i)*dot(vi,vi);
    end

    for i = 1:n-1
        ri = r_all(k,3*i-2:3*i);

        for j = i+1:n
            rj = r_all(k,3*j-2:3*j);
            rij = norm(ri-rj);

            PE(k) = PE(k) - G*m_bodies(i)*m_bodies(j)/rij;
        end
    end
end

E = KE + PE;

end

function r_interp = makeInterpolants(T,r_all,n)

r_interp = cell(n,1);

for i = 1:n
    ri = r_all(:,3*i-2:3*i);

    Fx = griddedInterpolant(T,ri(:,1),'spline');
    Fy = griddedInterpolant(T,ri(:,2),'spline');
    Fz = griddedInterpolant(T,ri(:,3),'spline');

    r_interp{i} = @(t) [Fx(t); Fy(t); Fz(t)];
end

end