clear all
close all
clc

%% DATA
% eart orbit
orbit_E = struct('mu', astroConstants(13), 'a', astroConstants(23) + 2000, ...
                  'e', 0.1, 'inc', deg2rad(10), 'G', astroConstants(1), ...
                  'theta0', 0, 'R', astroConstants(23));

orbit_E.n = sqrt(orbit_E.mu/orbit_E.a^3);
orbit_E.T = 2*pi / orbit_E.n;
orbit_E.w = 2*pi/(23*3600+56*60+4.09);   

% solar orbit
orbit_S = struct('mu', astroConstants(4), 'a', astroConstants(2), 'e', 0.0167, ...
                  'epsilon', deg2rad(23.45), 'G', astroConstants(1), ...
                  'theta0', 0  );

orbit_S.n =  sqrt(orbit_S.mu/orbit_S.a^3);


radiation = 1358; %(620 albedo@h +157 earth@h)
density = 2.72e-12;

m_t = 5.9722e24; %earth mass

% Spacecraft (inertia, initial angle, torque)
cubesat = struct('m', .8 * 12, 'a', .2, 'b', .2, 'c', .3405);

Ix = 1/12 * cubesat.m * (cubesat.b^2 + cubesat.c^2);
Iy = 1/12 * cubesat.m * (cubesat.a^2 + cubesat.c^2);
Iz = 1/12 * cubesat.m * (cubesat.a^2 + cubesat.b^2);

I = [Ix; Iy; Iz];
w0 = [0; 0; orbit_E.n];
M = [0;0;0];

% Wheel (inertia, direction, initial angle, torque)
Ir = 0.005;
k = [0;0;1];
wr0 = 0;
Mr = 0;

% Body definition, cubesat (6 surf) with two panels (4 surf) on the side
pos = zeros(3,10);
pos_CG = 0.015*[1;1;1]; % Displacement of CG from geometric center
n1=[1;0;0]; n2=[0;1;0]; n3=[0;0;1];

pos(:,1)=.1*n1; pos(:,4)=-.1*n1; pos(:,2)=.1*n2; pos(:,5)=-.1*n2; pos(:,3)=.15*n3; pos(:,6)=-.15*n3;
pos(:,7)=-.4*n2; pos(:,8)=-.4*n2; pos(:,9)=.4*n2; pos(:,10)=.4*n2;
r = pos-pos_CG;

solarPanel = struct( 'r', r, 'P', radiation/(3e8), ...
                    'rho_d', .1, 'rho_s', [.5, .5, .5, .5, .5 .5, .1, .1, .1, .1], ...
                    'area', [.06 .06 .04 .06 .06 .04 .12 .12 .12 .12]);

% Magnetic perturbation
m_par = [0.01;0.05;0.01];

mag_tilt = deg2rad(11.5);

H0 = sqrt((-29615e-9)^2+(-1728e-9)^2+(5186e-9)^2);


% Initial DCM
A0=[1 0 0; 0 1 0; 0 0 1];

% Sensor
gyro = struct('c', 15e-3, 'ARW', 0.15, 'sigma_b', deg2rad(0.3) / 3600, ...
              'freq', 262, 'accuracy', 2.4e-9); % the accurancy is taken from wikipedia [rad]
gyro.A_eps = [1, 0.2e-3, 0.2e-3;
                0.2e-3, 1, 0.2e-3
                0.2e-3, 0.2e-3, 1];

% Attitude Determination
sun_sensor    = struct('dvst', 3*0.2 ,'accuracy', deg2rad(0.2));
earth_horizon = struct('dvst', 3*0.2 ,'accuracy', deg2rad(0.2));

% PD parameters
K_p = Iz * (20 * orbit_E.n)^2;   % we assumed w = 20*n
K_d = 40 * Iz * 0.3 * orbit_E.n; % we assumed csi = 0.3 

%% COMPUTATION
% simu = sim("kde");





