%% Wrapper code with random fault parameters (step 1) for fault model by Okada (1992), v2.0

clear all; close all;

%% Define non-random variables

nfaults = 1e3;    % number of faults in model
zlevel  = 0e3;    % depth of internal deformation (m) -> 0 is free surface
boxsize = 1e5;    % dimension of box (m)
radius  = 6e4;    % radius of circle inside box to be populated with faults (m)
mu      = 1e9;    % shear modulus (Pa)
nu      = 0.25;   % Poisson ratio

delta = deg2rad(90);
phi   = deg2rad(90);

L  = 3e3;         % length along fault strike direction (m)
W  = 2e3;        % length along perpendicular direction to strike (m)
U  = 1.0;         % dislocation amount (m)
zt = 2e3;         % z coordinate of top left corner of fault, positive downwards (m)
z0 = zt + W*sin(delta);

%% Set up circle (to be populated with faults) inside gridded box

rotmat  = [sin(phi), cos(phi); -cos(phi), sin(phi)]; 
fmax    = max(rotmat\[L; W*cos(delta)]);

brange  = (-boxsize:fmax:boxsize+fmax);
[xg,yg] = meshgrid(brange, brange);
cinds   = find(xg.^2 + yg.^2 - radius^2 <= 0);

%% Populate circle with random distribtion of nfaults

if length(cinds) > nfaults
    frand = randsample(cinds, nfaults);
else
    error('Grid error. Please either: increase radius, decrease number of faults, or decrease fault dimensions.');
end

x0 = xg(frand);
y0 = yg(frand);

%% Plot faults on grid

figure(1); hold on;
set(gca, 'FontSize', 18)
mesh(xg, yg, zeros(size(xg))) 
plot(xg(cinds), yg(cinds), 'b.')
plot(x0, y0, 'k.')
xlim([-boxsize boxsize]); ylim([-boxsize boxsize]);
xlabel('X (m)'); ylabel('Y (m)');
title('Random Faults on Grid')

%% Iterate through faults and get displacements/derivatives for each

finefac = 5;

[xx, yy] = meshgrid(linspace(-boxsize, boxsize, finefac*length(brange)), linspace(-boxsize, boxsize, finefac*length(brange)));

[uu, vv, ww, duu_dx, dvv_dy, dww_dz, duu_dz, dww_dx, dvv_dz, dww_dy, duu_dy, dvv_dx] = deal(zeros(size(xx, 1), size(xx, 1)));

for n = 1:nfaults
    
    [u, v, w, du_dx, dv_dy, dw_dz, du_dz, dw_dx, dv_dz, dw_dy, du_dy, dv_dx] = ...
        okada92_kc(x0(n), y0(n), z0, xx, yy, zlevel, L, W, U, phi, delta, mu, nu);
    
    uu     = uu + u;            vv     = vv + v;            ww     = ww + w;
    duu_dx = duu_dx + du_dx;    dvv_dy = dvv_dy + dv_dy;    dww_dz = dww_dz + dw_dz;
    duu_dy = duu_dy + du_dy;    dvv_dx = dvv_dx + dv_dx;    dww_dx = dww_dx + dw_dx;
    duu_dz = duu_dz + du_dz;    dvv_dz = dvv_dz + dv_dz;    dww_dy = dww_dy + dw_dy;
    
end

%% Plot (displacement in z-direction surf plot)

figure(2); hold on;
set(gca, 'FontSize', 18)
surf(xx, yy, ww, 'FaceAlpha', .95) 
xlim([-boxsize boxsize]); ylim([-boxsize boxsize]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Displacement in Z-direction')
view(-65, 10)
grid on
