%% Wrapper code with random fault parameters (step 1) for fault model by Okada (1992), v2.0

clear all; close all;

%% Define non-random variables

nfaults = 5e3;    % number of faults in model
zlevel  = 0e3;    % depth of internal deformation (m) -> 0 is free surface
boxsize = 7e4;    % dimension of box (m)
radius  = 5e4;    % radius of circle inside box to be populated with faults (m)
mu      = 1e9;    % shear modulus (Pa)
nu      = 0.25;   % Poisson ratio

delta = deg2rad(90);
phi   = deg2rad(90);

L  = 1e3;         % length along fault strike direction (m)
W  = 5e2;        % length along perpendicular direction to strike (m)
U  = 1.0;         % dislocation amount (m)
zt = 3e3;         % z coordinate of top left corner of fault, positive downwards (m)
z0 = zt + W*sin(delta);

%% Set up circle (to be populated with faults) inside gridded box

rotmat  = [sin(phi), cos(phi); -cos(phi), sin(phi)]; 
fmax    = max(rotmat\[L; W*cos(delta)]);
% fmax    = sqrt(L^2 + W^2);

brange  = (-boxsize:fmax:boxsize+fmax);
[xg,yg] = meshgrid(brange, brange);
cinds   = find(xg.^2 + yg.^2 - radius^2 <= 0);

%% Populate circle with random distribtion of nfaults

if length(cinds) > nfaults
    frand = randsample(cinds, nfaults);
else
    error('Grid error. Please either: increase boxsize, decrease number of faults, or decrease fault dimensions.');
end

x0 = xg(frand);
y0 = yg(frand);

%% Plot fault placements inside circle on grid

figure(1); hold on;
set(gca, 'FontSize', 18)
mesh(xg, yg, zeros(size(xg)), 'LineWidth', 0.1);
plot(xg(cinds), yg(cinds), 'b.')
plot(x0, y0, 'k.')
xlim([-boxsize, boxsize+fmax]);  ylim([-boxsize, boxsize+fmax]);
title('Fault Placement on Grid')
xlabel('X (m)'); ylabel('Y (m)');

%% Iterate through faults and get displacements/derivatives for each

[xx, yy] = meshgrid(linspace(-boxsize, boxsize, 2*length(brange)), linspace(-boxsize, boxsize, 2*length(brange)));

[u, v, w, du_dx, dv_dy, dw_dz, du_dz, dw_dx, dv_dz, dw_dy, du_dy, dv_dx] = ...
    deal(zeros(size(xx, 1), size(xx, 1)));

[uu, vv, ww, duu_dx, dvv_dy, dww_dz, duu_dz, dww_dx, dvv_dz, dww_dy, duu_dy, dvv_dx] = ...
    deal(zeros(size(xx, 1), size(xx, 1)));

for n = 1:nfaults
    
%     figure(2); hold on;
%     set(gca, 'FontSize', 18)
%     plotfault(3, x0(n), y0(n), z0, L, W, phi, delta); 
%     plot3(x0(n), y0(n), -z0, 'ko', 'markerfacecolor', 'y')
%     view(-65, 10)
%     grid on
    
    for i = 1:size(xx, 1)
        for j = 1:size(xx, 1)
            [u(i,j), v(i,j), w(i,j), du_dx(i,j), dv_dy(i,j), dw_dz(i,j), du_dz(i,j), ...
                dw_dx(i,j), dv_dz(i,j), dw_dy(i,j), du_dy(i,j), dv_dx(i,j)] = ...
                okada92_kc(x0(n), y0(n), z0, xx(i,j), yy(i,j), zlevel, L, W, U*1e3, phi, delta, mu, nu);
        end
    end
    
    uu     = uu + u;            vv     = vv + v;            ww     = ww + w;
    duu_dx = duu_dx + du_dx;    dvv_dy = dvv_dy + dv_dy;    dww_dz = dww_dz + dw_dz;
    duu_dy = duu_dy + du_dy;    dvv_dx = dvv_dx + dv_dx;    dww_dx = dww_dx + dw_dx;
    duu_dz = duu_dz + du_dz;    dvv_dz = dvv_dz + dv_dz;    dww_dy = dww_dy + dw_dy;
    
end

%% Plot 2 (3D quiver plot of displacements + 3D fault plot)

figure(2); hold on;
quiver3(xx, yy, zlevel*ones(size(xx)), uu, vv, ww, 2); 
xlim([-boxsize boxsize]); ylim([-boxsize boxsize]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('3D Displacement Field with Fault')
view(-65, 10)
grid on

%% Plot 3 (displacement in z-direction surf plot)

figure(3); hold on;
set(gca, 'FontSize', 18)
surf(xx, yy, ww, 'FaceAlpha', .95, 'EdgeColor', 'None') 
xlim([-boxsize boxsize]); ylim([-boxsize boxsize]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Displacement in Z-direction')
view(-65, 10)
grid on

%% Plot 4 (strain plot)

area_dil = duu_dx + dvv_dy;
sstrain1 = duu_dx - dvv_dy;
sstrain2 = duu_dy + dvv_dx;

figure(4); hold on;

subplot(3,1,1); hold on;
set(gca, 'FontSize', 14)
surf(xx, yy, area_dil, 'FaceAlpha', .95, 'EdgeColor', 'None') 
xlabel('X (m)'); zlabel('du/dx + dv/dy');
view(-65, 10)
grid on

subplot(3,1,2); hold on;
set(gca, 'FontSize', 14)
surf(xx, yy, sstrain1, 'FaceAlpha', .95, 'EdgeColor', 'None') 
xlabel('X (m)'); zlabel('du/dx - dv/dy');
view(-65, 10)
grid on

subplot(3,1,3); hold on;
set(gca, 'FontSize', 14)
surf(xx, yy, sstrain2, 'FaceAlpha', .95, 'EdgeColor', 'None') 
xlabel('X (m)'); zlabel('du/dx + dv/dx');
view(-65, 10)
grid on
