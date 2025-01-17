%% Wrapper code for fault model by Okada (1992), v2.0

clear all; close all;

%% Initial settings

zlevel = 0e3;         % depth of internal deformation (m) -> 0 is free surface
mu     = 1e9;         % material properties, shear modulus (Pa)
nu     = 0.25;        % material properties, Poisson ratio

%% Define fault parameters:
%
%  L      = length along fault strike direction (m)
%  W      = length along perpendicular direction to strike (m)
%  U      = dislocation amount (m)
%  phi    = strike angle (rad) 
%  delta  = dip angle (rad)
%  x0, y0 = x and y coordinates of bottom left corner of fault (m)
%  zt     = z coordinate of top left corner of fault, positive downwards (m)
%  z0     = z coordinate of bottom left corner of fault (m)

L = 20e3;
W = 0.5*L;
U = 1.0;

phi   = deg2rad(90);
delta = deg2rad(90);

x0 = 0e3;
y0 = 0e3;
zt = 2e3;
z0 = zt + W*sin(delta);

% z0 = 4e3;
% zt = z0 - W*sin(delta);

%% Set up box of coordinates over which to get displacements

res  = 1e3;    % resolution of box 
bpad = 1e4;    % padding of box around fault

rotmat   = [sin(phi), cos(phi); -cos(phi), sin(phi)];
midpoint = rotmat\[0.5*L; 0.5*W*cos(delta)] + [x0; y0];

bmin_x = min(midpoint(1) - bpad, midpoint(1) + bpad);  
bmin_y = min(midpoint(2) - bpad, midpoint(2) + bpad);

bmax_x = max(midpoint(1) - bpad, midpoint(1) + bpad);  
bmax_y = max(midpoint(2) - bpad, midpoint(2) + bpad);

[xx, yy] = meshgrid(linspace(bmin_x, bmax_x, res), linspace(bmin_y, bmax_y, res));

%% Displacement from okada85 and okada92 

[uu, vv, ww, duu_dx, dvv_dy, dww_dz, duu_dz, dww_dx, dvv_dz, dww_dy, duu_dy, dvv_dx] = ...
    okada92_kc(x0, y0, z0, xx, yy, zlevel, L, W, U, phi, delta, mu, nu); 

%% Plot 1 (2D quiver plot of displacements + 2D fault plot)

figure(1); hold on;
set(gca, 'FontSize', 18)
quiver(xx, yy, uu, vv); hold on;
plotfault(3, x0, y0, z0, L, W, phi, delta);
xlim([bmin_x bmax_x]); ylim([bmin_y bmax_y]);
xlabel('X (m)'); ylabel('Y (m)');
title('2D Displacement Field with Fault')
grid on

%% Plot 2 (3D quiver plot of displacements + 3D fault plot)

figure(2); hold on;
set(gca, 'FontSize', 18)
quiver3(xx, yy, zlevel*ones(size(xx)), uu, vv, ww, 2); 
plotfault(3, x0, y0, z0, L, W, phi, delta); 
plot3(x0, y0, -z0, 'ko', 'markerfacecolor', 'y')
xlim([bmin_x bmax_x]); ylim([bmin_y bmax_y]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('3D Displacement Field with Fault')
view(-65, 10)
grid on

%% Plot 3 (displacement in z-direction surf plot)

figure(3); hold on;
set(gca, 'FontSize', 18)
surf(xx, yy, ww, 'EdgeColor', 'None') 
xlim([bmin_x bmax_x]); ylim([bmin_y bmax_y]);
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
surf(xx, yy, area_dil, 'EdgeColor', 'None') 
xlabel('X (m)'); zlabel('du/dx + dv/dy');
view(-65, 10)
grid on

subplot(3,1,2); hold on;
set(gca, 'FontSize', 14)
surf(xx, yy, sstrain1, 'EdgeColor', 'None') 
xlabel('X (m)'); zlabel('du/dx - dv/dy');
view(-65, 10)
grid on

subplot(3,1,3); hold on;
set(gca, 'FontSize', 14)
surf(xx, yy, sstrain2, 'EdgeColor', 'None') 
xlabel('X (m)'); zlabel('du/dx + dv/dx');
view(-65, 10)
grid on