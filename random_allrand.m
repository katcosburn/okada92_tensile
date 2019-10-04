%% Wrapper code with random parameters L, W, U, phi, delta, x0, y0

clear all; close all;

%% Define non-random variables

nevents = 10;     % number of intrusion events
radius  = 5e2;     % radius of disc source (m)

zlevel = 0e3;      % depth of internal deformation (m) -> 0 is free surface
mu     = 1e9;      % shear modulus (Pa)
nu     = 0.25;     % Poisson ratio

%% Set up circle (to be populated with dykes) and gridded box

boxsize    = 2e4;      % dimension of box over which to get solution (m)
fine_res   = 5e3;
coarse_res = 5e2;

[xx, yy] = meshgrid(linspace(-boxsize, boxsize, coarse_res), linspace(-boxsize, boxsize, coarse_res));
[xp, yp] = meshgrid(radius*linspace(-1, 1, fine_res), radius*linspace(-1, 1, fine_res));
cinds    = find(xp.^2 + yp.^2 - radius^2 <= 0);

%% For each event, populate circle with random distribution of dykes

dykelist = zeros(nevents, 1);
[uu, vv, ww, duu_dx, dvv_dy, dww_dz, duu_dz, dww_dx, dvv_dz, dww_dy, duu_dy, dvv_dx] = deal(zeros(size(xx)));

for n = 1:nevents 
                                                                            
    dyke1 = randsample(cinds, 1);  
    dykelist(n) = dyke1;
    
    x0 = xp(dyke1);
    y0 = yp(dyke1);
    z0 = 3100;
    
    L = 1000 + 4000*rand;
    W = 500  + 2500*rand;
    U = 0.25 + 0.5*abs(randn);
    
    delta = deg2rad(50 + 80*rand);
    phi   = 2*pi*rand;
           
    [u, v, w, du_dx, dv_dy, dw_dz, du_dz, dw_dx, dv_dz, dw_dy, du_dy, dv_dx] = ...
        okada92_kc(x0, y0, z0, xx, yy, zlevel, L, W, U, phi, delta, mu, nu);
    
    uu     = uu + u;            vv     = vv + v;            ww     = ww + w;
    duu_dx = duu_dx + du_dx;    dvv_dy = dvv_dy + dv_dy;    dww_dz = dww_dz + dw_dz;
    duu_dy = duu_dy + du_dy;    dvv_dx = dvv_dx + dv_dx;    dww_dx = dww_dx + dw_dx;
    duu_dz = duu_dz + du_dz;    dvv_dz = dvv_dz + dv_dz;    dww_dy = dww_dy + dw_dy;

end

%% Plot distribution of events

figure(1); hold on;
set(gca, 'FontSize', 18)
plot(radius*cosd(0:1:360), radius*sind(0:1:360), 'r--')
plot(xp(dykelist), yp(dykelist), 'ko', 'markerfacecolor', 'b')
xlabel('X (m)'); ylabel('Y (m)');
title('Random Faults on Grid')

%% Plot (displacement in z-direction surf plot)

figure(2); hold on;
set(gca, 'FontSize', 18)
surf(xx, yy, ww, 'EdgeColor', 'None') 
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Displacement in Z-direction')
view(-35, 5)
grid on
