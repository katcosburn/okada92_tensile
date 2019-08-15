%% Wrapper code changing zlevel for fault model by Okada (1992), v2.0

clear all; close all;

%% Initial settings

mu = 1e9;         % material properties
nu = 0.25;

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

L = 2e3;
W = 1e3;
U = 1.0;

phi   = deg2rad(90);
delta = deg2rad(90);

x0 = 0e3;
y0 = 0e3;
zt = 1e3;
z0 = zt + W*sin(delta);

zlevel = -1*[(0:20:zt-1), (zt+1:20:z0-1), (z0+1:20:z0+500)];

for zl = 1:length(zlevel)

    res  = 40;
    bpad = 5e3;
    
    rotmat   = [sin(phi), cos(phi); -cos(phi), sin(phi)];
    midpoint = rotmat\[x0 + 0.5*L; y0 + 0.5*W*cos(delta)];    
    endpoint = rotmat\[x0 + L; y0 + W*cos(delta)];
    
    bmin_x = min(midpoint(1) - bpad, midpoint(1) + bpad);
    bmin_y = min(midpoint(2) - bpad, midpoint(2) + bpad);
    
    bmax_x = max(midpoint(1) - bpad, midpoint(1) + bpad);
    bmax_y = max(midpoint(2) - bpad, midpoint(2) + bpad);
    
    xrange = linspace(bmin_x, bmax_x, res);
    yrange = linspace(bmin_y, bmax_y, res);
                  
    [xx, yy] = meshgrid(xrange, yrange);
    
    [uu, vv, ww, duu_dx, dvv_dy, dww_dz, duu_dz, dww_dx, dvv_dz, dww_dy, duu_dy, dvv_dx] = deal(zeros(size(xx)));
    
    %% Displacement from okada85 and okada92
    
    for i = 1:res
        for j = 1:res
            [uu(i,j), vv(i,j), ww(i,j), duu_dx(i,j), dvv_dy(i,j), dww_dz(i,j), duu_dz(i,j), dww_dx(i,j), dvv_dz(i,j), ...
                dww_dy(i,j), duu_dy(i,j), dvv_dx(i,j)] = okada92_kc(x0, y0, z0, xx(i,j), yy(i,j), zlevel(zl), L, W, U, phi, delta, mu, nu);
        end
    end
    
    %% Plot 1 (3D quiver plot of displacements + 3D fault plot)
    
    figure(1); clf; hold on;
    set(gca, 'FontSize', 18)
    quiver3(xx, yy, zlevel(zl)*ones(size(xx)), uu, vv, ww, 2);
    plotfault(3, x0, y0, z0, L, W, phi, delta);
    plot3(x0, y0, -z0, 'ko', 'markerfacecolor', 'y')
    xlim([bmin_x bmax_x]); ylim([bmin_y bmax_y]); zlim([-(z0+600) 600]);
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('3D Displacement Field with Fault')
    view(-85, 10)
    grid on

end