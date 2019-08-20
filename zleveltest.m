%% Wrapper code changing zlevel for fault model by Okada (1992), v2.0

clear all; close all;

%% Initial settings

create_gif = false;

mu = 1e9;
nu = 0.25;

%% Define fault parameters:

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

if create_gif
    
    h1 = figure(1);
    h2 = figure(2);

    fname1 = 'quiverplot_v1.gif';
    fname2 = 'surfplot_v1.gif';

end

for zl = 1:length(zlevel)

    res  = 20;
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
    
    [uu, vv, ww, duu_dx, dvv_dy, dww_dz, duu_dz, dww_dx, dvv_dz, dww_dy, duu_dy, dvv_dx] = ...
                                                    okada92_kc(x0, y0, z0, xx, yy, zlevel, L, W, U, phi, delta, mu, nu); 
    
    %% Plot 1 (3D quiver plot of displacements + 3D fault plot)
    
    figure(1); clf; hold on;
    set(gca, 'FontSize', 18)
    quiver3(xx, yy, zlevel(zl)*ones(size(xx)), uu, vv, ww, 2);
    plotfault(3, x0, y0, z0, L, W, phi, delta);
    plot3(x0, y0, -z0, 'ko', 'markerfacecolor', 'y')
    xlim([bmin_x bmax_x]); ylim([bmin_y bmax_y]); zlim([-(zlevel(end)+100) zlevel(1)+600]);
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('3D Displacement Field with Fault')
    view(-85, 10)
    grid on
    
    if create_gif
        
        frame1 = getframe(h1);
        im = frame2im(frame1);
        [imind,cm] = rgb2ind(im,256);
        
        if d == 1
            imwrite(imind, cm, fname1, 'gif', 'Loopcount', inf);
        else
            imwrite(imind, cm, fname1, 'gif', 'WriteMode', 'append');
        end
        
        
        frame2 = getframe(h2);
        im = frame2im(frame2);
        [imind,cm] = rgb2ind(im,256);
        
        if d == 1
            imwrite(imind, cm, fname2, 'gif', 'Loopcount', inf);
        else
            imwrite(imind, cm, fname2, 'gif', 'WriteMode', 'append');
        end
        
    end

end