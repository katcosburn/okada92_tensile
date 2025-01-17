%% Wrapper code with changing deltas for fault model by Okada (1992), v2.0

clear all; close all;

%% Initial settings

create_gif = false;

zlevel = 0e3;
mu     = 1e9;
nu     = 0.25;

%% Define fault parameters:

L = 2e3;
W = 1e3;
U = 1.0;

phi   = deg2rad((0:10:360));
delta = deg2rad(90);

x0 = 0e3;
y0 = 0e3;
zt = 3e3;
z0 = zt + W*sin(delta);

if create_gif
    
    h1 = figure(1);
    h2 = figure(2);

    fname1 = 'quiverplot_v1.gif';
    fname2 = 'surfplot_v1.gif';

end

for p = 1:length(phi)
    
    res  = 100;    % resolution of box
    bpad = 1e4;    % padding of box around fault
    
    rotmat   = [sin(phi(p)), cos(phi(p)); -cos(phi(p)), sin(phi(p))];
    midpoint = rotmat\[0.5*L; 0.5*W*cos(delta)] + [x0; y0];
    
    bmin_x = min(midpoint(1) - bpad, midpoint(1) + bpad);
    bmin_y = min(midpoint(2) - bpad, midpoint(2) + bpad);
    
    bmax_x = max(midpoint(1) - bpad, midpoint(1) + bpad);
    bmax_y = max(midpoint(2) - bpad, midpoint(2) + bpad);
    
    [xx, yy] = meshgrid(linspace(bmin_x, bmax_x, res), linspace(bmin_y, bmax_y, res));
    
    [uu, vv, ww, duu_dx, dvv_dy, dww_dz, duu_dz, dww_dx, dvv_dz, dww_dy, duu_dy, dvv_dx] = ...
                                                    okada92_kc(x0, y0, z0, xx, yy, zlevel, L, W, U, phi(p), delta, mu, nu); 
    
    %% Plot 1 (3D quiver plot of displacements + 3D fault plot)
    
    figure(1); clf; hold on;
    set(gca, 'FontSize', 18)
    quiver3(xx, yy, zlevel*ones(size(xx)), uu, vv, ww, 5);
    plotfault(3, x0, y0, z0, L, W, phi(p), delta);
    plot3(x0, y0, -z0, 'ko', 'markerfacecolor', 'y')
    xlim([bmin_x bmax_x]); ylim([bmin_y bmax_y]); zlim([-4e3 1e3]);
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title(sprintf('Strike Angle: %.1f\\circ', rad2deg(phi(p))))
    view(40, 10)
    grid on
    
    %% Plot 2 (displacement in z-direction surf plot)
    
    figure(2); clf; hold on;
    set(gca, 'FontSize', 18)
    surf(xx, yy, ww, 'FaceAlpha', .95, 'EdgeColor', [.3 .3 .3])
    xlim([bmin_x bmax_x]); ylim([bmin_y bmax_y]); zlim([-.01 .02]);
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title(sprintf('Strike Angle: %.1f\\circ', rad2deg(phi(p))))
    view(40, 10)
    grid on
    
    if create_gif
    
        frame1 = getframe(h1);
        im = frame2im(frame1);
        [imind,cm] = rgb2ind(im,256);
        
        if p == 1
            imwrite(imind, cm, fname1, 'gif', 'Loopcount', inf);
        else
            imwrite(imind, cm, fname1, 'gif', 'WriteMode', 'append');
        end
        
        
        frame2 = getframe(h2);
        im = frame2im(frame2);
        [imind,cm] = rgb2ind(im,256);
        
        if p == 1
            imwrite(imind, cm, fname2, 'gif', 'Loopcount', inf);
        else
            imwrite(imind, cm, fname2, 'gif', 'WriteMode', 'append');
        end
    
    end

end