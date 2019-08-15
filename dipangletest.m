%% Wrapper code with changing deltas for fault model by Okada (1992), v2.0

clear all; close all;

%% Initial settings

zlevel = 0e3;         % depth of internal deformation (0 = free surface) 
mu     = 1e9;         % material properties
nu     = 0.25;

%% Define fault parameters:

L = 2e3;
W = 1e3;
U = 1.0;

phi   = deg2rad(90);
delta = deg2rad((2:4:90));

x0 = 0e3;
y0 = 0e3;
zt = 3e3;
z0_arr = zt + W*sin(delta);

% h1     = figure(1);           h2     = figure(2);
% fname1 = 'quiverplot_v1.gif';    fname2 = 'surfplot_v1.gif';

for d = 1:length(delta)
    
    z0 = z0_arr(d);
    
    res  = 100;    % resolution of box
    bpad = 5e3;    % padding of box around fault
    
    rotmat    = [sin(phi), cos(phi); -cos(phi), sin(phi)];
    midpoint  = rotmat\[x0 + 0.5*L; y0 + 0.5*W*cos(delta(d))];
    
    bmin_x = min(midpoint(1) - bpad, midpoint(1) + bpad);
    bmin_y = min(midpoint(2) - bpad, midpoint(2) + bpad);
    
    bmax_x = max(midpoint(1) - bpad, midpoint(1) + bpad);
    bmax_y = max(midpoint(2) - bpad, midpoint(2) + bpad);
    
    [xx, yy] = meshgrid(linspace(bmin_x, bmax_x, res), linspace(bmin_y, bmax_y, res));
    
    [uu, vv, ww, duu_dx, dvv_dy, dww_dz, duu_dz, dww_dx, dvv_dz, dww_dy, duu_dy, dvv_dx] = deal(zeros(size(xx)));
    
    for i = 1:res
        for j = 1:res
            [uu(i,j), vv(i,j), ww(i,j), duu_dx(i,j), dvv_dy(i,j), dww_dz(i,j), duu_dz(i,j), dww_dx(i,j), dvv_dz(i,j), ...
                dww_dy(i,j), duu_dy(i,j), dvv_dx(i,j)] = okada92_kc(x0, y0, z0, xx(i,j), yy(i,j), zlevel, L, W, U, phi, delta(d), mu, nu);
        end
    end
    
    %% Plot 1 (3D quiver plot of displacements + 3D fault plot)
    
    figure(1); clf; hold on;
    set(gca, 'FontSize', 18)
    quiver3(xx, yy, zlevel*ones(size(xx)), uu, vv, ww, 5);
    plotfault(3, x0, y0, z0, L, W, phi, delta(d));
    plot3(x0, y0, -z0, 'ko', 'markerfacecolor', 'y')
    xlim([bmin_x bmax_x]); ylim([bmin_y bmax_y]); zlim([-4000 750]);
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title(sprintf('Dip Angle: %.1f\\circ', rad2deg(delta(d))))
%     view(-65, 5)
    view(-90, 0)
    grid on
    
%     frame1 = getframe(h1);
%     im = frame2im(frame1);
%     [imind,cm] = rgb2ind(im,256);
% 
%     if d == 1
%         imwrite(imind, cm, fname1, 'gif', 'Loopcount', inf);
%     else
%         imwrite(imind, cm, fname1, 'gif', 'WriteMode', 'append');
%     end
    
    %% Plot 2 (displacement in z-direction surf plot)
    
    figure(2); clf; hold on;
    set(gca, 'FontSize', 18)
    surf(xx, yy, ww, 'FaceAlpha', .95, 'EdgeColor', [.3 .3 .3])
    xlim([bmin_x bmax_x]); ylim([bmin_y bmax_y]); zlim([-.01 .1]);
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title(sprintf('Z Displacement; Dip Angle: %.1f\\circ', rad2deg(delta(d))))
%     view(-65, 10)
    view(-90, 0)
    grid on
    
%     frame2 = getframe(h2);
%     im = frame2im(frame2);
%     [imind,cm] = rgb2ind(im,256);
% 
%     if d == 1
%         imwrite(imind, cm, fname2, 'gif', 'Loopcount', inf);
%     else
%         imwrite(imind, cm, fname2, 'gif', 'WriteMode', 'append');
%     end

end