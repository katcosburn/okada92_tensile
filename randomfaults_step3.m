%% Wrapper code with random fault parameters (step 2) for fault model by Okada (1992), v2.0

clear all; close all;

%% Define non-random variables

nevents = 5e1;      % number of intrusion events
radius  = 2e3;      % radius of disc source (m)

zlevel  = 0e3;      % depth of internal deformation (m) -> 0 is free surface
mu      = 1e9;      % shear modulus (Pa)
nu      = 0.25;     % Poisson ratio

L  = 0.5*radius;    % length along fault strike direction (m)
W  = 0.5*L;         % length along perpendicular direction to strike (m)
U  = 1.0;           % dislocation amount (m)

%% Set up circle (to be populated with faults) and gridded box

boxsize  = 1e4;      % dimension of box over which to get solution (m)
res      = 3e2;
fmax     = sqrt(L^2 + W^2);

[xx, yy] = meshgrid(linspace(-boxsize, boxsize, res), linspace(-boxsize, boxsize, res));
cinds    = circfind(xx, yy, 0, 0, radius);
ccheck   = circfind(xx, yy, 0, 0, fmax);
cweights = (.125*ismember(cinds, ccheck) + .1);

[uu, vv, ww, duu_dx, dvv_dy, dww_dz, duu_dz, dww_dx, dvv_dz, dww_dy, duu_dy, dvv_dx] = deal(zeros(size(xx, 1), size(xx, 1)));

%% For each event, populate circle with random distribtion of faults

for n = 1:nevents 
                                                                            
    fault1 = randsample(cinds, 1, true, cweights);
    
    if ismember(fault1, ccheck)
        faultlist = fault1;
    else       
        f1inds = circfind(xx, yy, xx(fault1), yy(fault1), 2*fmax);
        cfdiff = setdiff(cinds, f1inds);
        
        randflist = [];
        faultinds = [f1inds];
        
        while length(cfdiff) > length(randflist)
            
            randfault = randsample(cfdiff, 1);
            randinds  = circfind(xx, yy, xx(randfault), yy(randfault), 2*fmax);
            
            randflist = [randflist; randfault];
            faultinds = [faultinds; randinds];
            cfdiff    = setdiff(cinds, faultinds);
        end
        
        if length(randflist) > 1
            faultlist = [fault1; randsample(randflist, randi(length(randflist)))];
        else
            faultlist = [fault1; randflist];
        end
    end
    
    x0 = xx(faultlist);
    y0 = yy(faultlist);
    
    phi   = 2*pi*rand(length(faultlist), 1);
    delta = deg2rad(70) + deg2rad(40).*rand(length(faultlist),1);
    
    z0 = 1e3 + 4e3.*rand(1,1);
    
    %% Iterate through faults and get displacements/derivatives for each
    
    for f = 1:length(faultlist)
        
        [u, v, w, du_dx, dv_dy, dw_dz, du_dz, dw_dx, dv_dz, dw_dy, du_dy, dv_dx] = ...
            okada92_kc(x0(f), y0(f), z0, xx, yy, zlevel, L, W, U, phi(f), delta(f), mu, nu);
        
        uu     = uu + u;            vv     = vv + v;            ww     = ww + w;
        duu_dx = duu_dx + du_dx;    dvv_dy = dvv_dy + dv_dy;    dww_dz = dww_dz + dw_dz;
        duu_dy = duu_dy + du_dy;    dvv_dx = dvv_dx + dv_dx;    dww_dx = dww_dx + dw_dx;
        duu_dz = duu_dz + du_dz;    dvv_dz = dvv_dz + dv_dz;    dww_dy = dww_dy + dw_dy;
   
    end
end

%% Plot (displacement in z-direction surf plot)

figure(2); hold on;
set(gca, 'FontSize', 18)
surf(xx, yy, ww, 'EdgeColor', 'None') 
xlim([-boxsize boxsize]); ylim([-boxsize boxsize]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Displacement in Z-direction')
view(-25, 50)
grid on

%% Function for finding indices within a given circle of radius r, centred at x0, y0

function circinds = circfind(x, y, x0, y0, r)

    circinds = find((x-x0).^2 + (y-y0).^2 - r^2 <= 0);

end
