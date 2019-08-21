function plotfault(dim, x0, y0, z0, L, W, phi, delta) 

rotmat  = [sin(phi), cos(phi); -cos(phi), sin(phi)];
LWprime = rotmat\[L; W*cos(delta)];

dz = W*sin(delta);
z1 = z0 - dz;

xy1 = LWprime + [x0; y0];

if dim == 2
    
    X = [xy0(1) xy1(1)];
    Y = [xy0(2) xy1(2)];
    
    plot(X, Y, 'r-', 'Linewidth', 2)
    
elseif dim == 3
    
    X = [x0 xy1(1) xy1(1) x0 x0];
    Y = [y0 xy1(2) xy1(2) y0 y0];
    Z = [z0 z0     z1     z1 z0];
    
    plot3(X, Y, -Z); hold on;
    
else
    error('Please enter either 2 or 3 for dimension variable.') 
end

end