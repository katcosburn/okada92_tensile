function plotfault(dim, x0, y0, z0, L, W, phi, delta) 

rotmatinv = [sin(phi), -cos(phi); cos(phi), sin(phi)]; 

dz = W*sin(delta);

z1 = z0 - dz;
xn = x0 + L;
yn = y0 + W*cos(delta);

xy1 = rotmatinv*[xn; y0];
xy2 = rotmatinv*[xn; yn];
xy3 = rotmatinv*[x0; yn];

if dim == 2
    
    X = [xy3(1) xy2(1)];
    Y = [xy3(2) xy2(2)];
    
    plot(X, Y, 'r-', 'Linewidth', 2)
    
elseif dim == 3
    
    X = [x0 xy1(1) xy2(1)  xy3(1) x0];
    Y = [y0 xy1(2) xy2(2)  xy3(2) y0];
    Z = [z0 z0     z1      z1     z0];
    
    plot3(X, Y, -Z); hold on;
    
else
    error('Please enter either 2 or 3 for dimension variable.') 
end

end