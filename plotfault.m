function plotfault(dim, x0, y0, z0, L, W, phi, delta) 

p1 = [x0, y0, z0] + L*[sin(phi), cos(phi), 0];
p2 = p1 + W*[-cos(delta)*cos(phi), cos(delta)*sin(phi), -sin(delta)];
p3 = [x0, y0, z0] + (p2 - p1);

if dim == 2
    
    X = [x0, p1(1), p2(1), p3(1), x0];
    Y = [y0, p1(2), p2(2), p3(2), y0];
    
    plot(X, Y); hold on;
    
elseif dim == 3
    
    X = [x0, p1(1), p2(1), p3(1), x0];
    Y = [y0, p1(2), p2(2), p3(2), y0];
    Z = [z0, p1(3), p2(3), p3(3), z0];
    
    plot3(X, Y, -Z); hold on;
    
else
    error('Please enter either 2 or 3 for dimension variable.') 
end

end