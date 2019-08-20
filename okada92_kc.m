function [u, v, w, dudx, dvdy, dwdz, dudz, dwdx, dvdz, dwdy, dudy, dvdx] = okada92_kc(x0, y0, z0, x, y, z, L, W, U, phi, delta, mu, nu)

lambda = 2*mu*nu/(1-2*nu);
alpha  = (lambda + mu)/(lambda + 2*mu);

% From the reference system of Fig. 3 (Okada, 1992): 
% -> centre at (x0,y0) and rotate by strike angle phi

rotmat = [sin(phi), cos(phi); -cos(phi), sin(phi)];
xyrot  = rotmat*[x(:) - x0, y(:) - y0]';
xc     = xyrot(1,:);
yc     = xyrot(2,:);

%% Collect fA, fB, fC functions from Table 6

d = [z0 - z, z0 + z];

for i = 1:length(d)
    
    p = yc*cos(delta) + d(i)*sin(delta);
    q = yc*sin(delta) - d(i)*cos(delta);
    
    if i == 1
        
        [fAa, dfAa_dx, dfAa_dy, dfAa_dz] = f_dfxyz('A', xc,   p,   q, z, delta, alpha);
        [fAb, dfAb_dx, dfAb_dy, dfAb_dz] = f_dfxyz('A', xc,   p-W, q, z, delta, alpha);
        [fAc, dfAc_dx, dfAc_dy, dfAc_dz] = f_dfxyz('A', xc-L, p,   q, z, delta, alpha);
        [fAd, dfAd_dx, dfAd_dy, dfAd_dz] = f_dfxyz('A', xc-L, p-W, q, z, delta, alpha);
        
        [fBa, dfBa_dx, dfBa_dy, dfBa_dz] = f_dfxyz('B', xc,   p,   q, z, delta, alpha);
        [fBb, dfBb_dx, dfBb_dy, dfBb_dz] = f_dfxyz('B', xc,   p-W, q, z, delta, alpha);
        [fBc, dfBc_dx, dfBc_dy, dfBc_dz] = f_dfxyz('B', xc-L, p,   q, z, delta, alpha);
        [fBd, dfBd_dx, dfBd_dy, dfBd_dz] = f_dfxyz('B', xc-L, p-W, q, z, delta, alpha);
        
        [fCa, dfCa_dx, dfCa_dy, dfCa_dz] = f_dfxyz('C', xc,   p,   q, z, delta, alpha);
        [fCb, dfCb_dx, dfCb_dy, dfCb_dz] = f_dfxyz('C', xc,   p-W, q, z, delta, alpha);
        [fCc, dfCc_dx, dfCc_dy, dfCc_dz] = f_dfxyz('C', xc-L, p,   q, z, delta, alpha);
        [fCd, dfCd_dx, dfCd_dy, dfCd_dz] = f_dfxyz('C', xc-L, p-W, q, z, delta, alpha);
        
    else
        
        [fATa, dfATa_dx, dfATa_dy, dfATa_dz] = f_dfxyz('A', xc,   p,   q, -z, delta, alpha);
        [fATb, dfATb_dx, dfATb_dy, dfATb_dz] = f_dfxyz('A', xc,   p-W, q, -z, delta, alpha);
        [fATc, dfATc_dx, dfATc_dy, dfATc_dz] = f_dfxyz('A', xc-L, p,   q, -z, delta, alpha);
        [fATd, dfATd_dx, dfATd_dy, dfATd_dz] = f_dfxyz('A', xc-L, p-W, q, -z, delta, alpha);
        
    end
end
     
%% From Equation 15 (Chinnery's expression)

A  = fAa  - fAb  - fAc  + fAd;
AT = fATa - fATb - fATc + fATd;
B  = fBa  - fBb  - fBc  + fBd;
C  = fCa  - fCb  - fCc  + fCd;
    
dAdx  = dfAa_dx  - dfAb_dx  - dfAc_dx  + dfAd_dx;
dATdx = dfATa_dx - dfATb_dx - dfATc_dx + dfATd_dx;
dBdx  = dfBa_dx  - dfBb_dx  - dfBc_dx  + dfBd_dx;
dCdx  = dfCa_dx  - dfCb_dx  - dfCc_dx  + dfCd_dx;

dAdy  = dfAa_dy  - dfAb_dy  - dfAc_dy  + dfAd_dy;
dATdy = dfATa_dy - dfATb_dy - dfATc_dy + dfATd_dy;
dBdy  = dfBa_dy  - dfBb_dy  - dfBc_dy  + dfBd_dy;
dCdy  = dfCa_dy  - dfCb_dy  - dfCc_dy  + dfCd_dy;

dAdz  = dfAa_dz  - dfAb_dz  - dfAc_dz  + dfAd_dz;
dATdz = dfATa_dz - dfATb_dz - dfATc_dz + dfATd_dz;
dBdz  = dfBa_dz  - dfBb_dz  - dfBc_dz  + dfBd_dz;
dCdz  = dfCa_dz  - dfCb_dz  - dfCc_dz  + dfCd_dz;

%% Calculate displacements in x, y, z directions using convention in Eqs. 16, 17

uP = (0.5*U/pi)*(A(1,:) - AT(1,:) + B(1,:) + z*C(1,:));
vP = (0.5*U/pi)*((A(2,:) - AT(2,:) + B(2,:) + z*C(2,:))*cos(delta) - (A(3,:) - AT(3,:) + B(3,:) + z*C(3,:))*sin(delta));
wP = (0.5*U/pi)*((A(2,:) - AT(2,:) + B(2,:) - z*C(2,:))*sin(delta) + (A(3,:) - AT(3,:) + B(3,:) - z*C(3,:))*cos(delta));

uvrot = rotmat\[uP; vP];

u  = reshape(uvrot(1,:), size(x));
v  = reshape(uvrot(2,:), size(y));
w  = reshape(wP, size(x));

%% Calculate x, y, z derivatives in x, y, z directions using convention in Eqs. 16, 17

rotmat_dvs = [sin(phi)^2  ,    cos(phi)^2    ,  0,  sin(phi)*cos(phi),    0     ,  sin(phi)*cos(phi),     0   ,     0    ,    0     ; ...
              cos(phi)^2  ,    sin(phi)^2    ,  0, -sin(phi)*cos(phi),    0     , -sin(phi)*cos(phi),     0   ,     0    ,    0     ; ...
                  0       ,        0         ,  1,          0        ,    0     ,          0        ,     0   ,     0    ,    0     ; ...
        -sin(phi)*cos(phi), sin(phi)*cos(phi),  0,     sin(phi)^2    ,    0     ,     -cos(phi)^2   ,     0   ,     0    ,    0     ; ...
                  0       ,        0         ,  0,          0        , sin(phi) ,          0        , cos(phi),     0    ,    0     ; ...
        -sin(phi)*cos(phi), sin(phi)*cos(phi),  0,    -cos(phi)^2    ,    0     ,      sin(phi)^2   ,     0   ,     0    ,    0     ; ...
                  0       ,        0         ,  0,          0        , -cos(phi),          0        , sin(phi),     0    ,    0     ; ...
                  0       ,        0         ,  0,          0        ,    0     ,          0        ,     0   ,  sin(phi), cos(phi) ; ...
                  0       ,        0         ,  0,          0        ,    0     ,          0        ,     0   , -cos(phi), sin(phi)];

dudxP = (0.5*U/pi)*(dAdx(1,:) - dATdx(1,:) + dBdx(1,:) + z*dCdx(1,:));
dvdxP = (0.5*U/pi)*((dAdx(2,:) - dATdx(2,:) + dBdx(2,:) + z*dCdx(2,:))*cos(delta) - (dAdx(3,:) - dATdx(3,:) + dBdx(3,:) + z*dCdx(3,:))*sin(delta));
dwdxP = (0.5*U/pi)*((dAdx(2,:) - dATdx(2,:) + dBdx(2,:) - z*dCdx(2,:))*sin(delta) + (dAdx(3,:) - dATdx(3,:) + dBdx(3,:) - z*dCdx(3,:))*cos(delta));

dudyP = (0.5*U/pi)*(dAdy(1,:) - dATdy(1,:) + dBdy(1,:) + z*dCdy(1,:));
dvdyP = (0.5*U/pi)*((dAdy(2,:) - dATdy(2,:) + dBdy(2,:) + z*dCdy(2,:))*cos(delta) - (dAdy(3,:) - dATdy(3,:) + dBdy(3,:) + z*dCdy(3,:))*sin(delta));
dwdyP = (0.5*U/pi)*((dAdy(2,:) - dATdy(2,:) + dBdy(2,:) - z*dCdy(2,:))*sin(delta) + (dAdy(3,:) - dATdy(3,:) + dBdy(3,:) - z*dCdy(3,:))*cos(delta));

dudzP = (0.5*U/pi)*(dAdz(1,:) + dATdz(1,:) + dBdz(1,:) + C(1,:) + z*dCdz(1,:));
dvdzP = (0.5*U/pi)*((dAdz(2,:) + dATdz(2,:) + dBdz(2,:) + C(2,:) + z*dCdz(2,:))*cos(delta) - (dAdz(3,:) + dATdz(3,:) + dBdz(3,:) + C(3,:) + z*dCdz(3,:))*sin(delta));
dwdzP = (0.5*U/pi)*((dAdz(2,:) + dATdz(2,:) + dBdz(2,:) - C(2,:) - z*dCdz(2,:))*sin(delta) + (dAdz(3,:) + dATdz(3,:) + dBdz(3,:) - C(3,:) - z*dCdz(3,:))*cos(delta));

dvs_rot = rotmat_dvs\[dudxP; dvdyP; dwdzP; dudyP; dudzP; dvdxP; dvdzP; dwdxP; dwdyP];

dudx = reshape(dvs_rot(1,:), size(x));
dvdy = reshape(dvs_rot(2,:), size(y));
dwdz = reshape(dvs_rot(3,:), size(x));
dudy = reshape(dvs_rot(4,:), size(y));
dudz = reshape(dvs_rot(5,:), size(x));
dvdx = reshape(dvs_rot(6,:), size(y));
dvdz = reshape(dvs_rot(7,:), size(x));
dwdx = reshape(dvs_rot(8,:), size(y));
dwdy = reshape(dvs_rot(9,:), size(x));

end