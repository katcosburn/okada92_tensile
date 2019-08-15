function [f, dfdx, dfdy, dfdz] = f_dfxyz(letter, xi, eta, q, z, delta, alpha)

epsilon = 1e-8;      % criterion for determining if a variable is zero

R       = sqrt(xi^2 + eta^2 + q^2);
X       = sqrt(xi^2 + q^2);
ytilde  = eta*cos(delta) + q*sin(delta);
dtilde  = eta*sin(delta) - q*cos(delta);
ctilde  = dtilde + z;

Rdtilde = R + dtilde;
alpha2  = 1 - alpha;
D11     = 1/(R*Rdtilde);

%% Check singularlity conditions

% condition (i)

check_i = abs(q) > epsilon;
theta   = check_i*atan(xi*eta/(q*R));

% condition (ii)

check_ii = abs(xi) > epsilon;

% condition (iii)

check_iii = abs(R + xi) > epsilon;
ln_Rxi    = ~check_iii*(-1)*log(R - xi) + check_iii*log(R + xi);
X11       = check_iii/R/(R + xi);
X32       = check_iii*(2*R + xi)/R^3/(R + xi)^2;
X53       = check_iii*(8*R^2 + 9*R*xi + 3*xi^2)/R^5/(R + xi)^3;

% condition (iv)

check_iv = abs(R + eta) > epsilon;
ln_Reta  = ~check_iv*(-1)*log(R - eta) + check_iv*log(R + eta);
Y11      = check_iv/R/(R + eta);
Y32      = check_iv*(2*R + eta)/R^3/(R + eta)^2;
Y53      = check_iv*(8*R^2 + 9*R*eta + 3*eta^2)/R^5/(R + eta)^3;

%% Define other variables from equation 14

h   = q*cos(delta) - z;
Z32 = sin(delta)/R^3 - h*Y32;
Z53 = 3*sin(delta)/R^5 - h*Y53;

Y0  = Y11 - Y32*xi^2;
Z0  = Z32 - Z53*xi^2;

%% Define I parameters (for displacement equations, Table 6)

if cos(delta) < epsilon
    
    I3 = 0.5*(eta/Rdtilde + ytilde*q/Rdtilde^2 - ln_Reta);
    I4 = check_ii*0.5*xi*ytilde/Rdtilde^2;
    
else
    
    I3 = ytilde/cos(delta)/Rdtilde - (ln_Reta - sin(delta)*log(Rdtilde))/cos(delta)^2;
    I4 = check_ii*(xi*tan(delta)/Rdtilde + 2*atan((eta*(X + q*cos(delta)) + X*(R + X)*sin(delta))...
                                                    /(xi*(R + X)*cos(delta)))/cos(delta)^2);  
end

I1 = -(xi/Rdtilde)*cos(delta) - I4*sin(delta);
I2 = log(Rdtilde) + I3*sin(delta);

%% Define J and K parameters (for x-derivative equations, Table 7)

J2 = xi*ytilde*D11/Rdtilde;
J5 = -D11*(dtilde + (ytilde^2)/Rdtilde);

if cos(delta) < epsilon
   
    K1 = xi*q*D11/Rdtilde;
    K3 = sin(delta)*(D11*xi^2 - 1.0)/Rdtilde;
    J3 = -xi*(D11*q^2 - 0.5)/Rdtilde^2;
    J6 = -ytilde*(D11*xi^2 - 0.5)/Rdtilde^2;
    
else
    
    K1 = xi*(D11 - Y11*sin(delta))/cos(delta);
    K3 = (q*Y11 - ytilde*D11)/cos(delta);
    J3 = (K1 - J2*sin(delta))/cos(delta);
    J6 = (K3 - J5*sin(delta))/cos(delta);
    
end

K2 = 1/R + K3*sin(delta);
K4 = xi*Y11*cos(delta) - K1*sin(delta);
J1 = J5*cos(delta) - J6*sin(delta);
J4 = -xi*Y11 - J2*cos(delta) + J3*sin(delta);

%% Define E, F, G etc. parameters (for y-derivative equations, Table 8)

E = sin(delta)/R - ytilde*q/R^3;
F = dtilde/R^3 + Y32*sin(delta)*xi^2;
G = 2*X11*sin(delta) - ytilde*q*X32;
H = dtilde*q*X32 + xi*q*Y32*sin(delta);
P = cos(delta)/R^3 + q*Y32*sin(delta);
Q = 3*ctilde*dtilde/R^5 - sin(delta)*(z*Y32 + Z32 + Z0);

%% Define E', F', G' etc. parameters (for z-derivative equations, Table 9)

Ep = cos(delta)/R + dtilde*q/R^3;
Fp = ytilde/R^3 + Y32*cos(delta)*xi^2;
Gp = 2*X11*cos(delta) + dtilde*q*X32;
Hp = ytilde*q*X32 + xi*q*Y32*cos(delta);
Pp = sin(delta)/R^3 - q*Y32*cos(delta);
Qp = 3*ctilde*ytilde/R^5 + q*Y32 - cos(delta)*(z*Y32 + Z32 + Z0);

%% Create f, df/dx, df/dy, df/dz arrays depending on letter 'A', 'B', or 'C'

if strcmp(letter, 'A')
    
    f1 = -0.5*alpha2*ln_Reta - 0.5*alpha*Y11*q^2;
    f2 = -0.5*alpha2*ln_Rxi - 0.5*alpha*X11*q^2;
    f3 = 0.5*(theta - alpha*q*(eta*X11 + xi*Y11));
    
    df1x = -0.5*alpha2*xi*Y11 + 0.5*alpha*xi*Y32*q^2;
    df2x = -0.5*alpha2/R + 0.5*alpha*(q^2)/R^3;
    df3x = -0.5*alpha2*q*Y11 - 0.5*alpha*Y32*q^3;
    
    df1y = -0.5*alpha2*(cos(delta)/R + q*Y11*sin(delta)) - 0.5*alpha*q*F;
    df2y = -0.5*alpha2*ytilde*X11 - 0.5*alpha*q*G;
    df3y = 0.5*alpha2*(dtilde*X11 + xi*Y11*sin(delta)) + 0.5*alpha*q*H;
    
    df1z = 0.5*alpha2*(sin(delta)/R - q*Y11*cos(delta)) - 0.5*alpha*q*Fp;
    df2z = 0.5*alpha2*dtilde*X11 - 0.5*alpha*q*Gp;
    df3z = 0.5*alpha2*(ytilde*X11 + xi*Y11*cos(delta)) + 0.5*alpha*q*Hp;
    
elseif strcmp(letter, 'B')
    
    f1 = Y11*q^2 - (alpha2*I3*sin(delta)^2)/alpha;
    f2 = X11*q^2 + (alpha2*xi*sin(delta)^2)/alpha/Rdtilde;
    f3 = q*(eta*X11 + xi*Y11) - theta - (alpha2*I4*sin(delta)^2)/alpha;
    
    df1x = -xi*Y32*q^2 - alpha2*J4*(sin(delta)^2)/alpha;
    df2x = -(q^2)/R^3 - alpha2*J5*(sin(delta)^2)/alpha;
    df3x = Y32*q^3 - alpha2*J6*(sin(delta)^2)/alpha;
    
    df1y = q*F - alpha2*J1*(sin(delta)^2)/alpha;
    df2y = q*G - alpha2*J2*(sin(delta)^2)/alpha;
    df3y = -q*H - alpha2*J3*(sin(delta)^2)/alpha;
    
    df1z = q*Fp + alpha2*K3*(sin(delta)^2)/alpha;
    df2z = q*Gp + alpha2*xi*D11*(sin(delta)^2)/alpha;
    df3z = -q*Hp + alpha2*K4*(sin(delta)^2)/alpha;
    
elseif strcmp(letter, 'C')
    
    f1 = -alpha2*(sin(delta)/R + q*Y11*cos(delta)) - alpha*(z*Y11 - Z32*q^2);
    f2 = 2*alpha2*xi*Y11*sin(delta) + dtilde*X11 - alpha*ctilde*(X11 - X32*q^2);
    f3 = alpha2*(ytilde*X11 + xi*Y11*cos(delta)) + alpha*q*(ctilde*eta*X32 + xi*Z32);
    
    df1x = alpha2*xi*sin(delta)/R^3 + xi*q*Y32*cos(delta) + alpha*xi*(3*ctilde*eta/R^5 - 2*Z32 - Z0);
    df2x = 2*alpha2*Y0*sin(delta) - dtilde/R^3 + alpha*ctilde*(1 - 3*(q/R)^2)/R^3;
    df3x = -alpha2*(ytilde/R^3 - Y0*cos(delta)) - alpha*(3*ctilde*eta*q/R^5 - q*Z0);
    
    df1y = alpha2*(q/R^3 + Y0*sin(delta)*cos(delta)) + alpha*(z*cos(delta)/R^3 + 3*ctilde*dtilde*q/R^5 - q*Z0*sin(delta));
    df2y = -2*alpha2*xi*P*sin(delta) - ytilde*dtilde*X32 + alpha*ctilde*(X32*(ytilde + 2*q*sin(delta)) - X53*ytilde*q^2);
    df3y = -alpha2*(xi*P*cos(delta) - X11 + X32*ytilde^2) + alpha*ctilde*(X32*(dtilde + 2*q*cos(delta)) - X53*ytilde*eta*q) + alpha*xi*Q;
    
    df1z = -eta/R^3 + Y0*cos(delta)^2 - alpha*(z*sin(delta)/R^3 - 3*ctilde*ytilde*q/R^5 - Y0*sin(delta)^2 + q*Z0*cos(delta));
    df2z = 2*alpha2*xi*Pp*sin(delta) - X11 + X32*dtilde^2 - alpha*ctilde*(X32*(dtilde - 2*q*cos(delta)) - X53*dtilde*q^2);
    df3z = alpha2*(xi*Pp*cos(delta) + ytilde*dtilde*X32) + alpha*ctilde*(X32*(ytilde - 2*q*sin(delta)) + X53*dtilde*eta*q) + alpha*xi*Qp;
    
else
    error('please enter letter A, B, or C');   
end

f    = [f1, f2, f3];
dfdx = [df1x, df2x, df3x];
dfdy = [df1y, df2y, df3y];
dfdz = [df1z, df2z, df3z];

end