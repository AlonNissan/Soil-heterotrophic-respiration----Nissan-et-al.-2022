clc; clear all; close all;
%% Single water patch

% Parameters 
C0  = 10;           % Atmospheric Oxygen concentration [mol/m^3]
Cs  = 1;            % Dissolved organic carbon concentration [mol/m^3]
r0  = 5e-3;         % water patch radius [m]
Por = 0.3;          % porosity
lc  = 1e-5;         % characteristic grain size
Vm  = 1e-7;         % Michaelis–Menten maximum rate
kMo = 0.5;          % Michaelis constant for Oxygen   
kMs = 1;            % Michaelis constant for Dissolved organic carbon
Dm  = 1e-10;        % Dissolved oxygen diffusion coefficient

% non-dimensional parameters
alpha = kMo./C0;
beta  = (3.*Vm.*Por.*Cs.*r0.^2)./(Dm.*lc.*C0.*(Cs+kMs));

% PDE solver -- https://ch.mathworks.com/help/matlab/ref/pdepe.html
m       = 2;                % Symmetry constant [2] for 3-D Spherical coordinates with azimuthal and zenith angular symmetries
xmesh   = linspace(0,r0,500);
tspan   = linspace(1,1e5,10);
Rc      = (3*Por*Vm/lc)*(Cs/(Cs+kMs));


pde_an      = @(x,t,u,dudx) pdefun(x,t,u,dudx,Dm,Rc,kMo);
bc_an       = @(xl,ul,xr,ur,t) bcfun(xl,ul,xr,ur,t,C0);
sol         = pdepe(m,pde_an,@icfun,bc_an,xmesh,tspan);
sol_u       = sol(end,:); % taking the last solution, ~steady state
sol_reac    = ((3.*Vm.*Por./lc)*(Cs./(Cs+kMs)) .* (sol_u./(sol_u+kMo)));% converting from oxygen concentation within the patch to reaction ---> based on Michaelis–Menten





plot(xmesh,sol_u)
xlabel('Radial distance [m]');
ylabel('Oxygen concentation')


function [c,f,s] = pdefun(x,t,u,dudx,Dm,Rc,kMo)
c = 1/Dm;
f = dudx;
s = -(1/Dm)*(Rc)*(u/(u+kMo));
end


function u0 = icfun(x)
u0 = 0;
end


function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t,C0)
pl = 0;
ql = 1;
pr = ur-C0;
qr = 0;
end
