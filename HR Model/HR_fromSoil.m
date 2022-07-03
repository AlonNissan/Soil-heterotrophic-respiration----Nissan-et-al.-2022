clc; clear all; close all;
%% Respiration from the soil (bulk)

% Parameters 
Theta   = 0.4;                  % soil saturation degree (relative,0 -- 1)
L       = 1;                    % characteristic dimension of the soil [m]
C0      = 10;                   % Atmospheric Oxygen concentration [mol/m^3]
Cs      = 1;                    % Dissolved organic carbon concentration [mol/m^3]
Por     = 0.3;                  % porosity
lc      = 1e-4;                 % characteristic grain size [m]
Vm      = 1e-10;                 % Michaelis–Menten maximum rate [mol/(m^3*s)]
kMo     = 0.1;                  % Michaelis constant for Oxygen  [mol/(m^3)] 
kMs     = 1;                    % Michaelis constant for Dissolved organic carbon [mol/(m^3)]
Dm      = 1e-9;                 % Dissolved oxygen diffusion coefficient [m^2/s]
gamma1  = 100;                  % see equation 3 in the manuscript
tau     = 0.218;                % see equation 3 in the manuscript
xi      = gamma1.*(1-Theta);    % see equation 3 in the manuscript


NP = ((L/lc).^3).^(1-Theta).*(((1./xi).^(tau-1)).*gamma(1-tau) - expint(1./xi)); % number of water patches per L^3 soil
SP = (Theta./NP);                                                                % characteristic size of water patches
r0 = ((3.*SP./1)./(pi.*4)).^(1./3);    % characteristic partch's radius [m]                                              % characteristic radius of the water patches




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


Rv      = linspace(xmesh(1),xmesh(end),length(xmesh)).^3;
volR    = (4/3).*pi.* diff(Rv);
Hr_wp   = mean([sol_reac(1:end-1);sol_reac(2:end)]);
solSum  = sum(Hr_wp.*volR);
fac     = (365*24*60*60);  % convert from sec to year
HR      = solSum.*NP.*fac; % [mol/(m^2*year)]  

str = sprintf('Respiration rate is: %f [mol/(m^2*year)]',HR)



% plot(xmesh,sol_u)
% xlabel('Radial distance [m]');
% ylabel('Recation')

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




