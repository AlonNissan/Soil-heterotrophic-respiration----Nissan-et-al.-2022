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


% import comsol libaries and run
import com.comsol.model.*
import com.comsol.model.util.*
mphopen('RDE_basedModeling_SS.mph');
model.param.set('alpha', strcat(num2str(alpha)));
model.param.set('beta', strcat(num2str(beta)));
ModelUtil.showProgress(true);
model.sol('sol1').runAll;




% plotting oxygen concentation inside the water patch
pd = mphplot(model,'pg1','rangenum',1);
title('Oxygen normalized concentation')
cb = colorbar;
cb.FontSize = 42;
x=get(cb,'Position');
set(cb,'Position',[0.75 0.1 0.025 0.8],'FontSize',52)
pbaspect([1 1 1])
