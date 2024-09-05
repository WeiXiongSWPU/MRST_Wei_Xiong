%% Introduction to Vapor-Liquid Equilibrium Calculations using eCPA
% The order of component is strict that is H2O, CO2, H2S, SO2, CH4, N2, O2, Ar
% Ref "Phase equilibrium modeling for carbondioxide Capture and Storage
% (CCS) fluids in brine using an electrolyte association equation of state"
mrstModule add compositional ad-core

% H2O-CO2-NaCl/KCl
P = 100.*barsa;
T = 273.15+45;
m = (0:6)';
x1 = 0.5./(2.*m./55.51+1);
x2 = m./55.51.*x1;
z = [x1,1-x1-2*x2,x2,x2];
names = {'Water','CarbonDioxide','Na+','Cl-'};
mixture = ECPATableCompositionalMixture(names);
eCPA = ECPAEquationOfStateModel([], mixture, 'eCPA');
[~, x, ~, Z_L, ~, rhoL] = eCPAstandaloneFlash(P, T, z, eCPA);
mu = eCPA.PropertyModel.computeViscosity(eCPA, P, x, Z_L, T, true);
mu = mu .* 1e3;

% brine
x1 = x(:,1)+x(:,2);
x3 = x(:,3);x4 = x(:,4);
z = [x1,1-x1-x3-x4,x3,x4];
[~, xp, ~, Z_Lp, ~, rhoLp] = eCPAstandaloneFlash(P, T, z, eCPA);
mup = eCPA.PropertyModel.computeViscosity(eCPA, P, xp, Z_Lp, T, true);
mup = mup .* 1e3;

msalt =55.51* x(:,3)./x(:,1);
res1 = [msalt,x(:,2), rhoL,rhoLp,mu,mup];

% H2O-CO2-CaCl2
x1 = 0.5./(3.*m./55.51+1);
x2 = m./55.51.*x1;
z = [x1,1-x1-3*x2,x2,2*x2];
names = {'Water','CarbonDioxide','Ca2+','Cl-'};
mixture = ECPATableCompositionalMixture(names);
eCPA = ECPAEquationOfStateModel([], mixture, 'eCPA');

[~, x, ~, Z_L, ~, rhoL] = eCPAstandaloneFlash(P, T, z, eCPA);
mu = eCPA.PropertyModel.computeViscosity(eCPA, P, x, Z_L, T, true);
mu = mu .* 1e3;

x1 = x(:,1)+x(:,2);
x3 = x(:,3);x4 = x(:,4);
z = [x1,1-x1-x3-x4,x3,x4];
[~, xp, ~, Z_Lp, ~, rhoLp] = eCPAstandaloneFlash(P, T, z, eCPA);
mup = eCPA.PropertyModel.computeViscosity(eCPA, P, xp, Z_Lp, T, true);
mup = mup .* 1e3;

msalt =55.51* x(:,3)./x(:,1);
res2 = [msalt,x(:,2), rhoL,rhoLp,mu,mup];