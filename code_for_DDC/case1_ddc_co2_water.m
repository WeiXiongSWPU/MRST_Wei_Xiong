%% Implementation of Dissolution-Diffusion-Convection
% Implementation of Fickian Diffusion in ECPAGenericOverallCompositionModel 
% and ECPAGenericNaturalVariablesModel
mrstModule add compositional ad-core ad-props mrst-gui
gravity reset on

% Grid and petrophysical data
nx = 0:0.01:1;
ny = 0:1:1;
nz = (0:0.001:0.01);
cellsize = max(nz);
while cellsize<1
    cellsize = cellsize + 0.01;
    nz = [nz,cellsize];
end

G = tensorGrid(nx, ny, nz);
G = computeGeometry(G);

% Generate a random Gaussian permeability field
n_x = numel(nx)-1; n_z = numel(nz)-1;
rng(0);
K = 10*ones(n_x*n_z,1);
K = K + 0.01.*(2.*rand(n_x*n_z,1)-1).*K;
rock = makeRock(G, K*darcy, 0.3);

figure; 
plotCellData(G, K);
colorbar;
view(0,0)

% diffusion
rock.Mechanisms.diffusion = 1;
rock.Di=[0,2]*10^-9;
rock.tau = 1;

% Start with fluid properties from a standard black-oil, water-gas model
f = initSimpleADIFluid('phases', 'wg', 'blackoil', false, 'rho', [1000, 700]);
mixture = ECPATableCompositionalMixture({'Water','Carbondioxide'});

% Construct models for both formulations. Same input arguments
arg = {G, rock, f, ...                              % Standard arguments
       mixture,...                                  % Compositional mixture
       'water', true, 'oil', false, 'gas', true,... % Water-Gas system
       'liquidPhase', 'W', 'vaporPhase', 'G'};      % Water=liquid, gas=vapor
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'rowMajor', true);
overall = ECPAGenericOverallCompositionModel(arg{:}, 'AutoDiffBackend', diagonal_backend);   % Overall mole fractions
natural = ECPAGenericNaturalVariablesModel(arg{:}, 'AutoDiffBackend', diagonal_backend);     % Natural variables

% Validate both models to initialize the necessary state function groups
overall = overall.validateModel();
natural = natural.validateModel();
for ii = 1:numel(rock.Di)
    Diff = makeRock(G, rock.Di(ii), NaN);
    T1 = getFaceTransmissibility(G, Diff);
    natural.operators.T_diff{ii} = T1(natural.operators.internalConn);
end
natural = setWENODiscretization(natural);

p0 = 100*barsa; T = 273.15 + 45;
z1 = [1-0.022,0.022]; z2 = [1,0];
eCPA = ECPAEquationOfStateModel([], mixture, 'eCPA');

[L1, x1, ~,Z_L, Z_V,rho] = eCPAstandaloneFlash(p0, T, z2, eCPA);  
s2 = Z_V*(1-L1)/(Z_V*(1-L1)+Z_L*L1);
mu0 = eCPA.PropertyModel.computeViscosity(eCPA, p0, x1, Z_L, T, true);

bc = pside([], G, 'Top', p0, 'sat', [1, 0]);            % Standard bc
bc.components = repmat(z1, numel(bc.face), 1);          % Boundary z

g = norm(gravity);
p = zeros(n_x, n_z);
p(:,1)=p0 + rho.*(nz(2)-nz(1)).*g./2;
[~, ~, ~,~,~,rho] = eCPAstandaloneFlash(p(:,1), T, z2, eCPA);
rho(2)=rho(1);
for i = 1:n_z-1
    for j = 1:2000
        p(:,1+i) = p(:,i) + (rho(i).*(nz(i+1)-nz(i))+rho(i+1).*(nz(i+2)-nz(i+1))).*g./2;
        [L2, ~, ~,~,~,rho2] = eCPAstandaloneFlash(p(1,1+i), T, z2, eCPA);
        if abs(rho2-rho(i+1)) < 1e-12
            rho(i+2)=rho2;
            break
        else
            rho(i+1)=rho2;
        end
    end
end

p=p(:);
z = repmat(z2,n_z*n_x,1);
s0 = repmat([1, 0],n_x*n_z,1);

state0 = initResSol(G, p, s0);
state0.T = repmat(T, G.cells.num, 1);
state0.components = z;

dt = rampupTimesteps(10*day, 0.1*day, 10);

schedule = simpleSchedule(dt, 'bc', bc);
nls = NonLinearSolver('useRelaxation', true);

[ws1, states1, reports1] = simulateScheduleAD(state0, natural, schedule, 'nonlinearsolver', nls);
% [ws2, states2, reports2] = simulateScheduleAD(state0, overall, schedule, 'nonlinearsolver', nls);

figure, 
plotToolbar(G, states1)
view(0,0);colorbar

pe = states1{end}.pressure;
ze = states1{end}.Z_L;
ve = ze.*8.314.*T./pe;
xe = states1{end}.x;
me = xe(:,1).*0.018+xe(:,2).*0.044;
re = me./ve;
figure
plotToolbar(G, re)
view(0,0);colorbar

diffusionFlux=[];
time = cumsum(dt);
for i = 1:numel(states1)
     Flux = abs(sum(states1{i}.diffusionFlux(:,2)));
    diffusionFlux = [diffusionFlux; time(i), Flux];
end
figure
plot(diffusionFlux(:,1), diffusionFlux(:,2))


%%
%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}