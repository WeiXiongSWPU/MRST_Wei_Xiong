%% Benchmark for long-term CO2 sequestration simulations
% The aquifer has a dip of 1%, and it is deep enough such that the
% injected CO2 is in supercritical phase, which is immiscible with the
% resident water. The top and bottom are no-flow boundaries. The initial
% in-situ pressures at the left and right boundaries are held constant. 
mrstModule add ad-core ad-props mrst-gui compositional deckformat linearsolvers
gravity reset on

% The numerical model consists of 150 grid blocks in horizontal direction 
% and 50 grid blocks in vertical direction.The grid blocks in horizontal 
% direction have the same size of 100 meter. The grid blocks in the
% vertical direction have the same size of 1.25 meter.
nx = 0:100:15000;
ny = 0:100:100;
nz = 0:1.25:50;
dz = [150:-1:0,150:-1:0];

G = tensorGrid(nx, ny, nz,'depthz', dz);
G = computeGeometry(G);
plotGrid(G); view(3); axis tight

% Assume an isotropic and homogenous aquifer with a permeability of 100 mD,
% porosity of 0.15
rock = makeRock(G, 100*milli*darcy, 0.15);
% diffusion
rock.Mechanisms.diffusion = 1;
rock.Di=[0,2,0,0]*10^-9;
rock.tau = 2;

% Capillary pressure model
rock.VG = struct('Pce', 0.2*barsa, 'sWi', 0.2, 'snt', 0.005);

f = initSimpleADIFluid('phases', 'wg', 'blackoil', false, 'rho', [1000, 700], 'n', [4, 2]);
f.pcWG = 1;
ECPAmixture = ECPATableCompositionalMixture({'Water','Carbondioxide', 'K+', 'Cl-'});

% Construct model
ECPAarg = {G, rock, f, ...                              % Standard arguments
       ECPAmixture,...                                  % Compositional mixture
       'water', true, 'oil', false, 'gas', true,...     % Water-Gas system
       'liquidPhase', 'W', 'vaporPhase', 'G'};          % Water=liquid, gas=vapor
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'rowMajor', true);
ECPAoverall = ECPAGenericOverallCompositionModel(ECPAarg{:}, 'AutoDiffBackend', diagonal_backend);   % Overall mole fractions
ECPAoverall = imposeRelpermScaling(ECPAoverall, 'SWCR', 0.2, 'KRG', 0.4,'SGU',0.8);

% Validate model
ECPAoverall = ECPAoverall.validateModel();
for ii = 1:numel(rock.Di)
    Diff = makeRock(G, rock.Di(ii), NaN);
    T1 = getFaceTransmissibility(G, Diff);
    ECPAoverall.operators.T_diff{ii} = T1(ECPAoverall.operators.internalConn);
end
%% Set up BC + initial conditions/initial guess
p0 = 300*barsa; T = 84.4+273.15;

% salinity
z = [0.784, 0, 0.108,0.108];

eCPA = ECPAEquationOfStateModel([], ECPAmixture, 'eCPA');
[~, ~, ~,~, ~,rho0] = eCPAstandaloneFlash(p0, T, z, eCPA);  
g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* rho0, [z_0, z_max], p0);
p = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil

% The simulations are continued for up to 4000 years after CO2 injection has
% stopped.
totTime = 4000*year;

% constant pressure boundary conditions
bc = pside([], G, 'xmax', p(G.cells.centroids(:,1)==14950), 'sat', [1, 0]);       % Standard bc
bc = pside(bc, G, 'xmin', p(G.cells.centroids(:,1)==50), 'sat', [1, 0]);       % Standard bc
bc.components = repmat(z, numel(bc.face), 1);      % Boundary z

% CO2 is injected at a rate of 9000 t/year for 20 years
[~, ~, ~,~, ~,rho] = eCPAstandaloneFlash(1*barsa, 298.15, [0, 1, 0, 0], eCPA);  
W = addWell([], G, rock, 5900, 'name', 'inj', 'type', 'rate', 'Compi', ...
[0 1], 'val', 9000*1000/rho/year, 'components', [0, 1, 0, 0]);

% Plot well
show = true([G.cells.num, 1]);
cellInx = 50:150:5900;
show(cellInx) = false;
plotCellData(G, convertTo(p, barsa), show, 'EdgeColor', 'k')
plotWell(G, W, 'height', 10)
view(0, 0), camproj perspective

% Initialize the problem
s0 = [1, 0];
state0 = initResSol(G, p, s0);
state0.T = repmat(T, G.cells.num, 1);
state0.components = repmat(z, G.cells.num, 1);

% CO2 is injected at a rate of 9000 t/year for 20 years
dt1 = rampupTimesteps(20*year, 0.5*year, 12);
dt2 = rampupTimesteps(totTime-20*year, 1*year, 20);
dt=[dt1;dt2];

schedule = simpleSchedule(dt,'W', W, 'bc', bc);
schedule.step.control(numel(dt1)+1:end) = 2;
schedule.control(2) = schedule.control(1);
schedule.control(2).W.val = 0; 

% Set up nonlinear solver with high report level to output intermediate
% states
nls = NonLinearSolver('useRelaxation', true);
[~, ECPAstateso, ECPAreporto] = simulateScheduleAD(state0, ECPAoverall, schedule, 'nonlinearsolver', nls);

%% Launch interactive plotting
figure; 
plotToolbar(G, ECPAstateso)
view(0,0)
colorbar

PV=[];
time = cumsum(dt)./year;
for i = 1:numel(ECPAstateso)
    sg = ECPAstateso{i}.s(:,2);
    vg = 12500*0.15.*sg;
    PV=[PV; time(i),sum(vg)];
end
figure
plot(PV(:,1),PV(:,2));

Tip=[];
for i = 1:numel(ECPAstateso)
    sg = ECPAstateso{i}.s(1:150,2);
    vg = 100*find(sg>0, 1, 'last' );
    if isempty(vg)
       vg=5000;
    end
    Tip=[Tip; time(i),vg-5000];
end
figure
plot(Tip(:,1),Tip(:,2));

fr=[];
for i = 1:numel(ECPAstateso)
    L =ECPAstateso{i}.L;
    x =ECPAstateso{i}.x(:,2);
    y =ECPAstateso{i}.y(:,2);
    zco2=x.*L+y.*(1-L);
    n=x.*L;
    X = sum(n)./sum(zco2);
    fr = [fr;X];
end
figure
plot(cumsum(dt./year), fr)
fra=[cumsum(dt./year), fr];

f=[];
pv = poreVolume(G, rock);
mw=[18.015, 44.01,39,35.45].*1e-3;
for i = 1:numel(ECPAstateso)
    p =ECPAstateso{i}.pressure;
    y =ECPAstateso{i}.y;
    x =ECPAstateso{i}.x;
    Z_V =ECPAstateso{i}.Z_V;
    Z_L =ECPAstateso{i}.Z_L;
    s =ECPAstateso{i}.s;
    rhg = eCPA.PropertyModel.computeDensity(eCPA, p, y, Z_V, T, false);
    massg = pv.*rhg.*s(:,2);
    mwg = y(:,2).*mw(2)./sum(y.*mw, 2);
    massg = sum(massg .* mwg);

    rhl = eCPA.PropertyModel.computeDensity(eCPA, p, x, Z_L, T, true);
    massl = pv.*rhl.*s(:,1);
    mwl = x(:,2).*mw(2)./sum(x.*mw, 2);
    massl = sum(massl .* mwl);

    f = [f;massg,massl];
end
figure
plot(cumsum(dt./year), f)
mass=[cumsum(dt./year), f];
hold on
plot(cumsum(dt./year), sum(f, 2))

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