%% calculate capillary pressure of BC, VG, and VG-IFT
mrstModule add ad-core ad-props mrst-gui compositional

% p, T, z
T = 84.4+273.15;
p = 300*barsa;
x = (0.5:0.01:0.99)';
z = [x,(1-x)];
nc = numel(x);
T = repmat(T, nc, 1);

% sW
ECPAmixture = ECPATableCompositionalMixture({'Water','Carbondioxide'});
eCPA = ECPAEquationOfStateModel([], ECPAmixture, 'eCPA');
[L1, x2, ~,Z_L, Z_V] = eCPAstandaloneFlash(p, T, z, eCPA);
sV = Z_V.*(1-L1)./(Z_V.*(1-L1)+Z_L.*L1);
sW = 1-sV;

% IFT
IFT_pure = eCPApureInterfacialTension(eCPA, T);
IFT = eCPAInterfacialTension(eCPA, T, x2, IFT_pure, L1);

% permeability, porosity
phi = 0.15;
k = 100*milli*darcy;
ratio = phi./k;
ratio = repmat(ratio, nc, 1);

% parameters of capillary pressure model
Pce = 0.64;
sWi = 0.2;
snt = 0.005;

% VG-IFT model
pc_VG_IFT = Pce .*((sW-sWi) ./ (1-sWi)).^(-0.5) .* IFT .* ratio.^0.5;
select = value(sW)>(1-snt);
if any(select)
    pc_VG_IFT(select) = Pce./snt .* ((1-snt-sWi)./(1-sWi)).^(-0.5).*(1-sW(select)).* IFT(select) .* ratio(select).^0.5;
end
pc_VG_IFT = 1e-5.*pc_VG_IFT;

% VG model
pc_VG = 0.2 .*((sW-sWi) ./ (1-sWi)).^(-0.5);
select = value(sW)>(1-snt);
if any(select)
    pc_VG(select) = 0.2./snt .* ((1-snt-sWi)./(1-sWi)).^(-0.5).*(1-sW(select));
end

% BC model
pc_BC = 0.2 .* ((sW-sWi) ./ (1-sWi)).^(-0.5);

% Launch interactive plotting
figure
plot(sW, pc_BC, 'or')
hold on
plot(sW, pc_VG, '-g')
hold on
plot(sW, pc_VG_IFT, '+b')
xlabel('sW')
ylabel('Capillary pressure (bar)');
legend('BC model','VG model','VG-IFT model')