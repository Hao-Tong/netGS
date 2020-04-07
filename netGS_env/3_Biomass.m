% This is the code for estimating flux distribution in steady-state adding the exchange reaction ratios. 
% Key: MOMA quadratic programming for flux distribution in metabolic network.
% Output: The final flux distribution including biomass predictions in netGS across environment.
% This code is ruuning in Matlab. 
% Model details see below.
% Contact: tong@mpimp-golm.mpg.de
 
%% add path and toolbox
%addpath(genpath('/opt/MATLAB/tomlab'));
addpath(genpath('/opt/MATLAB/glpk'));
addpath(genpath('/pot/MATLAB/glpkmex'));
addpath(genpath('/opt/MATLAB/opencobra-cobratoolbox-7be8e9b'));
changeCobraSolver('glpk');
addpath('/../netGS_env/');

cd /../netGS_env/cvx/
cvx_setup

cd /../netGS_env/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % STEP 3 % % % % %
% % % % % genotype flux distribution in steady-state % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% nonzero flux only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aramodel = readCbModel('model.xml');

S = full(aramodel.S);
c = aramodel.c.';

idnzero = csvread('nonzeroid.csv',0,0);

S = S(:,idnzero);
[Srow Scol] = size(S);
n = Scol;

c = c(idnzero);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% flux id 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% read read cross-validation fold id
fido = csvread('foldid.csv');

% read flux prediction accuracy
regr2o = csvread('fluxpredict_cor.csv');
regido = 1:336;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% exchange reaction ratios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Col0 flux in optimal N condition
fluxc = load('fluxcol0_optN.mat','fluxc');
fluxc = fluxc.fluxc;

% Col0 flux in low N condition
fluxc2 = load('fluxcol0_lowN.mat','fluxc2');
fluxc2 = fluxc2.fluxc2;

%% ID of exchange reactions
exid = 328:333;
exfluxc = fluxc(exid);
exfluxc2 = fluxc2(exid);

%% the ratio of exchange fluxes between E1 abnd E2
exratio = exfluxc2./exfluxc;
exratio =round(exratio,5);

ee = [-0.1,0.1];
exratio1 = exratio+ee(1);
exratio2 = exratio+ee(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% in cross-validation sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aidinall = [];

%% replications
for t = 1:50,

t

fid = fido(:,t);

fluxai = sprintf('fluxpredict_%d.csv', t);
fluxai = importdata(fluxai);
[fm fn] = size(fluxai);
fluxaiall = zeros(fm,n);
fluxaiall(:,regido) = fluxai;

%% three-fold cross-validation
for p = 1:3,

fidd = find(fid==p);

%%%%%% FLUXES PREDICTABILITY > biomass INCLUDED %%%%%%%%%%

regr2 = regr2o((3*t+p-3),:);
aidi = find(regr2>=regr2(336));
Z = zeros(1,n);
Z(aidi) = 1; 
[aidim aidin] = size(aidi);

fluxa = fluxaiall(fidd,:);
[am an] = size(fluxa);

exfluxx = fluxa(:,exid).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% quadratic programming by cvx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b = 1; 

carb = 6;
oxy = 61;
starch = 20;
suc = 31;

ccou = repelem(0,n);
ccou(carb) = 1;
ccou(oxy) = -3.81445;

ccol = repelem(0,n);
ccol(carb) = 1;
ccol(oxy) = -0.93815;

cssu = repelem(0,n);
cssu(starch) = 1;
cssu(suc) = -3.3694;

cssl = repelem(0,n);
cssl(starch) = 1;
cssl(suc) = -0.7898;

lb = [aramodel.lb(idnzero)];
ub = [aramodel.ub(idnzero)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

biofunco = csvread('biomreaction.csv',0,0);
biofunc = biofunco(:,fidd);


fluxall = [];
SS = S;

for i= 1:am,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb = [aramodel.lb(idnzero)];
ub = [aramodel.ub(idnzero)];

%%%% Here is the genotype-specific biomass function measured in E1!!!
SS(:,n) = biofunc(:,i);

A = 1./fluxa(i,:);
AA = Z.*A;
AAn = AA(find(AA~=0));

exfluxxr1 = exfluxx(:,i).*exratio1;
exfluxxr2 = exfluxx(:,i).*exratio2;

lb(exid([1 2 3 5])) = exfluxxr1([1 2 3 5]);
ub(exid([1 2 3 5])) = exfluxxr2([1 2 3 5]);

ub(exid([1])) = 1000;

x = [];

cvx_begin quiet

	variable x(n);
	minimize(norm(AAn*x(aidi)-b));

	subject to

	SS*x == repelem(0,Srow).';	
	lb <= x <= ub;

	ccou*x <= 0;
	ccol*x >= 0;
	cssu*x <= 0;
	cssl*x >= 0;

	c*x >= 0;

cvx_end


fluxall = [fluxall,x];

end

% Output: final flux distribution including biomass predicition for netGS across environment
fluxalli = sprintf('biomasspredict_lowN_r%d_f%d.csv',t,p);
csvwrite(fluxalli,fluxall)


end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



