% This is the code for estimating flux distribution in steady-state. 
% Key: MOMA quadratic programming for flux distribution in metabolic network.
% Output: The final flux distribution including biomass predictions in netGS.
% This code is ruuning in Matlab. 
% Model details see below.
% Contact: tong@mpimp-golm.mpg.de
 
%% add path and toolbox
%addpath(genpath('/opt/MATLAB/tomlab'));
addpath(genpath('/opt/MATLAB/glpk'));
addpath(genpath('/pot/MATLAB/glpkmex'));
addpath(genpath('/opt/MATLAB/opencobra-cobratoolbox-7be8e9b'));
changeCobraSolver('glpk');
addpath('/../netGS/');
	
cd /../netGS/cvx/
cvx_setup

cd /../netGS/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % STEP 5 % % % % %
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

%%%%%% only FLUXES PREDICTABILITY > biomass INCLUDED %%%%%%%%%%

regr2 = regr2o((3*t+p-3),:);
aidi = find(regr2>=regr2(336));
Z = zeros(1,n);
Z(aidi) = 1; 
[aidim aidin] = size(aidi);

fluxa = fluxaiall(fidd,:);
[am an] = size(fluxa);

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

SS(:,n) = biofunc(:,i);

A = 1./fluxa(i,:);
AA = Z.*A;
AAn = AA(find(AA~=0));

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

% Output: final flux distribution including biomass predicition for netGS
fluxalli = sprintf('biomasspredict_r%d_f%d.csv',t,p);
csvwrite(fluxalli,fluxall)

aidinall = [aidinall;aidin];


end

end

% Output: number of fluxes included
csvwrite('numfluxbyr2.csv',aidinall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



