% This is the code for estimating the reference and other genotype flux distributions. 
% Key: MOMA quadratic programming for flux distribution in metabolic network.
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
% % % % % STEP 1 % % % % % 
% % % % % Col0 flux distribution % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read model
model = readCbModel('model.xml');

% ID for carboxylation, oxygenation, starch and sucrose reactions
carb = 6;
oxy = 85;
starch = 22;
suc = 43;

% set the ratio according to meansurements
co_ratio = 10.7/3.72;
ss_ratio = 6.32/2.45;

model.S(end+1,[carb oxy]) = [1 -co_ratio];
model.b(end+1) = 0;
model.mets(end+1) = {'carb/oxy ratio'};
model.metNames(end+1) = {'carb/oxy ratio'};

model.S(end+1,[starch suc]) = [1 -ss_ratio];
model.b(end+1) = 0;
model.mets(end+1) = {'starch/suc ratio'};
model.metNames(end+1) = {'starch/suc ratio'};

% run FBA
sol = optimizeCbModel(model);
v=sol.x;

%%%%%% col0 flux %%%%%%
fluxc = v; 

idnzero = union(find(fluxc>1E-5),find(fluxc<-1E-5)); %% non-zero flux id
idzero = intersect(find(fluxc<1E-5),find(fluxc>-1E-5)); %% zero flux id

csvwrite('nonzeroid.csv',idnzero);
%csvwrite('zeroid.csv',idzero);

% Output: Col0 flux distribution
save('fluxcol0.mat','fluxc')
csvwrite('fluxcol0.csv',fluxc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % STEP 2 % % % % % 
% % % % % maximum biomass under constraints used to scale measured biomass % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% accession specific biomass function
%%%% relaxed CO and SS ratio 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read biomass reaction specific values
biofunc = csvread('biomreaction.csv',0,0);

zmaxall=[];

for i=1:67,

% read model
model = readCbModel('model.xml');

% replace biomass reaction to specific 
model.S(:,549) = biofunc(:,i);

S = full(model.S);
[Srow Scol] = size(S);
c = model.c.';
n = Scol;

f = c;

c1 = S;

carb = 6;
oxy = 85;
starch = 22;
suc = 43;

ox = repelem(0,n);
ox(carb) = 1;
c3 = ox;

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

c41 = ccou;
c42 = ccol;

c51 = cssu;
c52 = cssl;

Aeq = vertcat(c1,c3,c41,c42,c51,c52);

r1 = repelem(0,Srow);
beq = [r1,0,0,0,0,0];

ctype = [repelem('S',Srow),'L','U','L','U','L'];
vartype = repelem('C',n);
sense = -1;

lb = [model.lb.'];
ub = [model.ub.'];

%% max biomass in model with CO and SS ratio
[xmin, fmin, status, extra] = glpk (f, Aeq, beq, lb, ub, ctype, vartype, sense);

zmaxall = [zmaxall,fmin];
 
end

%%% average overall accessions 
bb = mean(zmaxall);
zmax = bb;
%zmax = 0.003503620892060; 
%%% zmax as a scale to ensure the feasibility in minimization program later. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % STEP 3 % % % % % 
% % % % % genotype flux distribution by constraining: % % % % % 
% % % % % 1.the biomass ratio % % % % % 
% % % % % 2.specific biomass functions % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% nonzero flux only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aramodel = readCbModel('model.xml');

fluxc = load('fluxcol0.mat','fluxc');
fluxc = fluxc.fluxc;

S = full(aramodel.S);
[Srow Scol] = size(S);
c = aramodel.c.';
n = Scol;

idnzero = csvread('nonzeroid.csv',0,0);

fluxc = fluxc(idnzero);
[n m] = size(fluxc);

S = S(:,idnzero);
[Srow Scol] = size(S);

c = c(idnzero);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scale biomass value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

biom = csvread('biomass_optN.csv',0,1);
biom = biom(:,1);

biomc_m = biom(15,:); %% Col0 in measurement
biomc_p = fluxc(find(c==1),:); %% Col0 in model

biomratio = biom*biomc_p/biomc_m;

biommax = max(biomratio); % max biomass in measurement

delta = 0.00011; %% add small error 
biomsall = (zmax-delta)*biomratio/biommax; %% scale biomass by maximum

biomsall = round(biomsall,5);

SEM = std(biomsall)/sqrt(length(biomsall));           
ts = [-1.645,1.645]; %%% 90% bimass interval
%CI = mean(biomsall)+ts*SEM;                      
ee = ts*SEM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% quadratic programming by cvx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = 1./fluxc.';
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

biofunc = csvread('biomreaction.csv',0,0);

fluxall = [];
optvalall = [];
SS = S;

for i=1:67,

SS(:,n) = biofunc(:,i);

x = [];

cvx_begin quiet

	variable x(n);
	minimize(norm(A*x-b));

	subject to

	SS*x == repelem(0,Srow).';	
	lb <= x <= ub;
	
	biomsall(i)+ee(1) <= c*x <= biomsall(i)+ee(2);

	ccou*x <= 0;
	ccol*x >= 0;
	cssu*x <= 0;
	cssl*x >= 0;

cvx_end


fluxall = [fluxall,x];
optvalall = [optvalall,cvx_optval];

end

% Output: genotype flux distribution
csvwrite('fluxgenotype.csv',fluxall);
save('fluxgenotype.mat','fluxall');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


