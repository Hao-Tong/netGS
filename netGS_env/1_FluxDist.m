% This is the code for estimating the reference flux distribution in optimal (E1) and low N (E2) conditions. 
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
addpath('/../netGS_env/');

cd /../netGS_env/cvx/
cvx_setup 

cd /../netGS_env/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % STEP 1 % % % % % 
% % % % % Col0 flux distribution in E1 % % % % % 
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

fluxc = fluxc(idnzero);

% Output: Col0 flux distribution for optimal N
save('fluxcol0_optN.mat','fluxc')
csvwrite('fluxcol0_optN.csv',fluxc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % STEP 2 % % % % % 
% % % % % Col0 flux distribution in E2 using 1. low N biomass function % % % % % 
% % % % % % % % % % % % % % % % % % % % % %  2. the biomass ratio % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% non-zero flux only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aramodel = readCbModel('model.xml');

S = full(aramodel.S);
[Srow Scol] = size(S);
c = aramodel.c.';
n = Scol;

%%% this the low N biomass function
SlN = S(:,548);

idnzero = csvread('nonzeroid.csv',0,0);

fluxc = load('fluxcol0_optN.mat','fluxc');
fluxc = fluxc.fluxc;

[n m] = size(fluxc);

S = S(:,idnzero);
[Srow Scol] = size(S);

c = c(idnzero);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scale biomass value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

biom = csvread('biomass_optNlowN.csv');
biom1 = biom(:,2); %% high N % E1
biom2 = biom(:,3); %% low N % E2

biomc_m = biom1(15,:); %% Col0 in measurement E1
biomc_p = fluxc(find(c==1),:); %% Col0 in model E1

biomc_me = biom2(15,:); %% Col0 in measurement E2

biomratio = biomc_me*biomc_p/biomc_m;
biomratio = round(biomratio,5); %% Col0 in model E1

biomsall = biomratio;

%SEM = std(biomratio)/sqrt(length(biomratio));           
%ts = [-2.576,2.576]; %%% 99% bimass interval
%CI = mean(biomsall)+ts*SEM;                      
ee = [-0.00001,0.00001];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% quadratic programming by cvx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

lb = aramodel.lb(idnzero);
ub = aramodel.ub(idnzero);

%%% replace the biomass function to low N biomass function 
S(:,336) = SlN;

x = []; 

cvx_begin quiet

	variable x(n);
	minimize(norm(A*x-b));

	subject to

	S*x == repelem(0,Srow).';	
	lb <= x <= ub;

	biomsall+ee(1) <= c*x <= biomsall+ee(2);

	ccou*x <= 0;
	ccol*x >= 0;
	cssu*x <= 0;
	cssl*x >= 0;

cvx_end

fluxc2 = x;

% Output: Col0 flux distribution for low N
csvwrite('fluxcol0_lowN.csv',fluxc2);
save('fluxcol0_lowN.mat','fluxc2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


