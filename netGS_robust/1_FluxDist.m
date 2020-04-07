% This is the code for robustness test of the reference flux distribution. 
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
addpath('/../netGS_robust/');

cd /../netGS_robust/cvx/
cvx_setup

cd /../netGS_robust/

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

% the max with col0 biomass function
c = model.c.';
bioc = fluxc(find(c==1),:);
%bioc = 0.003124746706309;
					
idnzero = union(find(fluxc>1E-5),find(fluxc<-1E-5)); %% non-zero flux id
idzero = intersect(find(fluxc<1E-5),find(fluxc>-1E-5)); %% zero flux id

csvwrite('nonzeroid.csv',idnzero);
%csvwrite('zeroid.csv',idzero);

% Output: Col0 flux distribution
save('fluxcol0.mat','fluxc')
csvwrite('fluxcol0.csv',fluxc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % STEP 2 % % % % % 
% % % % % sampled Col0 flux distribution % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fluxc = load('fluxcol0.mat','fluxc');
fluxc = fluxc.fluxc;

idnzero = csvread('nonzeroid.csv',0,0);
fluxc = fluxc(idnzero);

fluxrrall = [];
for i=1:336,
	
	fluxi = fluxc(i,1);
	
	%%%% sampling respect to the error of 20%
	a = fluxi*0.8;
	b = fluxi*1.2;

	%%% run for 50 random samples 
	rr = (b-a).*rand(50,1) + a;
	fluxrrall = [fluxrrall;rr.'];

end

% Output: random sampled Col0 flux distribution
save('fluxcol0_sample.mat','fluxrrall')
csvwrite('fluxcol0_sample.csv',fluxrrall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % STEP 3 % % % % % 
% % % % % sampled Col0 flux distribution in steady-state as robustness % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% nonzero flux only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aramodel = readCbModel('model.xml');

fluxcr = load('fluxcol0_sample.mat','fluxrrall');
fluxcr = fluxcr.fluxrrall;

S = full(aramodel.S);
[Srow Scol] = size(S);
c = aramodel.c.';
n = Scol;

idnzero = csvread('nonzeroid.csv',0,0);

[n m] = size(fluxcr);

S = S(:,idnzero);
[Srow Scol] = size(S);

c = c(idnzero);

%%% fixed biomass value
biomsall = round(bioc,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% quadratic programming by cvx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b = 1;

carb = 6;
oxy = 61;
starch = 20;
suc = 31;

ccou = repelem(0,n);
ccou(carb) = 1;
ccou(oxy) = -2.8764;

ccol = repelem(0,n);
ccol(carb) = 1;
ccol(oxy) = -2.8762;

cssu = repelem(0,n);
cssu(starch) = 1;
cssu(suc) = -2.5797;

cssl = repelem(0,n);
cssl(starch) = 1;
cssl(suc) = -2.5795;

lb = [aramodel.lb(idnzero)];
ub = [aramodel.ub(idnzero)];


fluxall = [];

SS = S;

for rr=1:50,

A = 1./fluxcr(:,rr).';


x = [];

cvx_begin quiet

	variable x(n);
	minimize(norm(A*x-b));

	subject to

	SS*x == repelem(0,Srow).';	
	lb <= x <= ub;
	
	c*x == biomsall;

	ccou*x <= 0;
	ccol*x >= 0;
	cssu*x <= 0;
	cssl*x >= 0;

cvx_end


fluxall = [fluxall,x];

end

% Output: random sampled Col0 flux distribution in steady-state as robustness
save('fluxcol0_robust.mat','fluxall')
csvwrite('fluxcol0_robust.csv',fluxall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


