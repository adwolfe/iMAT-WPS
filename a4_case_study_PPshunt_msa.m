%% story of PPshunt

% * CAN BE REMOVED AFTER PUBLICATION * 
% not used in the ms

% block one of the two and redo FVA of the other 
% 'RM00705'	0.372765491570059	'Coenzyme A + Nicotinamide adenine dinucleotide + Malonate semialdehyde  -> Acetyl-CoA + CO2 + Nicotinamide adenine dinucleotide - reduced '	'coa[m] + nad[m] + msa[m]  -> accoa[m] + co2[m] + nadh[m] '	-0.372765491570059	2.57143570112192	3.11009820243296	0
% 'RM00706'	0.00833566554063898	'Coenzyme A + Nicotinamide adenine dinucleotide phosphate + Malonate semialdehyde  -> Acetyl-CoA + CO2 + Nicotinamide adenine dinucleotide phosphate - reduced '	'coa[m] + nadp[m] + msa[m]  -> accoa[m] + co2[m] + nadph[m] '	-0.00833566554063898	0.920940280719251	0.368684648673686	0

% CHANGED THE CONFIDENCE LEVEL ANNOTATION IF THEY ARE BOUNDED AFTER BLOCKING

% parameters 
yield = 0.65;
relCap_minLow = 0.05;
relCap_metFit = 0.05;

%% add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/% the gurobi path
addpath scripts/
addpath scripts/PlotPub/lib/
initCobraToolbox(false);
%% Load model
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');
load('./input/model/epsilon_generic_withUptakes.mat'); % see walkthrough_generic.m for guidance on generating the epsilon values
load('input/WPS/categ_expression_and_WPS.mat');
branchTbl = readtable('input/WPS/final_branchPoint_table.csv'); % must have 'mets', 'rxn1','rxn2', and 'maxCosine'

% setup the model
model = configurateModel(model);
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

% SPECIAL TREATMENT FOR ANALYZING MSA
% block the nadph alternative to study valid flux variability 
model = changeRxnBounds(model,'RM00706',0,'b'); 


% Set parameters
modelType = 2; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
speedMode = 1;

% run the program

% analyze the FVA with candidate responsiveness information removed

% run iMAT++ with yeild constraint

% set up the yield constraint
% we use the constraint (disassimilation) rate to constrain the bacteria waste
% this is to force the nutrient to be efficiently used instead of wasted in
% bulk
% add the disassimilation constraints 
model_coupled = model;
% add the disassimilation constraints 
model_coupled.S(end+1,:) = zeros(1,length(model_coupled.rxns));
model_coupled.S(end, strcmp('EXC0050',model_coupled.rxns)) = yield; 
model_coupled.S(end, strcmp('BIO0010',model_coupled.rxns)) = 1; 
model_coupled.csense(end+1) = 'G';
model_coupled.b(end+1) = 0;
model_coupled.BiGG(end+1) = {'NA'};
model_coupled.metCharges(end+1) = 0;
model_coupled.metFormulas(end+1) = {'NA'};
model_coupled.metNames(end+1) = {'Pseudo-metabolite that represents a constraint'};
model_coupled.mets(end+1) = {['NonMetConst',num2str(length(model_coupled.mets))]};

% start program
myCSM = struct(); % myCSM: my Context Specific Model
[myCSM.OFD,...
myCSM.PFD,...
myCSM.N_highFit,...
myCSM.N_zeroFit,...
myCSM.minLow,...
myCSM.minMetBalanceLoss,...
myCSM.minTotal,...
myCSM.minTotal_OFD,...
myCSM.MILP,...
myCSM.MILP_PFD,...
myCSM.HGenes,...
myCSM.RLNames,...
myCSM.OpenGene,...
myCSM.latentRxn,...
myCSM.Nfit_latent,...
myCSM.wasteDW,...
myCSM.RHNames,...
myCSM.branchMets,...
myCSM.minimizedLowRxns]...
= IMATplusplus_wiring_triple_inetgration_final(model_coupled,epsilon_f,epsilon_r, ExpCateg, branchTbl, modelType,speedMode,...
0.05, 0.05, 1, 1, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);

myCSM_triple = myCSM;

%% quick plain FVA
targetRxns = {'RM01608','RM00705','RM00908'};

% basic FVA parameters
parforFlag = 0;
relMipGapTol = 1e-12; % for fast result, use 1e-3
verbose = true;

% run FVA by calling:
[minval, maxval] = FVA_MILP(myCSM_triple.MILP_PFD, model, targetRxns,parforFlag,relMipGapTol,verbose);

minval
maxval
% ==> as expected, reaction RM00705 has to carry flux when its nadph
% alternative was blocked.

