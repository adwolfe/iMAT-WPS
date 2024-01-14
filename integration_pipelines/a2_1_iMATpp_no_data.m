function a2_1_iMATpp_no_data(yield)
% Load model - this is the same across all five integrations 
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');
load('./input/model/epsilon_generic_withUptakes.mat'); % see walkthrough_generic.m for guidance on generating the epsilon values

% setup the model
model = configurateModel(model);
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

% Set parameters
modelType = 2; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
speedMode = 1;

% run iMAT++ with yeild constraint

% make an empty gene category (all as dynamic)
ExpCateg = struct();
ExpCateg.zero = {};
ExpCateg.low = {};
ExpCateg.dynamic = model.genes;
ExpCateg.high = {};
ExpCateg.responsive = {};
ExpCateg.nonresponsive = {};

% additionally, to setup the no-integration control, we need to impose a
% biomass constraint to make it an objective-based prediction
% we set the minimal biomass production arbituraryly at 0.8, which is
% roughly the biomass production in all other integrations
model.lb(strcmp(model.rxns,'BIO0010')) = 0.8;

% set up the yield constraint
% we use the constraint (disassimilation) rate to constrain the bacteria waste
% this is to force the nutrient to be efficiently used instead of wasted in
% bulk
% add the disassimilation constraints 
model_coupled = model;
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


% iMAT++
myCSM = struct(); % myCSM: my Context Specific Model
[myCSM.OFD,...
myCSM.PFD,...
myCSM.N_highFit,...
myCSM.N_zeroFit,...
myCSM.minLow,...
myCSM.minTotal,...
myCSM.minTotal_OFD,...
myCSM.MILP,...
myCSM.MILP_PFD,...
myCSM.HGenes,...
myCSM.RLNames,...
myCSM.OpenGene,...
myCSM.latentRxn,...
myCSM.Nfit_latent,...
myCSM.wasteDW]...
= IMATplusplus_wiring_dual_integration_final(model_coupled,epsilon_f,epsilon_r, ExpCateg, modelType,speedMode,...
1e-5, 1, 1, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10); % IMATplusplus_wiring_dual_integration_final is generally similar to original iMAT++ since it doesnt have the metabolite coupling integration step

myCSM_no_data = myCSM;

% NOTE: this flux distribution is equvelent to that obtained by directly 
% minimization of total flux on model_coupled (ATPm needs to be changed 
% to 10); we have confirmed this identity. We used the exact same 
% integration function to be consistent in the way that the simulation is
% set up but it is worth noting that this produced distribution is
% identical to a pFBA distribution under the same constraints (yield,
% biomass production and ATPm). 

save('output/integration_output/myCSM_no_data.mat','myCSM_no_data');


%% codes for manual interactions
% %% fast FVA
% addpath scripts/parfor_wait/
% targetRxns = model.rxns;
% parforFlag = 1;
% relMipGapTol = 1e-12; % this is to be redone with maximum precision, otherwsie the box will have problem
% verbose = true;
% 
% % we use the base model as a control, therefore, the FVA constraints were
% % given by the two reference model
% load('myCSM_merged.mat')
% load('myCSM_expression_only.mat')
% maxFlux = max(myCSM_exp.minTotal_OFD * 1.05, myCSM_merged.minTotal_OFD * 1.05);
% 
% load('OFD_MILP_myCSM_merged_FVA.mat')
% minBiomass = minval(strcmp(model.rxns,'BIO0010'));
% load('OFD_MILP_myCSM_expression_FVA.mat')
% minBiomass = min(minval(strcmp(model.rxns,'BIO0010')),minBiomass);
% 
% % set the FVA constraints
% myCSM_base.MILP_PFD.lb(strcmp(model.rxns,'BIO0010')) = minBiomass;
% 
% myCSM_base.MILP.lb(strcmp(model.rxns,'BIO0010')) = minBiomass;
% myCSM_base.MILP.b(myCSM_base.MILP.minTotalInd) = maxFlux;
% 
% 
% % run FVA by calling:
% [minval, maxval] = FVA_MILP(myCSM_base.MILP_PFD, model, targetRxns,parforFlag,relMipGapTol,verbose);
% save('PFD_MILP_myCSM_base_FVA.mat','minval','maxval');
% 
% [minval, maxval] = FVA_MILP(myCSM_base.MILP, model, targetRxns,parforFlag,relMipGapTol,verbose);
% save('OFD_MILP_myCSM_base_FVA.mat','minval','maxval');
% 
% 
% %% check by metabolite
% met = 'accoa[m]'; % 'accoa[c]' leu-L
% tbl3 = listRxn(model,myCSM_base.OFD,met); % myCSM_merged.OFD
% 
% %% more visualization and inspection
% fluxTbl = table(model.rxns, myCSM_base.OFD, printRxnFormula(model, model.rxns,0));
% fluxTbl.type = regexprep(fluxTbl.Var1,'[0-9]+','');
% sortedTbl = sortrows(fluxTbl, [4, 2]);
% 
% %% FVA of a single rxn
% targetRxns = {'RC02736'};
% parforFlag = 0;
% relMipGapTol = 1e-3; % this is to be redone with maximum precision, otherwsie the box will have problem
% verbose = false;
% 
% % run FVA by calling:
% [minval, maxval] = FVA_MILP(myCSM_base.MILP, model, targetRxns,parforFlag,relMipGapTol,verbose)
% 
% 
% %% inspection of the flux distribution and FVA 
% fluxTbl = table(model.rxns, myCSM_base.OFD, printRxnFormula(model, model.rxns,0));
% fluxTbl.type = regexprep(fluxTbl.Var1,'[0-9]+','');
% fluxTbl.lb = minval(1:length(model.rxns));
% fluxTbl.ub = maxval(1:length(model.rxns));
% fluxTbl.tightness = ismember(model.rxns,tightRxns);
% fluxTbl = sortrows(fluxTbl, [7, 4, 2],"descend");
% 
% fluxTbl_base = fluxTbl;
