function a2_3_iMATpp_dual1_integration(yield)
% this integrates absolute expression + responsiveness information.

% Load model - this is the same across all five integrations 
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');
load('./input/model/epsilon_generic_withUptakes.mat'); % see walkthrough_generic.m for guidance on generating the epsilon values
load('input/WPS/categ_expression_and_WPS.mat');

% setup the model
model = configurateModel(model);
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

% Set parameters
modelType = 2; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
speedMode = 1;

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
1e-5, 1, 1, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);
% note: using IMATplusplus_wiring_triple_inetgration_final with empty
% metabolite table gives equavilent fitting (like total flux and minlow)
% but different alternative OFD. We dont use
% IMATplusplus_wiring_triple_inetgration_final to avoid any potential bug
% related to an empty metabolite constraint in the MILP

myCSM_exp_resp = myCSM;
save('output/integration_output/myCSM_exp_resp.mat','myCSM_exp_resp');

%% codes for manual interactions
% % fast FVA
% addpath scripts/parfor_wait/
% targetRxns = model.rxns;
% parforFlag = 1;
% relMipGapTol = 1e-12; % this is to be redone with maximum precision, otherwsie the box will have problem
% verbose = true;
% 
% % loose the soft constraints
% myCSM_merged.MILP_PFD.b(myCSM_merged.MILP.minLowInd) = myCSM_merged.minLow * 1.05;
% 
% myCSM_merged.MILP.b(myCSM_merged.MILP.minTotalInd) = myCSM_merged.minTotal_OFD * 1.05;
% myCSM_merged.MILP.b(myCSM_merged.MILP.minLowInd) = myCSM_merged.minLow * 1.05;
% 
% % run FVA by calling:
% [minval, maxval] = FVA_MILP(myCSM_merged.MILP_PFD, model, targetRxns,parforFlag,relMipGapTol,verbose);
% save('PFD_MILP_myCSM_merged_FVA.mat','minval','maxval');
% 
% [minval, maxval] = FVA_MILP(myCSM_merged.MILP, model, targetRxns,parforFlag,relMipGapTol,verbose);
% save('OFD_MILP_myCSM_merged_FVA.mat','minval','maxval');
% 
% 
% %% check by metabolite
% met = 'pyr[m]'; % 'accoa[c]' leu-L
% tbl3 = listRxn(model,myCSM_merged.OFD,met); % myCSM_merged.OFD
% 
% %% more visualization and inspection
% fluxTbl = table(model.rxns, myCSM_merged.OFD, printRxnFormula(model, model.rxns,0));
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
% [minval, maxval] = FVA_MILP(myCSM_merged.MILP, model, targetRxns,parforFlag,relMipGapTol,verbose)
% 
% 
% %% make the inspection tabel 
% fluxTbl = table(model.rxns, myCSM_merged.OFD, printRxnFormula(model, model.rxns,0));
% fluxTbl.type = regexprep(fluxTbl.Var1,'[0-9]+','');
% fluxTbl.lb = minval(1:length(model.rxns))';
% fluxTbl.ub = maxval(1:length(model.rxns))';
% fluxTbl.tightness = ismember(model.rxns,tightRxns);
% fluxTbl = sortrows(fluxTbl, [7, 4, 2],"descend");
% 
% fluxTbl_merged = fluxTbl;
% 
