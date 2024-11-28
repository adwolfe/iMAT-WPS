function cluster_a2_sensitivity_analysis_sampling(tmpDir,batchID_i, batchID_j)
%% about 
% re-do the full iMAT-WPS modeling with subsampled input 

% we use the disassimilation (yield) rate to constrain the bacteria waste
disAssim = 0.65;

%% Step 1: make the OFDs
% add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/10.0.0/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
addpath integration_pipelines/

load([tmpDir,'/environment.mat']);
restoreEnvironment(environment);

%% Load model
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');
load('./input/model/epsilon_generic_withUptakes.mat'); % see walkthrough_generic.m for guidance on generating the epsilon values

% load('input/WPS/categ_expression_and_WPS.mat');
% branchTbl = readtable('input/WPS/final_branchPoint_table.csv'); % must have 'mets', 'rxn1','rxn2', and 'maxCosine'

% setting up the model 
model = configurateModel(model);
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

% Set parameters
modelType = 2; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
speedMode = 1;

%% regenerate the WPS input based on the subset data 
conditionInfo = sampleData{str2num(batchID_i), str2num(batchID_j)};

ExpCateg = extract_responsiveness_constraints(conditionInfo);
branchTbl = extract_DEsimilarity_constraints(conditionInfo);


%% rerun the entire WPS-IMAT analysis

model_coupled = model;
% add the disassimilation constraints 
model_coupled.S(end+1,:) = zeros(1,length(model_coupled.rxns));
model_coupled.S(end, strcmp('EXC0050',model_coupled.rxns)) = disAssim; 
model_coupled.S(end, strcmp('BIO0010',model_coupled.rxns)) = 1; 
model_coupled.csense(end+1) = 'G';
model_coupled.b(end+1) = 0;
model_coupled.BiGG(end+1) = {'NA'};
model_coupled.metCharges(end+1) = 0;
model_coupled.metFormulas(end+1) = {'NA'};
model_coupled.metNames(end+1) = {'Pseudo-metabolite that represents a constraint'};
model_coupled.mets(end+1) = {['NonMetConst',num2str(length(model_coupled.mets))]};

    
% run integration
if samplingDepth(str2num(batchID_i)) ~= 0
    % WPS data exist
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
    myCSM.branchMets]...
    = IMATplusplus_wiring_triple_inetgration_final(model_coupled,epsilon_f,epsilon_r, ExpCateg, branchTbl, modelType,speedMode,...
    0.05, 0.05, 1, 1, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);
else
    % when there is no WPS data, we aggressively fit gene expression (as
    % done before) - use 1e-5 absolute minlow cap
    % just modify the format of the empty cell
    ExpCateg.responsive = {};
    ExpCateg.nonresponsive = {};

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
end

myCSM_ori = myCSM;

% FVA
% set FVA parameters
targetRxns = model_coupled.rxns; % for all reactions
parforFlag = 0;
relMipGapTol = 1e-12; % use most stringent MILP gap tolerence in increase numerical precision 
verbose = false; % minimal output
% to assess reasonable parameter sensitive range (for FVA), we allow 5%
% flexibility of the fitting optimality
FVA_rel_cap = 1.05; % we always allow 5% relative fitting cap for each fitting objective


% although the model is constrained with 5% caps already, we re-impose it to
% enable parameter control in the future.
% loose the soft constraints
% minlow - PFM
myCSM.MILP_PFD.b(myCSM.MILP.minLowInd) = myCSM.minLow * FVA_rel_cap; % same as current (already set)
% metFit - PFM
if samplingDepth(str2num(batchID_i)) ~= 0
    myCSM.MILP_PFD.b(myCSM.MILP.minMetBalanceLossInd) = myCSM.minMetBalanceLoss * FVA_rel_cap; % loosen it
end

% minTotal and minLow - OFM
myCSM.MILP.b(myCSM.MILP.minTotalInd) = myCSM.minTotal_OFD * FVA_rel_cap; % same as current (already set)
myCSM.MILP.b(myCSM.MILP.minLowInd) = myCSM.minLow * FVA_rel_cap; % same as current (already set)
% metFit - OFM
if samplingDepth(str2num(batchID_i)) ~= 0
    myCSM.MILP.b(myCSM.MILP.minMetBalanceLossInd) = myCSM.minMetBalanceLoss * FVA_rel_cap; % loosen it
end

% run FVA by calling:
[minval_PFM, maxval_PFM] = FVA_MILP(myCSM.MILP_PFD, model_coupled, targetRxns,parforFlag,relMipGapTol,verbose);

[minval_OFM, maxval_OFM] = FVA_MILP(myCSM.MILP, model_coupled, targetRxns,parforFlag,relMipGapTol,verbose);

save([tmpDir,'/iMAT_WPS_subsample_',batchID_i,'_', batchID_j,'.mat'],'myCSM_ori','minval_PFM','maxval_PFM','minval_OFM','maxval_OFM');
end