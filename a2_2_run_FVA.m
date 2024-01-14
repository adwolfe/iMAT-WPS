%% load dependency
% add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
addpath scripts/parfor_wait/

initCobraToolbox(false);
%% run FVA of all integrations 
% load model
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');
model = configurateModel(model);

% set FVA parameters
targetRxns = model.rxns;
parforFlag = 1;
relMipGapTol = 1e-12; % this is to be redone with maximum precision, otherwsie the box will have problem
verbose = false;
% to assess reasonable parameter sensitive range (for FVA), we allow 5%
% flexibility of the fitting optimality
FVA_rel_cap = 1.05;

%% the expression only integration  
% load the CSM
load('output/integration_output/myCSM_expression_only.mat')

% loose the soft constraints
% minlow - PFD
myCSM_exp_only.MILP_PFD.b(myCSM_exp_only.MILP.minLowInd) = myCSM_exp_only.minLow * FVA_rel_cap;
% minTotal and minLow - OFD
myCSM_exp_only.MILP.b(myCSM_exp_only.MILP.minTotalInd) = myCSM_exp_only.minTotal_OFD * FVA_rel_cap;
myCSM_exp_only.MILP.b(myCSM_exp_only.MILP.minLowInd) = myCSM_exp_only.minLow * FVA_rel_cap;

% run FVA by calling:
[minval, maxval] = FVA_MILP(myCSM_exp_only.MILP_PFD, model, targetRxns,parforFlag,relMipGapTol,verbose);
save('output/integration_output/FVA_PFD_expression_only.mat','minval','maxval');

[minval, maxval] = FVA_MILP(myCSM_exp_only.MILP, model, targetRxns,parforFlag,relMipGapTol,verbose);
save('output/integration_output/FVA_OFD_expression_only.mat','minval','maxval');

%% absolute expression + responsiveness
% load the CSM
load('output/integration_output/myCSM_exp_resp.mat')

% loose the soft constraints
% minlow - PFD
myCSM_exp_resp.MILP_PFD.b(myCSM_exp_resp.MILP.minLowInd) = myCSM_exp_resp.minLow * FVA_rel_cap;
% minTotal and minLow - OFD
myCSM_exp_resp.MILP.b(myCSM_exp_resp.MILP.minTotalInd) = myCSM_exp_resp.minTotal_OFD * FVA_rel_cap;
myCSM_exp_resp.MILP.b(myCSM_exp_resp.MILP.minLowInd) = myCSM_exp_resp.minLow * FVA_rel_cap;

% run FVA by calling:
[minval, maxval] = FVA_MILP(myCSM_exp_resp.MILP_PFD, model, targetRxns,parforFlag,relMipGapTol,verbose);
save('output/integration_output/FVA_PFD_exp_resp.mat','minval','maxval');

[minval, maxval] = FVA_MILP(myCSM_exp_resp.MILP, model, targetRxns,parforFlag,relMipGapTol,verbose);
save('output/integration_output/FVA_OFD_exp_resp.mat','minval','maxval');

%% absolute expression + DE similarity 
% load the CSM
load('output/integration_output/myCSM_exp_simi.mat')

% loose the soft constraints
% minlow - PFD
myCSM_exp_simi.MILP_PFD.b(myCSM_exp_simi.MILP.minLowInd) = myCSM_exp_simi.minLow * FVA_rel_cap; % same as current (already set)
% metFit - PFD
myCSM_exp_simi.MILP_PFD.b(myCSM_exp_simi.MILP.minMetBalanceLossInd) = myCSM_exp_simi.minMetBalanceLoss * FVA_rel_cap; % loosen it

% minTotal and minLow - OFD
myCSM_exp_simi.MILP.b(myCSM_exp_simi.MILP.minTotalInd) = myCSM_exp_simi.minTotal_OFD * FVA_rel_cap; % same as current (already set)
myCSM_exp_simi.MILP.b(myCSM_exp_simi.MILP.minLowInd) = myCSM_exp_simi.minLow * FVA_rel_cap; % same as current (already set)
% metFit - OFD
myCSM_exp_simi.MILP.b(myCSM_exp_simi.MILP.minMetBalanceLossInd) = myCSM_exp_simi.minMetBalanceLoss * FVA_rel_cap; % loosen it

% run FVA by calling:
[minval, maxval] = FVA_MILP(myCSM_exp_simi.MILP_PFD, model, targetRxns,parforFlag,relMipGapTol,verbose);
save('output/integration_output/FVA_PFD_exp_simi.mat','minval','maxval');

[minval, maxval] = FVA_MILP(myCSM_exp_simi.MILP, model, targetRxns,parforFlag,relMipGapTol,verbose);
save('output/integration_output/FVA_OFD_exp_simi.mat','minval','maxval');

%% absolute expression + responsiveness + DE similarity 
% load the CSM
load('output/integration_output/myCSM_exp_resp_simi.mat')

% although the model is constrained with 5% caps already, we reimpose it to
% enable parameter control in th
% loose the soft constraints
% minlow - PFD
myCSM_exp_resp_simi.MILP_PFD.b(myCSM_exp_resp_simi.MILP.minLowInd) = myCSM_exp_resp_simi.minLow * FVA_rel_cap; % same as current (already set)
% metFit - PFD
myCSM_exp_resp_simi.MILP_PFD.b(myCSM_exp_resp_simi.MILP.minMetBalanceLossInd) = myCSM_exp_resp_simi.minMetBalanceLoss * FVA_rel_cap; % loosen it

% minTotal and minLow - OFD
myCSM_exp_resp_simi.MILP.b(myCSM_exp_resp_simi.MILP.minTotalInd) = myCSM_exp_resp_simi.minTotal_OFD * FVA_rel_cap; % same as current (already set)
myCSM_exp_resp_simi.MILP.b(myCSM_exp_resp_simi.MILP.minLowInd) = myCSM_exp_resp_simi.minLow * FVA_rel_cap; % same as current (already set)
% metFit - OFD
myCSM_exp_resp_simi.MILP.b(myCSM_exp_resp_simi.MILP.minMetBalanceLossInd) = myCSM_exp_resp_simi.minMetBalanceLoss * FVA_rel_cap; % loosen it

% run FVA by calling:
[minval, maxval] = FVA_MILP(myCSM_exp_resp_simi.MILP_PFD, model, targetRxns,parforFlag,relMipGapTol,verbose);
save('output/integration_output/FVA_PFD_exp_resp_simi.mat','minval','maxval');

[minval, maxval] = FVA_MILP(myCSM_exp_resp_simi.MILP, model, targetRxns,parforFlag,relMipGapTol,verbose);
save('output/integration_output/FVA_OFD_exp_resp_simi.mat','minval','maxval');


%% The no-data integration
% load the CSM
load('output/integration_output/myCSM_no_data.mat')

% we use the base model as a control, therefore, the FVA constraints were
% given by the two reference model

% specifically, we set the total flux constraint of OFD FVA with the
% max total flux of the four integrations
load('output/integration_output/myCSM_expression_only.mat')
load('output/integration_output/myCSM_exp_resp.mat')
load('output/integration_output/myCSM_exp_simi.mat')
load('output/integration_output/myCSM_exp_resp_simi.mat')
maxFlux = max([myCSM_exp_only.minTotal_OFD;
                myCSM_exp_resp.minTotal_OFD;
                myCSM_exp_simi.minTotal_OFD;
                myCSM_exp_resp_simi.minTotal_OFD])* 1.05;

% to better align the FVA interval with the flux scale in the data
% integration groups, we need to set a minimal biomass production
% constraint properly (instead of keeping the 0.8 as in the iMAT++)
% we use the minimal biomass production in PFD FVA of all four data 
% integrations
load('output/integration_output/FVA_PFD_expression_only.mat')
minvals = minval(strcmp(model.rxns,'BIO0010'));
load('output/integration_output/FVA_PFD_exp_resp.mat')
minvals = [minvals; minval(strcmp(model.rxns,'BIO0010'))];
load('output/integration_output/FVA_PFD_exp_simi.mat')
minvals = [minvals; minval(strcmp(model.rxns,'BIO0010'))];
load('output/integration_output/FVA_PFD_exp_resp_simi.mat')
minvals = [minvals; minval(strcmp(model.rxns,'BIO0010'))];
minBiomass = min(minvals);

% set the FVA constraints
% min biomass production - PFD
myCSM_no_data.MILP_PFD.lb(strcmp(model.rxns,'BIO0010')) = minBiomass;

% min biomass and max total flux - PFD
myCSM_no_data.MILP.lb(strcmp(model.rxns,'BIO0010')) = minBiomass;
myCSM_no_data.MILP.b(myCSM_no_data.MILP.minTotalInd) = maxFlux;

% run FVA by calling:
[minval, maxval] = FVA_MILP(myCSM_no_data.MILP_PFD, model, targetRxns,parforFlag,relMipGapTol,verbose);
save('output/integration_output/FVA_PFD_no_data.mat','minval','maxval');

[minval, maxval] = FVA_MILP(myCSM_no_data.MILP, model, targetRxns,parforFlag,relMipGapTol,verbose);
save('output/integration_output/FVA_OFD_no_data.mat','minval','maxval');

% note: this FVA is essentially the same FVA on a LP constrained with the
% same min biomass and max total flux 
