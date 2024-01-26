%% Summary
% this is a wrapper function to run the iMAT-WPS by using various
% combination of input dataset, including the full triple data integration.

%% load dependency and paths
% add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);

%% generate epsilon - only need to run once per model setup 
% regenerate the epsilon when model is modified 
% load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');
% model = configurateModel(model);
% model = changeRxnBounds(model,'EXC0050',-1,'l');% we allow 1 unit of bacteria for flux threshold (epsilon) calculation
% [epsilon_f, epsilon_r,capacity_f,capacity_r] = makeEpsilonSeq_server(model, model.rxns, 0.01, 0.5);
% save('input/model/epsilon_generic_withUptakes.mat','epsilon_f', 'epsilon_r');
% save('input/model/capacity_generic_withUptakes.mat','capacity_f', 'capacity_r');

%% run the five integrations 
% note: we didn't parameterize the iCEL1314 model setup, including 1% storage,
% 1% side metabolites, 0.1% individual side metabolites, and 0.01/50%Vmax
% epsilon (w/ 1 unit of bacteria in epsilon calculation), 10 unit of flux 
% of ATPm and a non-critical iMAT++ parameter, 5% total flux cap for latent
% flux calculation. It will be overcomplicated to parameterize these
% trivial default setups and to test their parameter sensitivity.

% the three critial parameter in the data integration, yield, minLow cap
% and metFit cap were parameterized in the following codes. 

% set parameters 

% the first critical paramter is the biomass yeild. We found that
% maintaining a high yield (e.g., > 0.5) is important to force the model to
% avoid massive wasting of bacterial biomass. The 0.65 is the middle point
% from experimental data range, and sub-optimal compared with theoretical
% maximum yeild (~0.75). 
yield = 0.65;

% run no-data integration
a2_1_iMATpp_no_data(yield)

% run absolute expression integration only
a2_2_iMATpp_only_expression(yield)

% run absolute expression + responsiveness integration
a2_3_iMATpp_dual1_integration(yield)

% run absolute expression + WPS similarity integration
relCap_minLow = 0.05; % this is the relative tolerance of low flux minimization. Default is 5%.
relCap_metFit = 1e-6; % this is the relative tolerance of Loss function fitting for flux split (similarity integration). Default is 5% but here use a sufficiently small value to push for rigid fitting in the dual-integration scheme.
a2_3_iMATpp_dual2_integration(yield, relCap_minLow, relCap_metFit)

% run absolute expression + WPS responsiveness + WPS similarity integration
relCap_minLow = 0.05; % this is the relative tolerance of low flux minimization. Default is 5%.
relCap_metFit = 0.05; % this is the relative tolerance of Loss function fitting for flux split (similarity integration). Default is 5%.
a2_4_iMATpp_triple_integration(yield, relCap_minLow, relCap_metFit)

