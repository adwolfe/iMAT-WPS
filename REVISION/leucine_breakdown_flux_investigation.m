%% note
% this script has to be executed in the root folder (/iMAT_WPS/), not in
% /REVISOION/ folder


%% load dependency and paths
% add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);

%% run the integration

yield = 0.65;

% run absolute expression + WPS responsiveness + WPS similarity integration
relCap_minLow = 0.05; % this is the relative tolerance of low flux minimization. Default is 5%.
relCap_metFit = 0.05; % this is the relative tolerance of Loss function fitting for flux split (similarity integration). Default is 5%.

%% regular integration

% Load model - this is the same across all five integrations 
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat'); % model 
load('./input/model/epsilon_generic_withUptakes.mat'); % precalculated epsilons (see a2_2_run_FVA.m)
load('input/WPS/categ_expression_and_WPS.mat'); % pre-saved gene category data (see a1_analyze_responsiveness_constraints.m)

% change epsilon for Q9 cofactors 
coQ_reactions = {'RMC0028','RM08781','RM08775','RM08774','RM08773','RM08772','RM08771','RM08770','RMC0026',...
                    'RM07273','RC01301','RM01301','RM08767','RC08767','RM01616','RC01616','RC08766','RM08766',...
                    'RM03336','RC03336'};
dolp_reactions = {'RC05556','RCC0029','RCC0030','RCC0031','RC01018'};


epsilon_f(ismember(model.rxns, [coQ_reactions,dolp_reactions])) = 0.001;
epsilon_r(ismember(model.rxns, [coQ_reactions,dolp_reactions])) = 0.001;


branchTbl = readtable('input/WPS/final_branchPoint_table.csv'); % pre-saved flux split data (see a1_analyze_DEsimilarity_constraints.m); must have 'mets', 'rxn1','rxn2', and 'maxCosine'

% setup the model
model = configurateModel(model); % setup model constraints and modifications 
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

% Set parameters
modelType = 2; % 2 for generic C. elegans model. 
speedMode = 1; % lowest speed but highest numerical precision

% run iMAT-WPS with yeild constraint

% set up the yield constraint
% we use the biomass yield rate to constrain the bacteria waste
% this is to force the nutrient to be efficiently used instead of wasted in
% bulk
% yield * V(EXC0050) + V(BIO0010) >= 0 (V(EXC0050) is a negative number)
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


% iMAT-WPS
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
relCap_minLow, relCap_metFit, 1, 1, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);

% note, since we integrated WPS data in both minLow and metFit steps, we
% give them equal priority as well as flexibility in fitting, thus, 5%
% relative cap for both (relCap_minLow, relCap_metFit).

% other parameters were kept default for general iMAT++ integration. 

myCSM_exp_resp_simi_low_eps = myCSM;

%% integratiion with reduced epsilon
check = listRxn(model,myCSM_exp_resp_simi_low_eps.OFD,'hmgcoa[m]');

% make the figure to visualize
load('output/integration_output/myCSM_exp_resp_simi.mat')

% plot leucine breakdown flux in simple barplot

% Define the values to be plotted
values = [myCSM_exp_resp_simi.OFD(strcmp(model.rxns,'RM01360')), ...
            myCSM_exp_resp_simi_low_eps.OFD(strcmp(model.rxns,'RM01360'))];

% Define the labels for the bars
labels = {'default epsilon', 'rescaled epsilon'};

% Create the bar plot
figure;
bar(values);

% Set the labels for the x-axis
set(gca, 'XTickLabel', labels);

% Add title and labels
title(printRxnFormula(model,'RM01360'));
xlabel('Epsilon of coQ and dolp cofactor synthesis reactions');
ylabel('Leucine breakdown flux');


plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [3.5, 3];
plt.LineWidth = 0.5;
plt.FontSize = 8;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/leucine_breakdown_flux_tuning.pdf']);




