%% Identify the key data that drives the cyclic PPP 

% This script is just to repeat the LOO analysis with targeted reaction
% (the gspd-1 reaction) and check for the genes whose responsiveness have
% strongest impact on the cyclic PPP. This is motivated by seeing only
% integrating WPS responsiveness + expression levels is sufficient to
% narrow down the solution space of PPP to high cyclic flux. 

% parameters 
yield = 0.65;

% the reaction of interest
targetRxns = {'RC02736'};
load('input/WPS/categ_expression_and_WPS.mat');

queryGenes = [ExpCateg.nonresponsive;ExpCateg.responsive];
mytype = 1; % 1 for LOO and 2 for LOI

%% add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
addpath ./../../MetabolicLibrary/9_FBA_modeling/PlotPub/lib/
initCobraToolbox(false);
%% Load model
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');
load('./input/model/epsilon_generic_withUptakes.mat');
load('input/WPS/categ_expression_and_WPS.mat');

% setup the model
model = configurateModel(model);
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

% SPECIAL TREATMENT FOR ANALYZING PPP - actually doesnt matter for
% analyzing gspd-1 reaction
% we noticed that the FVA of PPP back flux is confounded by the two
% alternatives going to g6p-A and g6p-B. therefore, we block one of the two
% reaction to reveal actual results.
% block the f6p to g6p-A to study valid flux variability 
model = changeRxnBounds(model,'RC02740',0,'b'); 

% Set parameters
modelType = 2; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
speedMode = 1;

% analyze the FVA with candidate responsiveness information removed

% run iMAT-WPS with yeild constraint

% set up the yield constraint
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

% basic FVA parameters
parforFlag = 0;
relMipGapTol = 1e-12; % this is to be redone with maximum precision, otherwsie the box will have problem
verbose = false;

% start program - perform LOO analysis for each responsiveness annotation
h_overall = waitbar(0,'Starting analyzing responsiveness dependency...');
Maxs = nan(length(queryGenes), length(targetRxns));
Mins = nan(length(queryGenes), length(targetRxns));

OFD_ubs = nan(length(queryGenes), length(targetRxns));
OFD_lbs = nan(length(queryGenes), length(targetRxns));

OFDs = nan(length(queryGenes), length(targetRxns));

for zz = 1:length(targetRxns)
    
    targetRxn = targetRxns{zz}; 
  
    for yy = mytype
        environment = getEnvironment();
        WaitMessage = parfor_wait(length(queryGenes), 'Waitbar', true); 
        parfor ii = 1:(length(queryGenes)+1)
            restoreEnvironment(environment);
            
            ExpCateg_tmp = ExpCateg;
            responsive = ExpCateg_tmp.responsive;
            nonresponsive = ExpCateg_tmp.nonresponsive;
            if ii > 1 % the frist one is to reproduce FVA of original model
                if yy == 1 % perform LOO analysis
                    responsive(ismember(responsive,queryGenes{ii-1})) = [];
                    nonresponsive(ismember(nonresponsive,queryGenes{ii-1})) = [];
                elseif yy == 2 % perform Leave One In (LOI) analysis
                    responsive(~ismember(responsive,queryGenes{ii-1})) = [];
                    nonresponsive(~ismember(nonresponsive,queryGenes{ii-1})) = [];
                end
            end
            ExpCateg_tmp.responsive = responsive;
            ExpCateg_tmp.nonresponsive = nonresponsive;
        
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
            = IMATplusplus_wiring_dual_integration_final(model_coupled,epsilon_f,epsilon_r, ExpCateg_tmp, modelType,speedMode,...
            1e-5, 1, 1, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);
            % we only focus on analyzing dual integration here

            % FVA
            % loose the soft constraint
            myCSM.MILP_PFD.b(myCSM.MILP.minLowInd) = myCSM.minLow * 1.05;
            [minval1, maxval1] = FVA_MILP(myCSM.MILP_PFD, model_coupled, {targetRxn},parforFlag,relMipGapTol,verbose);
            
            myCSM.MILP.b(myCSM.MILP.minTotalInd) = myCSM.minTotal_OFD * 1.05;
            myCSM.MILP.b(myCSM.MILP.minLowInd) = myCSM.minLow * 1.05;
            [minval2, maxval2] = FVA_MILP(myCSM.MILP, model_coupled, {targetRxn},parforFlag,relMipGapTol,verbose);


            OFDs(ii,zz) = myCSM.OFD(strcmp(model_coupled.rxns, targetRxn));
            Maxs(ii,zz) = maxval1;
            Mins(ii,zz) = minval1;
            OFD_ubs(ii,zz) = maxval2;
            OFD_lbs(ii,zz) = minval2;

            
            WaitMessage.Send; 
        
        end
        WaitMessage.Destroy
       
    end
   
    waitbar(zz/length(targetRxns),h_overall)
end


%% plot the result
OFDs = OFDs';
Maxs = Maxs';
Mins = Mins';
OFD_ubs = OFD_ubs';
OFD_lbs = OFD_lbs';

% these genes caused a big change in solution space (i.e., 50% reduction)
queryGenes(Mins(2:end) < Mins(1)*0.5)
% this drives us to focus on K07E3.4 and idh-1, together with inspecting
% the heatmap of systematic LOO and LOI analysis. 





