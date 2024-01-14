%% story of PPP

% dissect the prediction of PPP 
% top flux reactions and novel wiring pattern --> driven by responsiveness
% --> two nonresponsive genes blocking alternative pathways -->
% surprisingly, independently supported by similarity --> confirmed by
% modeling (removing the two on top of similarity constraint does not
% remove the OFD and range stays high) --> high confience for validation

% next, we analyze it under the full model to see the effect from
% similarity

% parameters 
yield = 0.65;
relCap_minLow = 0.05;
relCap_metFit = 0.05;

% the reaction of interest
targetRxns = {'RC02736','RC02035','RC01528','RC01529','RC01641','RC01830','RC01827','RC03321'};
queryGenes = {{'K07E3.4'},{'idh-1'},{'K07E3.4', 'idh-1'}};
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
load('./input/model/epsilon_generic_withUptakes.mat'); % see walkthrough_generic.m for guidance on generating the epsilon values
load('input/WPS/categ_expression_and_WPS.mat');
branchTbl = readtable('input/WPS/final_branchPoint_table.csv'); % must have 'mets', 'rxn1','rxn2', and 'maxCosine'

% setup the model
model = configurateModel(model);
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

% SPECIAL TREATMENT FOR ANALYZING PPP
% we noticed that the FVA of PPP back flux is confounded by the two
% alternatives going to g6p-A and g6p-B. therefore, we block one of the two
% reaction to reveal actual results.
% block the f6p to g6p-A to study valid flux variability 
model = changeRxnBounds(model,'RC02740',0,'b'); 

% Set parameters
modelType = 2; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
speedMode = 1;
allMets = unique(branchTbl.mets);

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

% basic FVA parameters
parforFlag = 0;
relMipGapTol = 1e-12; % this is to be redone with maximum precision, otherwsie the box will have problem
verbose = false;
    
% start program
h_overall = waitbar(0,'Starting analyzing similarity dependency...');
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
            = IMATplusplus_wiring_triple_inetgration_final(model_coupled,epsilon_f,epsilon_r, ExpCateg_tmp, branchTbl, modelType,speedMode,...
            relCap_minLow, relCap_metFit, 1, 1, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);
            % the above caps (relCap_minLow, relCap_metFit) is used
            % directly in FVA

            [minval1, maxval1] = FVA_MILP(myCSM.MILP_PFD, model_coupled, {targetRxn},parforFlag,relMipGapTol,verbose);
            
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

% checkTbl = branchTbl(ismember(branchTbl.mets, driverMets),:);

%% plot the result
mytext = strcat(targetRxns,':',{' '}, printRxnFormula(model,targetRxns,0));
OFDs = OFDs';
Maxs = Maxs';
Mins = Mins';
OFD_ubs = OFD_ubs';
OFD_lbs = OFD_lbs';

% Create the horizontal bar plot
h = figure;

num_plots = size(OFDs,1); % Number of plots
colors = [0.8500 0.3250 0.0980;
         0.9290 0.6940 0.1250;
         0.4660 0.6740 0.1880;
         0.4940 0.1840 0.5560;
         0 0.4470 0.7410;
         0.3010 0.7450 0.9330;
         0.6350 0.0780 0.1840];


% We need to adjust the bar and the errorbar calls to loop over each group
for p = 1:num_plots % Loop over plots
    subplot(num_plots, 1, p);
    hold on;
    % xline(0,'-b', 'LineWidth',2)

    for i = 1:size(OFDs, 2)
        fill( [Mins(p, i); Maxs(p, i); Maxs(p, i); Mins(p, i)], [i-0.4 i-0.4 i+0.4 i+0.4], colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 1);
        fill( [OFD_lbs(p, i); OFD_ubs(p, i); OFD_ubs(p, i); OFD_lbs(p, i)], [i-0.4 i-0.4 i+0.4 i+0.4], nan, 'EdgeColor', [0,0,0],'LineWidth', 0.5);
        line([OFDs(p, i) OFDs(p, i)], [i-0.4 i+0.4],   'Color', 'r', 'LineWidth', 2); % Line from lower bound to height        
    end
    text(0, 3.5, '*', 'Color', 'blue', 'HorizontalAlignment', 'center','FontSize',25)

    hold off;
    title(mytext(p))

    % set the x-axis
    xlim([-1.4, 1.4])
    tmp = cellfun(@(x) strjoin(x, '\\Delta'), queryGenes, 'UniformOutput', false);
    set(gca,'YTickLabel', [{'original'},strcat('\Delta', tmp)],'YTick', 1:size(OFDs, 2))

end
xlabel('Flux')

h.Position(4) = 100*length(targetRxns);  % for example, set the height to 800
saveas(h, 'figures/Case_study_PPP_similarity_interaction_mechanisms.pdf');



