%% add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
addpath ./../../MetabolicLibrary/9_FBA_modeling/PlotPub/lib/
initCobraToolbox(false);

%% about - part I
% To justify the responsiveness analysis, the simplest test is to see if
% responsive genes are conflicting with nonresponsive genes at GPR level;
% the following code performs such randomization

% eg. A & B where A is responsive but B is not.

%% randomization
% to measure the conflict level at GPR, we count the number of 
% non-responsive genes (has influence on the reaction) in the responsive 
% reaction (using the GPR parsing script) and do a randomization. 

clear

% first load data
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');
load('input/WPS/categ_expression_and_WPS.mat');

% find out all rxns associated with responsive genes
% in a greedy manner, we consider all reactions associated with a
% responsive gene as the responsive reaction (ignore GPR parsing);
% this is the simplest way to parse responsive genes; otherwise it will be
% complex to determine what reaction should be high reaction to explain a
% high gene, which might be overinterpreting the data. This is the same as
% how responsive gene is used in iMAT++.
resp_rxns = model.rxns(any(model.rxnGeneMat(:,ismember(model.genes, ExpCateg.responsive)),2),1);

% find out all nonresponsive genes that conflict with a responsive reaction
% (gene)
% conflictedGene = {};
N_conflict = 0; % the number of nonresponsive genes that are conflict with at least one responsive gene
for i = 1:length(ExpCateg.nonresponsive)
    [~,haseffect, influencedRxn] = deleteModelGenes(model, ExpCateg.nonresponsive{i});
    if haseffect
        if any(ismember(influencedRxn, resp_rxns))
            N_conflict = N_conflict + 1;
            % conflictedGene = [conflictedGene; ExpCateg.nonresponsive{i}];
        end
    end
end
N_conflict_obs = N_conflict;
% N_conflict = length(conflictedGene);

% randomization 
figure;
genesetValid = union(ExpCateg.responsive, ExpCateg.nonresponsive);
rng(19951126);
nrand = 10000;
N_conflict_rand = [];
WaitMessage = parfor_wait(nrand, 'Waitbar', true); 
parfor zz = 1:nrand
    ExpCateg_rand = ExpCateg;
    permInd = randperm(length(genesetValid));
    ExpCateg_rand.responsive = genesetValid(permInd(1:length(ExpCateg.responsive)));
    ExpCateg_rand.nonresponsive = genesetValid(permInd((1+length(ExpCateg.responsive)):(length(ExpCateg.responsive)+length(ExpCateg.nonresponsive))));
    
    resp_rxns = model.rxns(any(model.rxnGeneMat(:,ismember(model.genes, ExpCateg_rand.responsive)),2),1);
    % conflictedGene = {};
    N_conflict = 0;
    for i = 1:length(ExpCateg_rand.nonresponsive)
        [~,haseffect, influencedRxn] = deleteModelGenes(model, ExpCateg_rand.nonresponsive{i});
        if haseffect
            if any(ismember(influencedRxn, resp_rxns))
                N_conflict = N_conflict + 1;
                % conflictedGene = [conflictedGene; ExpCateg_rand.nonresponsive{i}];
            end
        end
    end
    N_conflict_rand(zz) = N_conflict;
%     if mod(zz,100) == 0
%         fprintf('%d\n',zz)
%     end
    WaitMessage.Send; 

end
WaitMessage.Destroy
save('output/randomization/Randomization_GPR_responsiveness_conflict.mat','N_conflict_rand','N_conflict')

figure;
hold on
histogram(N_conflict_rand./length(ExpCateg.nonresponsive)*100, 'NumBins',34)
xline(N_conflict_obs/length(ExpCateg.nonresponsive)*100,'--','LineWidth',2,'Color','red')
xlabel('Conflicting nonresponsive gene (%)');
ylabel('Frequency');
title(['p < ',num2str((sum( (N_conflict_rand ./ length(ExpCateg.nonresponsive)) <= N_conflict_obs/length(ExpCateg.nonresponsive)) + 1) / (length(N_conflict_rand) + 1))])
plt = Plot(); % create a Plot object and grab the current figure
xlim([-0.5 15])
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
% plt.LegendLoc = 'SouthEast';
plt.export(['figures/Randomization_GPR_responsiveness_conflict.pdf']);


%% about - part II
% this is to test the randomization of responsiveness annotation; we use
% the total flux through nonresponsive reactions (defined by the reactions
% whose flux can be directly impacted by the in silico deletion of any 
% nonresponsive gene) as a readout to assess the level of conflict between
% high genes (responsive and high expression) and the nonresponsive genes. 

clear;
close all;
% parameters 
yield = 0.65;
relCap_minLow = 0.05;
relCap_metFit = 0.05;

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
myCSM.wasteDW]...
= IMATplusplus_wiring_triple_inetgration_final(model_coupled,epsilon_f,epsilon_r, ExpCateg, branchTbl, modelType,speedMode,...
relCap_minLow, relCap_metFit, 1, 1, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);

% to strictly control the randomization, we use the setup for full modeling
% integration and keep all parameters the same;

% calculate the total nonresponsive flux
relatedRxns = [];
for i = 1:length(ExpCateg.nonresponsive)
    [~,haseffect, influencedRxn] = deleteModelGenes(model, ExpCateg.nonresponsive{i});
    if haseffect
        relatedRxns = union(relatedRxns, influencedRxn);
    end
end
totalNonrespFlux_obs = sum(abs(myCSM.OFD(ismember(model_coupled.rxns,relatedRxns))));
totalFlux_obs = sum(abs(myCSM.OFD));

%% randomization of the responsiveness annotation
genesetValid = union(ExpCateg.responsive, ExpCateg.nonresponsive); % all tested genes in WPS (pass filter)
rng(19951126);
nrand = 10000;

WaitMessage = parfor_wait(nrand, 'Waitbar', true); 
environment = getEnvironment();
parfor i = 1:nrand
    restoreEnvironment(environment);
    
    ExpCateg_rand = ExpCateg;
    % We only randomly assign responsive/non-responsive within
    % the tested conditions (i.e., randomly permute the original responsive 
    % and nonreponsive genes)
    permInd = randperm(length(genesetValid));
    ExpCateg_rand.responsive = genesetValid(permInd(1:length(ExpCateg.responsive)));
    ExpCateg_rand.nonresponsive = genesetValid(permInd((1+length(ExpCateg.responsive)):(length(ExpCateg.responsive)+length(ExpCateg.nonresponsive))));
    try % in case of infeasible models
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
                myCSM.wasteDW]...
                = IMATplusplus_wiring_triple_inetgration_final(model_coupled,epsilon_f,epsilon_r, ExpCateg_rand, branchTbl, modelType,speedMode,...
                relCap_minLow, relCap_metFit, 1, 1, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);
       
        % calculate the total nonresponsive flux
        relatedRxns = [];
        for kk = 1:length(ExpCateg_rand.nonresponsive)
            [~,haseffect, influencedRxn] = deleteModelGenes(model, ExpCateg_rand.nonresponsive{kk});
            if haseffect
                relatedRxns = union(relatedRxns, influencedRxn);
            end
        end
        totalNonrespFlux(i) = sum(abs(myCSM.OFD(ismember(model_coupled.rxns,relatedRxns))));
        totalFlux(i) = sum(abs(myCSM.OFD));

    catch 
        totalNonrespFlux(i) = NaN;
        totalFlux(i) = NaN;

    end
    WaitMessage.Send; 
end
WaitMessage.Destroy
save('output/randomization/randomization_responsiveness_reaction_level.mat',"totalNonrespFlux","totalFlux");


%% plots 
% check if there is any infeasible (dont really expect so)
any(isnan(totalNonrespFlux))

figure;
hold on
histogram(totalNonrespFlux, 'NumBins',200)
xline(totalNonrespFlux_obs,'--','LineWidth',2,'Color','red')
xlabel('Total conflicting flux through nonresponsive reactions');
ylabel('Frequency');
title(['p < ',num2str((sum(totalNonrespFlux <= totalNonrespFlux_obs) + 1) / (length(totalNonrespFlux) + 1))])
plt = Plot(); % create a Plot object and grab the current figure
xlim([0 160]) % some extreme values were cutoff
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
% plt.LegendLoc = 'SouthEast';
plt.export(['figures/Randomization_responsiveness_nonresponsive_reaction_absolute_flux.pdf']);

figure;
hold on
histogram(totalNonrespFlux ./ totalFlux * 100)
xline(totalNonrespFlux_obs/totalFlux_obs * 100,'--','LineWidth',2,'Color','red')
xlabel('Relative total flux through nonresponsive reactions (total flux %) (%)');
ylabel('Frequency');
title(['p < ',num2str((sum( (totalNonrespFlux ./ totalFlux) <= totalNonrespFlux_obs/totalFlux_obs) + 1) / (length(totalNonrespFlux) + 1))])
plt = Plot(); % create a Plot object and grab the current figure
xlim([-1 30])
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
% plt.LegendLoc = 'SouthEast';
plt.export(['figures/Randomization_responsiveness_nonresponsive_reaction_normalized_flux.pdf']);

