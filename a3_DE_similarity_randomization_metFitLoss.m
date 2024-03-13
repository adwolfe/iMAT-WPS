%% summary
% To assess the consistency between metabolite split constraints and the
% other constriants (responsiveness and abs exp), we randomize the
% metabolite split by generating random couplings for the constrained 
% metabolites and then redo everything again. We use the metabolite loss
% of fit (the objective function in the fitting) as a readout.

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

% Set parameters
modelType = 2; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
speedMode = 1;

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

% fitting metrics in real data 
minLow = myCSM.minLow;
highFit = myCSM.N_highFit/length(myCSM.HGenes);
zeroFit = myCSM.N_zeroFit/length(myCSM.RLNames);
FluxBurdenFit = myCSM.minTotal;
minMet = myCSM.minMetBalanceLoss;

% calculate the total weighted flux (this metric is not used in the
% published ms)
MetFlux = 0;
branchMets = unique(branchTbl.mets);
for ii = 1:length(branchMets)
    % first locate the reactions in the module 
    myrxns = unique([branchTbl.rxn1(strcmp(branchTbl.mets, branchMets{ii})); branchTbl.rxn2(strcmp(branchTbl.mets, branchMets{ii}))]);
    % calculate the total fitness ( cosine(met1) * total absolute flux
    % (met1) )
    MetFlux = MetFlux + max(branchTbl.maxCosine(strcmp(branchTbl.mets, branchMets{ii}))) * sum(abs(model.S(strcmp(model.mets, branchMets{ii}), ismember(model.rxns, myrxns)) .* myCSM.PFD(ismember(model.rxns, myrxns))'));
end
minMet_norm =  myCSM.minMetBalanceLoss / MetFlux; % relative fit
% the relative fit may not be necessary because we kept the metabolite
% being constrianed in the randomization so the direct raw total fit should
% be comparable. (although still different rxns will have diff. flux)


%% randomization of the DE similarity 

% We would like to test if the DE similarity is readily consistent with
% both basal expression and DE responsiveness (in comparason, we only
% tested the consistence between responsiveness and the basal expression
% when randomizaing the responsiveness). This is because only similarity
% itself may be esay to fit and similarity itself is also a weaker
% hypothesis for flux wiring, therefore, we focus on testing its
% consistence with all others together. In this case, we first integrate
% MILP as well as minLow given 5% cap, then calculate the minMet. 

% How to properly randomize the data to form a null that tests the
% chance of random observation is tricky. We reasoned that, the
% randomness being tested here is if the reaction couplings selected by
% the DE similarity (for these metabolites) are meaningful. The
% 'meaningful' is defined by that they presents higher consistency with
% other constriants (eg responsiveness) than random. A good
% definition of such by chance is to randomly associate reaction
% couplings for these metabolites, such that we know the selected
% couplings by DE similarity is more coherent than random. 

% therefore, we randomly choose same number of coupled reaction pairs for 
% EACH metabolite in the final metabolite modeling table, and randomly 
% assign cosine values (without replacement). 

% first generate all possible pair set for each metabolite of interest
load('input/model/capacity_generic_withUptakes.mat');
model_capped = model;
% update the boundary -- so indirectly inreversible reactions are now noted
% in its boundary 
model_capped.ub(~capacity_f) = min(model_capped.ub(~capacity_f),zeros(sum(~capacity_f),1));
model_capped.lb(~capacity_r) = max(model_capped.lb(~capacity_r),zeros(sum(~capacity_r),1));
allMets = unique(branchTbl.mets);
% produce a table of all possible flux coupling pairs 
pairSet = {};
for i = 1:length(allMets)
    % initialize variables
    rxn1 = {};
    rxn2 = {};
    % all associated rxns
    myrxns = model_capped.rxns(any(model_capped.S(strcmp(model_capped.mets,allMets{i}),:),1));
    metInd = strcmp(model_capped.mets,allMets{i});
    % curate each pair
    for ii = 1:(length(myrxns)-1)
        for jj = (1+ii): length(myrxns)
            if (  (model_capped.lb(strcmp(model_capped.rxns,myrxns(ii))) < 0 && model_capped.ub(strcmp(model_capped.rxns,myrxns(ii))) > 0) || ... 
                      (model_capped.lb(strcmp(model_capped.rxns,myrxns(jj))) < 0 && model_capped.ub(strcmp(model_capped.rxns,myrxns(jj))) > 0) ) % any reaction runs in both direction
                rxn1 = [rxn1; myrxns(ii)];
                rxn2 = [rxn2; myrxns(jj)];
            elseif (  ( model_capped.S(metInd,strcmp(model_capped.rxns,myrxns(ii))) * model_capped.lb(strcmp(model_capped.rxns,myrxns(ii))) > 0 || ... 
                            model_capped.S(metInd,strcmp(model_capped.rxns,myrxns(ii))) * model_capped.ub(strcmp(model_capped.rxns,myrxns(ii))) > 0 ...
                            ) && ... % ii producing
                          ( model_capped.S(metInd,strcmp(model_capped.rxns,myrxns(jj))) * model_capped.lb(strcmp(model_capped.rxns,myrxns(jj))) < 0 || ... 
                            model_capped.S(metInd,strcmp(model_capped.rxns,myrxns(jj))) * model_capped.ub(strcmp(model_capped.rxns,myrxns(jj))) < 0 ...
                            )  ... % jj consuming;;
                               ) || ... 
                       (  ( model_capped.S(metInd,strcmp(model_capped.rxns,myrxns(jj))) * model_capped.lb(strcmp(model_capped.rxns,myrxns(jj))) > 0 || ... 
                            model_capped.S(metInd,strcmp(model_capped.rxns,myrxns(jj))) * model_capped.ub(strcmp(model_capped.rxns,myrxns(jj))) > 0 ...
                            ) && ... % jj producing
                          ( model_capped.S(metInd,strcmp(model_capped.rxns,myrxns(ii))) * model_capped.lb(strcmp(model_capped.rxns,myrxns(ii))) < 0 || ... 
                            model_capped.S(metInd,strcmp(model_capped.rxns,myrxns(ii))) * model_capped.ub(strcmp(model_capped.rxns,myrxns(ii))) < 0 ...
                            ) ... % ii consuming;;
                               )
                % one consume and one produce 
                rxn1 = [rxn1; myrxns(ii)];
                rxn2 = [rxn2; myrxns(jj)];
            end
        end
    end
    pairSet{i} = [rxn1, rxn2];
end

% randomization
rng(19951126);
nrand = 10000;

WaitMessage = parfor_wait(nrand, 'Waitbar', true); 
environment = getEnvironment();
parfor i = 1:(nrand+1)
    restoreEnvironment(environment);
    
    % randomize the couplings

    branchTbl_rand = branchTbl;
    % reinitialize it as an empty table
    branchTbl_rand.rxn1 = repmat({''},size(branchTbl_rand,1),1);
    branchTbl_rand.rxn2 = repmat({''},size(branchTbl_rand,1),1);
    branchTbl_rand.maxCosine = nan(size(branchTbl_rand,1),1);
    branchTbl_rand.formula1 = repmat({''},size(branchTbl_rand,1),1);
    branchTbl_rand.formula2 = repmat({''},size(branchTbl_rand,1),1);
    branchTbl_rand.degree = nan(size(branchTbl_rand,1),1);
    % assign permuated values
    permutedCosines = branchTbl.maxCosine(randperm(length(branchTbl.maxCosine)));
    currentCosinePick = 1;
    for zz = 1:length(allMets)
        inds = strcmp(allMets{zz}, branchTbl_rand.mets);
        mypool = pairSet{zz}; % the pool of pairs to pick randomly from
        randomPairs =  mypool(randperm(size(mypool,1), sum(inds)),:);
        branchTbl_rand.rxn1(inds) = randomPairs(:,1);
        branchTbl_rand.rxn2(inds) = randomPairs(:,2);
        branchTbl_rand.maxCosine(inds) = permutedCosines(currentCosinePick:(currentCosinePick + sum(inds) - 1));
        currentCosinePick = currentCosinePick + sum(inds);
    end
    
    % run iMAT-WPS
    try
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
                = IMATplusplus_wiring_triple_inetgration_final(model_coupled,epsilon_f,epsilon_r, ExpCateg, branchTbl_rand, modelType,speedMode,...
                relCap_minLow, relCap_metFit, 1, 0, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);
        
        % since we only look at the fitting metrics, we skip the OFD (by
        % setting doLatent = 0
        % fitting and stops at PFD

        % save the fitting metrics
        fit_rate_high(i) = myCSM.N_highFit/length(myCSM.HGenes);
        fit_rate_zero(i) = myCSM.N_zeroFit/length(myCSM.RLNames);
        fit_rate_minLow(i) = myCSM.minLow;
        fit_rate_minMet(i) = myCSM.minMetBalanceLoss;
        fit_rate_minPFD(i) = myCSM.minTotal;

        % calculate the total fitting fitness of metabolite coupling
        MetFlux = 0;
        branchMets = unique(branchTbl_rand.mets);
        for ii = 1:length(branchMets)
            % first locate the reactions in the module 
            myrxns = unique([branchTbl_rand.rxn1(strcmp(branchTbl_rand.mets, branchMets{ii})); branchTbl_rand.rxn2(strcmp(branchTbl_rand.mets, branchMets{ii}))]);
            % calculate the total fitness ( cosine(met1) * total flux
            % (met1) )
            MetFlux = MetFlux + max(branchTbl_rand.maxCosine(strcmp(branchTbl_rand.mets, branchMets{ii}))) * sum(abs(model.S(strcmp(model.mets, branchMets{ii}), ismember(model.rxns, myrxns)) .* myCSM.PFD(ismember(model.rxns, myrxns))'));
        end
        fit_rate_minMet_norm(i) =  myCSM.minMetBalanceLoss / MetFlux;

    catch 
        fit_rate_high(i) = NaN;
        fit_rate_zero(i) = NaN;
        fit_rate_minLow(i) = NaN;
        fit_rate_minMet(i) = NaN;
        fit_rate_minPFD(i) = NaN;
        fit_rate_minMet_norm(i) = NaN;

    end
    WaitMessage.Send; 

end
WaitMessage.Destroy
save('output/randomization/randomization_DE_similarity_fitting.mat',"fit_rate_high","fit_rate_zero","fit_rate_minLow","fit_rate_minPFD","fit_rate_minMet","fit_rate_minMet_norm");

%% plot figures
% As a QC, make sure there is no NA in fitting (meaning no randomization
% errored out)
any(isnan(fit_rate_minMet))

figure;
hold on
histogram(fit_rate_minMet)
xline(minMet,'--','LineWidth',2,'Color','red')
xlabel('Loss of fit (flux * weight)');
ylabel('Frequency');
title(['p < ',num2str((sum(fit_rate_minMet <= minMet) + 1) / (length(fit_rate_minMet) + 1))])
plt = Plot(); % create a Plot object and grab the current figure
xlim([0 100])
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
% plt.LegendLoc = 'SouthEast';
plt.export(['figures/Randomization_DEsimilarity_fitting.pdf']);

figure;
hold on
histogram(fit_rate_minMet_norm * 100)
xline(minMet_norm * 100,'--','LineWidth',2,'Color','red')
xlabel('Relative loss of fit (%)');
ylabel('Frequency');
title(['p < ',num2str((sum( fit_rate_minMet_norm <= minMet_norm) + 1) / (length(fit_rate_minMet_norm) + 1))])
plt = Plot(); % create a Plot object and grab the current figure
xlim([-1 100])
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
% plt.LegendLoc = 'SouthEast';
plt.export(['figures/Randomization_DEsimilarity_fitting_normalized.pdf']);


figure
histogram(fit_rate_high)
highFit

figure
histogram(fit_rate_zero)
zeroFit

figure
histogram(fit_rate_minPFD)
FluxBurdenFit

figure
histogram(fit_rate_minLow)
minLow


%% NOTE ABOUT PERMUTATION METHOD 

% instead of the permutation set up as above, another idea is to shuffle 
% the cosine matrix (i.e., rows and cols simontanously) and regenerate the
% metabolite table by cutoff at 0.2 etc. This one seems reasonable but is
% flawed, because after shuffling, the constrianed metabolites can be very
% different so it is not testing if WPS data gives consistent coupling,
% rather, it is testing whether structure of the cosine relationship is
% specifically giving a significance higher than random.

% in fact, this methods gives insignificant (weakly sig) p values. The
% interpretation becomes unclear because of the large changes and how the
% fitting metric can be compared to that of real data is questionable
% (given different metabolites being modeled). so we dont use this unclear
% idea.

%% outdated codes to be removed upon publication 
% % load the data needed to produce randomized split table
% DEsim = readtable('./../../MetabolicLibrary/2_DE/output/cosineSimilarity_FC_denoised_stationery_metabolic.csv','ReadRowNames',true);
% lookupTbl = readtable('./input/model/IDtbl.csv');
% condInfo = readtable('./../../MetabolicLibrary/2_DE/output/RNAi_condition_metaInfo.csv');
% load('input/model/capacity_generic_withUptakes.mat');
% model_capped = model;
% % update the boundary -- so indirectly inreversible reactions are now noted
% % in its boundary 
% model_capped.ub(~capacity_f) = min(model_capped.ub(~capacity_f),zeros(sum(~capacity_f),1));
% model_capped.lb(~capacity_r) = max(model_capped.lb(~capacity_r),zeros(sum(~capacity_r),1));
% % generate the metabolite split table
% RNAiGeneNameList = repmat({'NA'},size(DEsim,1),1);
% [A B] = ismember(DEsim.Properties.VariableNames, regexprep(regexprep(condInfo.RNAiID,' ','_'),'\.','_'));
% RNAiGeneNameList(A) = condInfo.RNAi_WBID(B(A));
% [A B] = ismember(RNAiGeneNameList, lookupTbl.WormBase_Gene_ID);
% RNAiGeneNameList(A) = lookupTbl.ICELgene(B(A));

% allBranches = makeSplitTable(model_capped, DEsim, RNAiGeneNameList);

    % randomize table 
%     DEsim_rand = DEsim;
% 
%     if i > 1
%         randInd = randperm(size(DEsim_rand,1));
%         DEsim_rand = DEsim(randInd, randInd);
%     end
%     branchTbl_rand = makeSplitTable(model_capped, DEsim_rand, RNAiGeneNameList);
%     
%     % the table should only contain P/C type 
%     branchTbl_rand = branchTbl_rand(branchTbl_rand.maxCosine > 0.2,:);     
    % NOTE FOR check for metabolites with more than one pairs 
    % in real data analysis, we mannually excluded one metabolite that has
    % multiple coupling modules. This is rare and is minor modification; we
    % ignore it in the randomization (because manual curation is not
    % available in randomization). It is OK becuase the inclusion of such pairs
    % usually help with fitting becuase the more reactions included, the easier
    % to balance the mass flux of a metabolite. But this indeed increase
    % the number of branches to fit for the randomization (178-->187)
    % (we need to do something about this
    % only if the difference in randomization is very small)
    % (if we dont filter for real data, the minMet is 14.84 as compared to
    % the used data 14.54. No big diff.
