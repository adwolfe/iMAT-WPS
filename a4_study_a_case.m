%% about the mechanism 
% there are four types of mechansim analyses, each with its own pros and
% cons in understanding the prediction, due to the complex interactions
% between exp, responsiveness and similarity constraints. The four are: LOO
% and LOI on a exp + resp. model; LOO and LOI on a exp + sim model; LOO and
% LOI of resp. on a exp+resp+sim model; and LOO and LOI of sim on a exp+resp+sim model;

% the first two reveals the information that resp and sim can add into a
% basic exp model, respectively. The last two reveals that for a triple
% model with interaction contexts. For example, sim and resp data may
% constrian the same reaction in the same way, therefore, the LOO of any
% will be buffered unless only use a dual model (first two). LOI may not
% suffer from the problem but it lost context from other genes or mets so
% can be misleading (also interacting with the cap). 

% so there is no perfect way to systemactially reveal mechanisms for all
% cases. The only thing we can do is to chose a way for tentative mechanism
% for a set of reactions; The only way to figure out details is to do case
% study of a reaction for all possible tests (including more contexted
% model and more than LOO (eg leave multiple out) or educated guess. 


%% note about FVA analysis

% issue: the OFD struture overestimates solution space because of latent
% constriants (minor driver) and total flux cap (major driver). flux cap is
% strongly dependent on both MILP and minLow cap, and even possibly the
% metBalance cap. 

% solution: assessing solution space independent of flux cap 
% FVA performed at OFD MILP level is dependent on the latent rxn
% constriants which are dependent on the total flux minimalization at PFD;
% therefore, the different integration (single, dual, triple) can be
% different and thus not comparable. The strict FVA to assess data-driven
% solution space is the FVA on PFD level that is independent of flux
% minimization and any related constraint (PFD MILP doesnt have a mintotal
% constraint in place). This FVA is fully data driven (no flux minimization
% caps) and is comparable across different integrations (the MILP structure
% stays the same; essentially we can use the triple integration MILP and
% set free the minLow or/and MetFlux constriants to run FVA on differnet
% integration; it will return identical result as that from a seperate
% iMAT++ run).

% OFD --> good for getting single flux solution (best guess of flux)
% PFD MILP --> good for FVA (data driven narrrowing of solution space)

%% about LOO and LOI
% LOO is more deterministic about prediction mechanism since if an effect
% is lost after LOO, it has to be an direct effect of the constraint;
% However, LOI is often suggestive and can be misleading. This is because
% how a contraint is contributing to the solution space is dependent on the
% other contriant contexts; for example, a nonresponsive gene may be
% overide by another responsive gene in whole data modeling but being
% impactful in LOI. same for the metabolite balance, a metabolite split can
% be simply tolerated when other constaint exit and given the 5% cap,
% however, cannot be tolerated when only fitting one constraint. So LOI is
% only suggestive and may indicates some mechanism that may not be true. 

% the low coverage of LOO mechanisms indicate the constraints are redundant
% and can be interacting.

%% PART I: THE APPLICATION TO THE GENERIC C. ELEGANS MODEL
% run through the PMNM integration of iMAT++ and FVA/network building
%% Step 1: make the OFDs
% add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);
%% Load model
load('./../input/makeWormModel/iCEL1314_withUptakes.mat');

% setting up the model 
% block a reaction that is a misannotation in the model and should be not
% existing 
model = changeRxnBounds(model,'TCM1071',0,'b'); % remove the nadp version

% (1) F45H10.3 is down and highly expressed, but has no response.
% considering this is the only exception of ETC gene (conflict with
% data),we remove this gene from model (so it might be an error)
% [model, affectedrxns] = removeGenesFromModel_XL(model,{'F45H10.3'});  -- no longer needed

% block artificial loops that can convert nadh to nadph without energy
% cost 
model = changeRxnBounds(model,'RC01758',0,'u'); % only allow backward flux per BiGG
model = changeRxnBounds(model,'RC01759',0,'u'); % only allow backward flux per BiGG
model = changeRxnBounds(model,'EX00532',-0.01,'l'); % resuce blocked flux
model = changeRxnBounds(model,'RC08379',0,'u'); % arbiturary based on map
model = changeRxnBounds(model,'RC02124',0,'u'); % arbiturary based on map
model = changeRxnBounds(model,'RC01093',0,'l'); % arbiturary based on map
model = changeRxnBounds(model,'RC01095',0,'l'); % arbiturary based on map
model = changeRxnBounds(model,'RC01036',0,'u'); % arbiturary based on map
model = changeRxnBounds(model,'RC01041',0,'u'); % arbiturary based on map
model = changeRxnBounds(model,'RC02235',0,'u'); % assume dhf is always being comsumed
model = changeRxnBounds(model,'RC02236',0,'u'); % assume dhf is always being comsumed
model = changeRxnBounds(model,'RC00936',0,'l'); % assume dhf is always being comsumed
model = changeRxnBounds(model,'RC00939',0,'l'); % assume dhf is always being comsumed
model = changeRxnBounds(model,'RC01224',0,'l'); % arbiturary based on folate/sam cycle
model = changeRxnBounds(model,'RC07168',0,'l'); % arbiturary based on folate/sam cycle
model = changeRxnBounds(model,'RC03302',0,'l'); % arbiturary based on map
model = changeRxnBounds(model,'RC00100',0,'l'); % only forward is feasible
model = changeRxnBounds(model,'RC08539',0,'l'); % only forward is feasible
model = changeRxnBounds(model,'RC01218',0,'l'); % arbiturary to avoid nadh --> nadph
model = changeRxnBounds(model,'RC02695',0,'l'); % arbiturary based on map
model = changeRxnBounds(model,'RC02697',0,'l'); % arbiturary based on map
model = changeRxnBounds(model,'RC07140',0,'b'); % remove the nadp version
model = changeRxnBounds(model,'RM03293',0,'l'); % only forward is feasible
model = changeRxnBounds(model,'RM03291',0,'l'); % only forward is feasible
model = changeRxnBounds(model,'RM00248',0,'l'); % assume glutamate dehydrogenase only in forward direction (making nh4) in our study
model = changeRxnBounds(model,'RM00243',0,'l'); % assume glutamate dehydrogenase only in forward direction (making nh4) in our study

% block the thermodynamically infeasible loops that can use the low-energy
% bound in ATP as a high-energy bound

% (1) the massive PPI produced in tRNA synthesis can be converted to GTP in a
% loop that further support protein synthesis; this is an infeasible loop
% to bypass real energy demand; we prevent this loop
model = changeRxnBounds(model,'RCC0139',0,'l'); % although BRENDA supports reversible, a sig. reverse flux is not likely feasible and this is the setting in human model


% rescue the fitting of a few genes 
model = changeRxnBounds(model,'DMN0033',0.01,'u'); % to fit responsive gene F39B2.3, demand q2 has to be active; we assume minimal flux
model = changeRxnBounds(model,'EX00535',-0.01,'l'); % to fit responsive gene F42F12.3, demand tststerone has to be active; we assume minimal flux for hormone


% setup parameters
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
% Prepare flux thresholds (epsilon values)

load('input/epsilon_generic_withUptakes.mat'); % see walkthrough_generic.m for guidance on generating the epsilon values

% Set parameters
modelType = 2; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
speedMode = 1;

% gather stuffs to analyze mechanisms for 
predTbl = readtable('output/fluxTable.csv');
resp_depend = predTbl.rxns(predTbl.merged_PFD_bounded & ~predTbl.expression_PFD_bounded);
simi_depend = predTbl.rxns(predTbl.triple_PFD_bounded & ~predTbl.merged_PFD_bounded);
predTbl.related_responsiveness_constraints = repmat({'ND'}, size(predTbl,1),1);
predTbl.related_similarity_constraints = repmat({'ND'}, size(predTbl,1),1);



load('input/metabolic/categ_merged.mat');
branchTbl = readtable('input/metabolic/final_branchPoint_table.csv'); % must have 'mets', 'rxn1','rxn2', and 'maxCosine'
allMets = unique(branchTbl.mets);

sigFlux = 1e-5; 

% we use the disassimilation rate to constrain the bacteria waste
disAssim = 0.65;

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


% analyze the mechanisms dependent on responsiveness
allGenes = union(ExpCateg.responsive, ExpCateg.nonresponsive);


%% check mechanism of a gene 
mygene = 'atic-1';
targetRxn = 'TCM1285';
parforFlag = 0;
relMipGapTol = 1e-12;
verbose = 0;

ExpCateg_tmp = ExpCateg;
responsive = ExpCateg_tmp.responsive;
nonresponsive = ExpCateg_tmp.nonresponsive;
responsive(strcmp(responsive,mygene)) = [];
nonresponsive(strcmp(nonresponsive,mygene)) = [];
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
1e-5, 0, 0, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);
% loose the soft constraint
myCSM.MILP_PFD.b(myCSM.MILP.minLowInd) = myCSM.minLow * 1.05;
[minval, maxval] = FVA_MILP(myCSM.MILP_PFD, model_coupled, {targetRxn},parforFlag,relMipGapTol,verbose)
MILP_delta = myCSM.MILP_PFD;

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
1e-5, 0, 0, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);
% loose the soft constraint
myCSM.MILP_PFD.b(myCSM.MILP.minLowInd) = myCSM.minLow * 1.05;
[minval0, maxval0] = FVA_MILP(myCSM.MILP_PFD, model_coupled, {targetRxn},parforFlag,relMipGapTol,verbose)
MILP_ori = myCSM.MILP_PFD;

checkFor = -1; % 1 for lb and -1 for ub

MILPproblem_minFlux = MILP_delta;
FluxObj = find(ismember(model_coupled.rxns,targetRxn)); 
% create a new objective function
c = zeros(size(MILPproblem_minFlux.A,2),1);
c(FluxObj) = 1;
MILPproblem_minFlux.c = c;
% reverse direction (lb)
MILPproblem_minFlux.osense = checkFor;
% solve the MILP
solution_delta = autoTuneSolveMILP(MILPproblem_minFlux,1,relMipGapTol,targetRxn);   
        

MILPproblem_minFlux = MILP_ori;
FluxObj = find(ismember(model_coupled.rxns,targetRxn)); 
% create a new objective function
c = zeros(size(MILPproblem_minFlux.A,2),1);
c(FluxObj) = 1;
MILPproblem_minFlux.c = c;
% reverse direction (lb)
MILPproblem_minFlux.osense = checkFor;
% solve the MILP
solution_ori = autoTuneSolveMILP(MILPproblem_minFlux,1,relMipGapTol,targetRxn);   
        
% mechanism for weird prediction of TCM1285 by art-1 and likely others (eg
% mel-32): the low flux of these gene limits the biomass production; since
% biomass is coupled with bacteria uptake, with the minimal yield; the
% bacteria uptake is furtehr limited; thus, the available side is limited,
% that further resulted in failure of uptaking cit and icit from side to
% reverse the flux of the transporter. In fact, these biomass reactions
% should not be limited thus the low flux of art-1 or mel-32 is more of a
% false positive constraints (they should be responsive gene).

% however, it is very important to realize that such FVA calculation exams
% the extreme flux distribution that are almost physiologically impossible.
% (for example, in the case of atic-1, the extreme condition uses this
% reaction to make 10fthf, instead of consuming 10fthf to make purine);
% therefore, it is more about the theoraticaly minimum or maximum instead
% of an interpretable distribution; the OFD FVA provides more interpretable
% information.

% in sum, the unintuitive predictions are often associated with
% quantitative constraints, the most common one is the
% biomass-yield-uptake related reactions; minimizations that influence
% biomass capacity can indirectly influence the bacteria uptake the
% constriant system and in turn change the FVA of unintuitive reactions.
% this is can be an artifect prediction, but can also be real. for
% instance, it indicates given limited capacity of producing biomass and
% optimizing the use of diet to produce biomass, there is not enough
% reservior to support some other fluxes such as the reverse citrate
% transport which requires a large amout cit from side metabolite
% (unrealistic). so we just dont overinterpret this part of result (too
% quantitative-dependent).

%% check by metabolite
met = '10fthf[c]'; % 'accoa[c]' leu-L
ori = listRxn(model,solution_ori.full(1:length(model_coupled.rxns)),met); % myCSM_merged.OFD
delta = listRxn(model_coupled,solution_delta.full(1:length(model_coupled.rxns)),met); % myCSM_merged.OFD


%% check mechanism of a metabolite 
myMet = 'ivcoa[m]';
targetRxn = 'RM04746';

branchTbl_tmp = branchTbl;
branchTbl_tmp(~strcmp(branchTbl_tmp.mets,myMet),:) = []; % check to model LOO or LOI

branchTbl_tmp2 = branchTbl;
branchTbl_tmp2(:,:) = []; % check to model LOO or LOI

% remove the responsiveness
ExpCateg.responsive = {};
ExpCateg.nonresponsive = {};

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
= IMATplusplus_wiring_triple_inetgration_final(model_coupled,epsilon_f,epsilon_r, ExpCateg, branchTbl_tmp2, modelType,speedMode,...
0.05, 0, 1, 1, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);
[minval0, maxval0] = FVA_MILP(myCSM.MILP_PFD, model_coupled, {targetRxn},parforFlag,relMipGapTol,verbose)
MILP_ori =  myCSM.MILP_PFD;


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
= IMATplusplus_wiring_triple_inetgration_final(model_coupled,epsilon_f,epsilon_r, ExpCateg, branchTbl_tmp, modelType,speedMode,...
0.05, 0.05, 1, 1, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);
[minval, maxval] = FVA_MILP(myCSM.MILP_PFD, model_coupled, {targetRxn},parforFlag,relMipGapTol,verbose)
MILP_added =  myCSM.MILP_PFD;


checkFor = -1; % 1 for lb and -1 for ub

MILPproblem_minFlux = MILP_added;
FluxObj = find(ismember(model_coupled.rxns,targetRxn)); 
% create a new objective function
c = zeros(size(MILPproblem_minFlux.A,2),1);
c(FluxObj) = 1;
MILPproblem_minFlux.c = c;
% reverse direction (lb)
MILPproblem_minFlux.osense = checkFor;
% solve the MILP
solution_added = autoTuneSolveMILP(MILPproblem_minFlux,1,relMipGapTol,targetRxn)      

MILPproblem_minFlux = MILP_ori;
FluxObj = find(ismember(model_coupled.rxns,targetRxn)); 
% create a new objective function
c = zeros(size(MILPproblem_minFlux.A,2),1);
c(FluxObj) = 1;
MILPproblem_minFlux.c = c;
% reverse direction (lb)
MILPproblem_minFlux.osense = checkFor;
% solve the MILP
solution_ori = autoTuneSolveMILP(MILPproblem_minFlux,1,relMipGapTol,targetRxn)      


% flux fixation effect 
% the metabolite constriants fixed the solution space for corresponding reactions to
% around OFD range. This caused a cascade effect such that some other
% reations has to carry flux to support the fixed reaction.(For example,
% the fixed reaction flux could have gone up to allow the other reaction to
% be zero, but being fixed prevents this to happen). this is also related
% to biomass-bacteria coupling (fixed flux for biomass related reactions
% prevents biomass to go up, in turn the bacteria uptake; This caused
% nutrient starvation such as cys-L has to be synthesized because the
% dietary cys-L is not enough. 

% the fixed flux can cause contraints to biomass and bacteria uptake in
% turn in different ways. for example, the gtp metabolite pushes
% gluconeogensis use of gtp to be as close to atp--> gtp conversion as
% posible, therefore limiting the flux through both. This further limits
% the protein synthesis rate that needs gtp and in turn bacteria uptake given the biomass constraints
% so the actual mechanism can be very indirect, and it is even hard to say
% if it is reasonable or not. However, they anyways seem to constrian the
% flux space in a reasonable region even for PFD FVA. so it looks good.

% the same biomass-limiting mechanism also applies to ala-L; this effect
% can be understood as the coupling between fluxes disallow unlimited use
% of some cheap flux, uptaking more bacteria and wasting nutrient. This makes sense
% and is aligned with the ideal of selective wiring. More generally, the metabolite
% constraint usually caps a set of flux, which is similar to the total flux cap; this 
% set of capped flux indirectly forces the flux through other reactions, an example
% is g3p where sugar usage was limited such that aa is broken down.
% another similar example is that the ivcoa balance forced coupling between 
% bcaa deg (ivcoa produced and comsumed in BCAA deg). This limits the export
% flux of ivcoa from [m] to [c] to fuel the BCFA synthesis because this
% flux caused the imbalance between BCAA deg upstream and downstream rxn.
% Therefore, this export flux was minimized thus limiting the biomass and
% bacteria uptake flux. This limiting step forces the model to use energy
% efficient pathways such as mitochondria beta oxidation, which makes this
% pathway bounded. Allowing free ATP generation can reset this bounded
% rxns. 

% in sum, the metabolite balance limiting the global solution space by
% indirectly minimizing flux through uncoupled rxns (if an uncoupled flux
% is competeing with a coupled one); sometime when the coupling cannot be
% balanced, their fluxes will be minimized. This converged to the
% limitation of flux levels through key reactions in biomass production,
% which in turn limits bacteria uptake. The effect works like a total flux
% constriant on specific reactions to greatly narrow down their solution
% space. In general it is not an intuitive effect, rather, it is a
% quantitative prediction (flux competition, yeild constriant and bacteria
% uptake (nutrient availability)). An intuitive understanding of this is
% that to enforce couplings around metabolites, the upstream and downstream
% flux are forced to be euqal, which limits how much the one with higher
% capacity can go up. Such upper limits of flux levels usually synergy with
% the lower limits imposed by yield, high genes, and minLow constraints,
% tightening up the solution space. This tightenning up usually has a
% cascade effect when it's on key reactions in biomass production route,
% which limits the nutrient flux and biomass flux, in turn changes the
% solution space for other indirect reactions (like a long range
% competition). Such effects are unexpected (becuase it is not intuitive)
% but not unreasonable. So we will see how it goes. 



%% check by metabolite
load('myCSM_expression_only.mat')
%%
met = 'ivcoa[m]'; % 'accoa[c]' leu-L
exp = listRxn(model,solution_ori.full(1:length(model_coupled.rxns)),met); % myCSM_merged.OFD
LOI = listRxn(model_coupled,solution_added.full(1:length(model_coupled.rxns)),met); % myCSM_merged.OFD

%% more visualization and inspection
fluxTbl = table(model.rxns, myCSM_triple.OFD, printRxnFormula(model, model.rxns,0));
fluxTbl.type = regexprep(fluxTbl.Var1,'[0-9]+','');
sortedTbl = sortrows(fluxTbl, [4, 2]);



%% outdated codes and notes
% %% expression only and no integration control 
% load('input/metabolic/categ_expression_only.mat');
% 
% myCSM = struct(); % myCSM: my Context Specific Model
% [myCSM.OFD,...
% myCSM.PFD,...
% myCSM.N_highFit,...
% myCSM.N_zeroFit,...
% myCSM.minLow,...
% myCSM.minTotal,...
% myCSM.minTotal_OFD,...
% myCSM.MILP,...
% myCSM.MILP_PFD,...
% myCSM.HGenes,...
% myCSM.RLNames,...
% myCSM.OpenGene,...
% myCSM.latentRxn,...
% myCSM.Nfit_latent,...
% myCSM.wasteDW]...
% = IMATplusplus(model_coupled,epsilon_f,epsilon_r, ExpCateg, modelType,speedMode,...
% 1e-5, 1, 1, 0.05, 0.01,0.05,10);
% 
% myCSM_exp = myCSM;
% save('myCSM_exp.mat','myCSM_exp');
% % %% merged
% % load('input/metabolic/categ_DE_only.mat');
% % 
% % myCSM = struct(); % myCSM: my Context Specific Model
% % [myCSM.OFD,...
% % myCSM.PFD,...
% % myCSM.N_highFit,...
% % myCSM.N_zeroFit,...
% % myCSM.minLow,...
% % myCSM.minTotal,...
% % myCSM.minTotal_OFD,...
% % myCSM.MILP,...
% % myCSM.MILP_PFD,...
% % myCSM.HGenes,...
% % myCSM.RLNames,...
% % myCSM.OpenGene,...
% % myCSM.latentRxn,...
% % myCSM.Nfit_latent,...
% % myCSM.wasteDW]...
% % = IMATplusplus(model,epsilon_f,epsilon_r, ExpCateg, modelType,speedMode,...
% % 1e-5, 1, 1, 0.05, 0.01,0.05,10);
% % 
% % myCSM_DE = myCSM;
% %% compare results
% myCSM_exp
% %myCSM_DE
% myCSM_merged
% %%
% myCSM_exp.N_highFit/length(myCSM_exp.HGenes)
% myCSM_exp.N_zeroFit/length(myCSM_exp.RLNames)
% 
% % myCSM_DE.N_highFit/length(myCSM_DE.HGenes)
% % myCSM_DE.N_zeroFit/length(myCSM_DE.RLNames)
% 
% myCSM_merged.N_highFit/length(myCSM_merged.HGenes)
% myCSM_merged.N_zeroFit/length(myCSM_merged.RLNames)
% % the fitting is not greatly compromised after DE information integration
% %% inspect the conflicts 
% tmp = model.rxns(ismember(model.rxns,myCSM_merged.RLNames));
% unfitted_zero = tmp(myCSM_merged.OFD(ismember(model.rxns,myCSM_merged.RLNames))~=0);
% unfitted_high = setdiff(myCSM_merged.HGenes, myCSM_merged.OpenGene);
% %%
% DEtbl = readtable('./../../MetabolicLibrary/2_DE/output/DE_merged_clean_pcutoff_0.005_master_table_FDR2d0.1_FC1.5_FCtype_log2FoldChange_raw_ALL.csv');
% N_DE = tabulate(DEtbl.RNAi);
% ID_table = readtable('./../../MetabolicLibrary/1_QC_dataCleaning/outputs/RNAi_summary_final_dataset_modified.csv');
% [A,B] = ismember(N_DE(:,1),ID_table.RNAi_ID);
% lookupTbl = readtable('./input/IDtbl.csv');
% N_DE(A,1) = ID_table.WBID(B(A));
% N_DE(~A,2) = {'NA'};
% [A,B] = ismember(N_DE(:,1),lookupTbl.WormBase_Gene_ID);
% N_DE(A,1) = lookupTbl.ICELgene(B(A));
% N_DE = N_DE(A,:);
% 
% 
% %% check case-by-case
% mygene = 'fat-1';% polg-1
% load('input/metabolic/categ_expression_only.mat');
% if(any(strcmp(ExpCateg.high, mygene)))
%     fprintf('the category in expression is high\n');
% elseif (any(strcmp(ExpCateg.zero, mygene)))
%     fprintf('the category in expression is zero\n');
% elseif (any(strcmp(ExpCateg.low, mygene)))
%     fprintf('the category in expression is low\n');
% else
%     fprintf('the category in expression is dynamic\n');
% end
% load('input/metabolic/categ_DE_only.mat');
% if(any(strcmp(ExpCateg.high, mygene)))
%     fprintf('the category in DE is high (DE = %d)\n',N_DE{strcmp(N_DE(:,1),mygene),2});
% elseif (any(strcmp(ExpCateg.zero, mygene)))
%     fprintf('the category in DE is zero (DE = %d)\n',N_DE{strcmp(N_DE(:,1),mygene),2});
% else
%     fprintf('the category in DE is unknown\n');
% end
% load('input/metabolic/categ_merged.mat');
% if(any(strcmp(ExpCateg.responsive, mygene)))
%     fprintf('the category in merged is responsive\n');
% elseif (any(strcmp(ExpCateg.nonresponsive, mygene)))
%     fprintf('the category in merged is nonresponsive\n');
% else
%     fprintf('the category in merged is unknonw\n');
% end
% 
% %% resolve conflicts
% blockgenes = {};%'piki-1'
% tmp = deleteModelGenes(model,blockgenes);
% model_adj = model;
% model_adj.lb = tmp.lb;
% model_adj.ub = tmp.ub;
% %openGene = {'F09E5.3','F42F12.3','F42F12.4','acl-14'};
% openGene = {};
% 
% load('input/metabolic/categ_merged.mat');
% % refine some category assignment or constraints to reconcile conflicts
% 
% % (a) GPR conflicts:
% % (1) F58H1.3 methionine salvage pathway. 
% % related to RCC0010 that is a lumped reaction. the original reactions have
% % two step, first step by F58H1.3. This is a hard conflict where F58H1.3 is
% % nonresponsive but F42F12.4 is responsive. However importantly, in the
% % inspection of network conflicts, we found that other genes (all three
% % upstream reactions) are nonresponsive, so together supporting the pathway
% % to be inactive. It is possibe that F42F12.4 is a misannotation or this
% % pathway starts from 2kmb from bacteria. Anyways, the evidence supports
% % the pathway to be off. so we resolve F42F12.4 conflicts by removing F42F12.4
% % (2) polg-1 mito pol gamma: the gene is moderately expressed and unknown
% % about the RNAi efficiency (target not down based on DE). We consider it
% % as unknown 
% % (3) fntb-1: it is known to be an essential subunit of the fnt complex,
% % the RNAi is not responsive regardless of down reg of target for 2Fold. It
% % may indicate some biology or the enzyme level is still sufficient.
% % Any case, we know this reaction should be on. so we put fntb-1 to unknown
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'polg-1','fntb-1'});
% 
% % (b) network conflicts
% % (1) moco related - maybe non-transcriptionally buffered
% % one reaction supports no flux while one supports no/weak flux. One uncertain  for
% % reconciling the conflicts, we remove these from nonresponsive
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'moc-3','T27A3.6'});
% 
% % (2) serine syn related
% %: a high expressed gene (F26H9.5) conflited with a
% % no-DE moderately expressed gene (C31C9.2). in the manual inspection of
% % the pathway, the downstream gene Y62E10A.13 is also low. improtantly,
% % both low are confident low (target down), so we change the high expressed
% % gene
% ExpCateg.high = setdiff(ExpCateg.high, {'F26H9.5'});
% ExpCateg.dynamic = [ExpCateg.dynamic; {'F26H9.5'}];
% 
% % (3) F39B2.3
% % this gene cannot be open since DMN0033 is prohibited to avoid loop. we
% % skip this.
% 
% % (4) 5adtststerone: F42F12.3
% % we open the exchange for this as the met is not in the side nutrient pool
% model_adj = changeRxnBounds(model_adj,'EX00535',-1000,'l');
% 
% % (5) F58H1.3 methionine salvage pathway.
% % see notes above 
% ExpCateg.responsive = setdiff(ExpCateg.responsive, {'F42F12.4'});
% 
% % (6) glycogen breakdown
% % T22F3.3 is highly expressed but unknown for the responsiveness. conflict
% % with the downstream reaction(agl-1) that is nonresponsive. we noticed
% % that genes in glycogene synthesis is also nonresponsive. so we think we
% % should follow the DE evidance when DE evidance conflicts with expression 
% ExpCateg.high = setdiff(ExpCateg.high, {'T22F3.3'});
% ExpCateg.dynamic = [ExpCateg.dynamic; {'T22F3.3'}];
% 
% % (7) phosphatidylglycerol synthesis (acl-14) 
% % to open acl-14 reaction, phosphatidylglycerol is allowed to be taken up. 
% % it represents direct use of phosphatidylglycerol in bacterial
% % biomass.Since the downstream product of phosphatidylglycerol in the model
% % can only go to biomass, the free exchange of phosphatidylglycerol wouldnt
% % cause a major energy or carbon imbalance 
% model_adj = changeRxnBounds(model_adj,'EX18126',-1000,'l');
% 
% % (8) THF mediated nadph production #INTERESTING, related with purine#
% % the two key enzyme, dao-3 and K07E3.4 (homo of MTHFD1/2) are not
% % responsive. this conflicts with the high expresison of alh-3. however,
% % alh-3 is undetermined for responsiveness as its RNAi is short and target
% % is not down. so, as a rule of thumb, we accept the known DE evidence.
% ExpCateg.high = setdiff(ExpCateg.high, {'alh-3'});
% ExpCateg.dynamic = [ExpCateg.dynamic; {'alh-3'}];
% 
% % (9) folate biosynthesis 
% % the key gene in making formate, cat-4 is responsive, however, conflits 
% % with the two downstream branches, ptps-1 (nonresp) Y105C5B.3 (zero
% % expression). We noticed that ptps-1 is low expressed (dynamic) and the
% % target gene RNAi effect is unknown (likely not down). so as a thumb of 
% % rule, when a evidenced open is conflicted with uncertain down, we
% % consider the open 
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'ptps-1'});
% 
% % (10) coq metabolism #INTERESTING#
% % coq-5 is highly expressed (unknown for responsiveness) and coq-6 is both 
% % high and responsive. these two genes conflict with
% % the nonresponsive genes coq-2 and clk-1. so we choose to follow DE
% % evidence. This is an interesting funding that the coq synthsis may not be
% % high-flux due to the direct uptake from bacteria, or, the KD could be fully
% % compensated from the bacterial coq
% % as a general rule, we consider such cases as flux-carrying by removing
% % the non-responsive. This flux-minimization will make them carrying
% % minimal flux
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'coq-2','clk-1'});
% 
% % (11) Methylglyoxal detoxification 
% % the high expressed gene glod-4 is conflicting with nonresponsive Y17G7B.3
% % however, Y17G7B.3 expression is low and target is not likely down.
% % Considering the methylglyoxal detoxi is important, and responsiveness of
% % glod-4 is not known, we allow the flux in this pathway
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'Y17G7B.3'});
% 
% % (12) propinate shunt canonical 
% % the high expression mmcm-1 conflicted with nonresponsive others (mce-1,
% % cpt-2, mlcd-1). the mmcm-1 is also not responsive, but the target is not
% % sig down. interestingly, we noticed that the FC metric indicates the
% % target is down. Therefore, we remove it from high category (but not put
% % it back to nonresp to be consistent with general rule)
% ExpCateg.high = setdiff(ExpCateg.high, {'mmcm-1'});
% ExpCateg.dynamic = [ExpCateg.dynamic; {'mmcm-1'}];
% 
% % (13) inositol phosphate conversion 
% % responsive gene mtm-1 is conflicting with all others. however, these
% % genes are likely neuron specific so RNAi efficciency very low. Indeed, we
% % didnt see any target down in related genes (mtm-3, mtm-6,ppk-3, vps-34)
% % even for some highly expressed one (ppk-3), so we allow flux in this
% % pathway (and we know its likely a neuronal specific)
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'vps-34','ppk-3'});
% 
% % (14) de novo purine biosynthesis #INTERESTING#
% % the responsive (phenotype) gene pacs-1 is conflicted with nonresp genes
% % ppat-1 and F10F2.2 and atic-1. however, the target of ppat-1 and f10f2.2 are all up
% % in seq. clearly it is contaminated (because target lowly expressed). the
% % real RNAi efficiency is uncertain, but likely efficient. However, to turn
% % on the purine pathway, folate (nadph producing) cycle must be on, which
% % is conflicted with nonresponsive dao-3 and K07E3.4. In addition, the
% % downstream atic-1 is nonresponsive with target down. Therefore, together
% % we think both purine de novo and folate is off! *surprising!*
% % probably we should check pacs-1 clone!
% ExpCateg.responsive = setdiff(ExpCateg.responsive, {'pacs-1'});
% 
% %（15) sec-trna biosynthesis 
% % seld-1 (8DE) are responsive and high, but conflict with 
% % secs-1 (moderate) and pstk-1 (low) are nonresponsive (no DE at all).
% % considering the DE number of seld-1 is also low, we consider the sec-trna
% % syn to be off
% ExpCateg.responsive = setdiff(ExpCateg.responsive, {'seld-1'});
% ExpCateg.high = setdiff(ExpCateg.high, {'seld-1'});
% ExpCateg.dynamic = [ExpCateg.dynamic; {'seld-1'}];
% 
% % (16) sphinganine degradation 
% % spl-1 is highly responsive but conflict with nonresponsive sphk-1. sphk-1
% % is moderately expression and target is not significantly down (but down
% % for a little by FC). In such situation, we consider the harder evidence
% % based on highly responsive spl-1
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'sphk-1'});
% 
% % the followings reconcile the unfitted zero
% 
% % (17) glycolysis 
% % low DE pgk-1 conflicts with gpd-1/2/3/4. however, looks like gpd is
% % actually also nonresponsive 
% % gpd-4/1/3 is 'responsive' but looks actually the offtarget effect as in 
% % the RNAi of 1/2/3/4, all gpd1-4 are down, but no other DE in 1/3/4. 
% % gpd-2 has weak response (DE valid =5) but still quite suspicious. could
% % be some tissue specificity that we couldnt address.
% % we move all them to nonresponsive
% ExpCateg.responsive = setdiff(ExpCateg.responsive, {'gpd-1','gpd-2','gpd-3','gpd-4'});
% ExpCateg.nonresponsive = union(ExpCateg.nonresponsive, {'gpd-1','gpd-2','gpd-3','gpd-4'});
% ExpCateg.high = setdiff(ExpCateg.high, {'gpd-1','gpd-2','gpd-3','gpd-4'});
% ExpCateg.dynamic = union(ExpCateg.dynamic, {'gpd-1','gpd-2','gpd-3','gpd-4'});
% 
% % (18) histidine degradation
% % T12A2.1 is responsive that conflicts with non-responsive upstream
% % (Y51H4A.7) and downstream (cpin-1). however, T12A2.1 is a low DE RNAi (N
% % =2) with one very highly changed gene. as it is a weak evidence, we turn
% % off this pathway
% ExpCateg.responsive = setdiff(ExpCateg.responsive, {'T12A2.1'});
% ExpCateg.nonresponsive = union(ExpCateg.nonresponsive, {'T12A2.1'});
% 
% % (19) pentose phosphate: 2dr5p #interesting#
% % a hard conflict between F09E5.3 (clear responsive) and F07A11.5 (target
% % down and nonrespon) and Y43F4B.5 (target down and only one DE (although
% % look relevant). for this special case, we have to allow flux and we dont
% % know which nonresponsive gene is actually carrying flux. In addition, it
% % could be due to incomplete model as 2dr5p may be directly from bacteria.
% % so to be conservative, we release both nonresponse genes
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'F07A11.5','Y43F4B.5'});
% 
% % (20) CL systhesis #interesting# #against basic hypothesis# 
% % the last step gene crls-1 is not responsive (target down, no DE).
% % However, it is required for biomass production, so conflicting with many
% % high gene. The flux should be on. based on PMID:29251771 and 22174409p,
% % this gene has tissue specificity and lof causes mito mophorlogy change,
% % but looks consistent that the loss of its flux is not deterious in a
% % broad sense. this gene is against our basic hypothsis. 
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'crls-1'});
% 
% % (21) leucine degradation #interesting#
% % ech-5 is nonresponsive while leucine degradation is clearly on based on
% % other genes. the ech-5 is moderately expressed and target down was not
% % observed. so not sure if RNAi is successful. plus this gene seems to have
% % tissue specificity per wormbase. so we remove it.
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'ech-5'});
% 
% % (22) cept-2 (phosphotidyl chorline) #interesting# #against basic hypothesis# 
% % cept-2 is nonresponsive (high expression, only 2DE, target is down).
% % maybe tissue, maybe wrong annotation (ie F54D7.2 is similar to cept-2 and
% % has phenotype). but it is a hard conflict with our basic hypothesis 
% % it's likely the model is wrong per Fei's literature search
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'cept-2'});
% 
% % (23) idi-1 and fdps-1 #interesting# #against basic hypothesis# 
% % idi-1 (maybe also fdps-1) is essential for biomass production and also experimentally 
% % (PMID:15765206). this RNAi is not responsive and target is not down
% % (moderately expressed). we suspect the RNAi is failed (tissue or technical). 
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'idi-1','fdps-1'});
% 
% % (24) M01E11.1 and icmt-1 #interesting# #against basic hypothesis# 
% % icmt-1 is zero expressed and both icmt-1 and M01E11.1 (target down) are 
% % not responsive. this is likely the wrong model reconstruction. for
% % example, the protSfarCmet should not be in biomass_minor since they are
% % not essential component; or the methyl transfering reaction can also be
% % catalyzed by some other gene. To reconcile for now, we just allow the
% % flux
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'icmt-1','M01E11.1'});
% 
% % （25）cpt-2: cpt in mitochondria #interesting# #against basic hypothesis# 
% % cpt-2 is a confident nonresponsive gene (target down, no DE). however, it
% % is required for other confident high gene (dif-1) and high gene (cpt-1)
% % to carry flux. the cpt-2 annotation is suspicious as it is the only cpt
% % in mitochondria. it may be a candidate for followup. For now, we allow
% % its flux [per Fei, annotattion is wrong, cpt-2 should be only for long
% % chain fa]
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'cpt-2'});
% 
% % (26) glycan synthsis;  #interesting# #against basic hypothesis# 
% % K09E4.2 essential for biomass production, but confidantly low. will be a good
% % candidate for followup. likely model has some problem.
% % C08B11.8,tag-179,gly-20 essentual for biomass, nonresponsive but target down unknown
% % a set of four nonresponsive genes in this linear pathway may indicate
% % something. we should look at it. (but indeed some gene like algn-2 is
% % highly responsive)
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'K09E4.2','C08B11.8','tag-179','gly-20'});
% 
% % (27) glycogen synthesis #interesting#
% % T04A8.7 is not responsive conflicts with highly expressed gsy-1 and
% % gyg-1/2. to note, the gyg-1/2 and gsy-1 are also nonresponsive although
% % their target is not observed down (highly expressed gene). so we decided
% % to turn off this pathway 
% ExpCateg.high = setdiff(ExpCateg.high, {'gsy-1','gyg-1','gyg-2'});
% ExpCateg.dynamic = union(ExpCateg.dynamic, {'gsy-1','gyg-1','gyg-2'});
% 
% % (28) bile acid synthesis  #interesting# #against basic hypothesis# 
% % ZK892.4 and C24A3.4 : C24A3.4 (nonresp) C24A3.4(zero exp) acs-20 conflicted with the confident
% % high acs-20 and T08H10.1. since C24A3.4 is unknown for target down or
% % not, we allow active flux
% % hsd-1 hsd-2 hsd-3: all of them are zero expression, however this
% % reactions are required for making bile acid. this is a good example where
% % key assumption about expression level falls apart.
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'C24A3.4'});
% ExpCateg.zero = setdiff(ExpCateg.zero, {'hsd-1','hsd-2','hsd-3'});
% ExpCateg.dynamic = union(ExpCateg.dynamic, {'hsd-1','hsd-2','hsd-3'});
% 
% % (29) Y71H2AM.6 in coa biosyn
% % Y71H2AM.6 is conflicting with high confi responsive gene F25H9.6 and
% % pnk-1. since Y71H2AM.6 target down or not is unknown, we allow flux
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'Y71H2AM.6'});
% 
% % (30) gale-1 and ZK1058.3 
% % gale-1 is high conf non-resp (target down);ZK1058.3 is low confi nonresp 
% % (target down unknown). At least one of the gale-1 and ZK1058.3
% % needs to be on to support glycan sysn (via reactions of gly-20 and bre-4)
% % since gale-1 is high confidence nonresponsive, we keep it and release the
% % other one
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'ZK1058.3'});
% 
% % (31) F21D5.1 in chitin syn  #interesting# #against basic hypothesis# 
% % confident nonresponsive F21D5.1 (target down) conflicts with confident
% % responsive C36A4.4. cannot be resolved but since it is producing chitin,
% % we allow the flux
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'F21D5.1'});
% 
% % (32) daf-18 in inositol
% % confident nonresponsive daf-18 conflicts with low-conf responsive aap-1
% % (5 DE but target down). in addition, another related gene age-1 is also nonresponsive,
% % yet target is not observed down. so it is tricky and may be tissue
% % relevant. also, pail345p should be drained, which will resolve the
% % conflict. so, as for now, we just allow flux
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'daf-18'});
% 
% % (33) mito tRNA transferases #interesting# #against basic hypothesis# 
% % sars-2 target down uncertain
% % Y41D4A.6 target down, but no other DE
% % C39B5.6 target down uncertain.
% % if it is just one gene, then likely some error. but there are three
% % genes, all of them are mito tRNA syntase. it must be some biology there.
% % but to allow the model producing biomass, we allow flux
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'sars-2','Y41D4A.6','C39B5.6'});
% 
% % (34) ftn-1 and ftn-2
% % the ftn-2 is confidently nonresponsive and ftn-1 is zero expression. this
% % reaction is required for biomass prodcution. but it looks this reaction
% % can be spontanous. so we allow the flux
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'ftn-2'});
% 
% % special cases
% % (35) fat-1
% % fat-1 is considered high conf non-responsive since target is down and
% % only one other DE (fat-2). so it may be no response or weak response. we
% % found that the Rtotal pools (rtotal3) is dependent on fat-1 flux in the
% % model because of stoichemotric balance. so, it is possible that
% % the C20:4 dependent on fat-1 is not essential (so they should be removed
% % from the pool). note the annotation is quite conveniencing per PMID:
% % 11972048. therefore, we think the fat-1 originally should carry flux, but
% % likely the loss of flux can be fully buffered in terms of transcription.
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'fat-1'});
% 
% 
% % NOTE: we may consider recal epsilons because of the constrains and model change
% 
% 
% myCSM = struct(); % myCSM: my Context Specific Model
% [myCSM.OFD,...
% myCSM.PFD,...
% myCSM.N_highFit,...
% myCSM.N_zeroFit,...
% myCSM.minLow,...
% myCSM.minTotal,...
% myCSM.minTotal_OFD,...
% myCSM.MILP,...
% myCSM.MILP_PFD,...
% myCSM.HGenes,...
% myCSM.RLNames,...
% myCSM.OpenGene,...
% myCSM.latentRxn,...
% myCSM.Nfit_latent,...
% myCSM.wasteDW]...
% = IMATplusplus_eng(model_adj,epsilon_f,epsilon_r, ExpCateg, modelType,speedMode,...
% 1e-5, 1, 1, 0.05, 0.01,0.05,10,0,openGene);
% 
% myCSM_merged
% myCSM
% %% zero is compromised 
% tmp = model.rxns(ismember(model.rxns,myCSM_merged.RLNames));
% fitted_zero1 = tmp(myCSM_merged.OFD(ismember(model.rxns,myCSM_merged.RLNames))==0);
% tmp = model_adj.rxns(ismember(model_adj.rxns,myCSM.RLNames));
% fitted_zero2 = tmp(myCSM.OFD(ismember(model_adj.rxns,myCSM.RLNames))==0);
% 
% intersect(setdiff(fitted_zero1, fitted_zero2),myCSM.RLNames)
% %% high is compromised 
% intersect(setdiff(myCSM_merged.OpenGene, myCSM.OpenGene),myCSM.HGenes)
% 
% %% update the merged control 
% myCSM_merged = myCSM;
% tmp = model.rxns(ismember(model.rxns,myCSM_merged.RLNames));
% unfitted_zero = tmp(myCSM_merged.OFD(ismember(model.rxns,myCSM_merged.RLNames))~=0);
% unfitted_high = setdiff(myCSM_merged.HGenes, myCSM_merged.OpenGene);
% myCSM_merged.N_highFit/length(myCSM_merged.HGenes)
% myCSM_merged.N_zeroFit/length(myCSM_merged.RLNames)
% save('myCSM_merged.mat','myCSM_merged');
% 
% %% compare flux distribution
% flux_cmp = table();
% flux_cmp.rxns = model.rxns;
% flux_cmp.exp = myCSM_exp.OFD;
% flux_cmp.exp_DE = myCSM_merged.OFD;
% flux_cmp.formula = printRxnFormula(model, model.rxns,0);
% flux_cmp.deltaflux = abs(flux_cmp.exp - flux_cmp.exp_DE);
% flux_cmp.Properties.RowNames = flux_cmp.rxns;
% flux_cmp('BIO0010',:)
% flux_cmp('EXC0050',:)
% flux_cmp('RC00754',:)
% %% FLUX SOLUTION SPACE
% p = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(p)
%     % Define the par cluster
%     myCluster = parcluster('local');
%     myCluster.NumWorkers = 128;
%     saveProfile(myCluster);
%     parpool(2,'SpmdEnabled',false);% adjust according to your computing environment
% end
%     
% %% setup FVA inputs
% targetRxns = {'RMC0011'};
% parforFlag = 1;
% relMipGapTol = 1e-3;
% verbose = false;
% 
% myFVA_exp = struct(); % save FVA in a structure variable
% % run FVA by calling:
% [myFVA_exp.lb, myFVA_exp.ub] = FVA_MILP(myCSM_exp.MILP, model, targetRxns,parforFlag,relMipGapTol,verbose);
% myFVA_exp
% 
% myFVA_merged = struct(); % save FVA in a structure variable
% % run FVA by calling:
% [myFVA_merged.lb, myFVA_merged.ub] = FVA_MILP(myCSM_merged.MILP, model_adj, targetRxns,parforFlag,relMipGapTol,verbose);
% myFVA_merged
% 
% %% FVA on custom objectives 
% % total bcaa breakdown
% 
% %% we can also consider to further constrain the model in addition to total
% % flux, as total flux is biased by the large fluxes (also interesting to
% % first inspect the large fluxes)
% 
% %% also give the random objective sampling a try 
% % try gapsplit sampler 
% addpath gapsplit/matlab/
% 
% 
% mysample = gapsplit(myCSM_merged.MILP,1000,'secondaryFrac',0,'minval',minval, 'maxval',maxval);% ,'secondaryFrac',0,'vars',1:length(model.rxns),'debug',1
% %%
% mysample2 = gapsplit(myCSM_exp.MILP,100,'secondaryFrac',0);% ,'secondaryFrac',0,'vars',1:length(model.rxns),'debug',1
% save('mysample2.mat','mysample2');
% %%
% model = changeRxnBounds(model,'EXC0050',-1,'l');% free bacteria uptake for integration
% tmp = optimizeCbModel(model);
% model.lb(strcmp(model.rxns,'BIO0010')) = tmp.obj * 0.9;
% mysample4 = gapsplit(model,10000,'secondaryFrac',0);% ,'secondaryFrac',0,'vars',1:length(model.rxns),'debug',1
% save('oriModel_09biomass_flux_samples.mat','mysample4');
% %% pentophosphate flux looks a good example
% histogram(mysample.samples(:,find(strcmp(model.rxns,'RC02035'))))
% %%
% histogram(mysample2.samples(:,find(strcmp(model.rxns,'RC02035'))))
% %%
% load('oriModel_flux_samples.mat');
% mysample_ori = mysample3;
% load('merged_flux_samples.mat');
% mysample_merged = mysample;
% load('exp_flux_samples.mat');
% mysample_exp = mysample;
% load('oriModel_09biomass_flux_samples.mat');
% mysample_ori_bio = mysample4;
% %%
% % ppp RC02035 -0.1:0.025:1.1;
% % bcaa RM01214  -0.1:0.025:0.5;
% % shunt RM04432 0:0.025:0.8;
% target = 'RC01218';
% BinEdges = -1000:100:1000;
% histogram(mysample_ori.samples(:,find(strcmp(model.rxns,target))),'Normalization','probability', 'BinEdges',BinEdges)
% hold on
% histogram(mysample_ori_bio.samples(:,find(strcmp(model.rxns,target))),'Normalization','probability', 'BinEdges',BinEdges)
% histogram(mysample_exp.samples(:,find(strcmp(model.rxns,target))),'Normalization','probability', 'BinEdges',BinEdges)
% histogram(mysample_merged.samples(:,find(strcmp(model.rxns,target))),'Normalization','probability', 'BinEdges',BinEdges)
% hold off
% legend({'intact','intact+biomass','expression','merged'})
% xlabel('Flux');
% ylabel('Probability')
% title(target)
% plt = Plot(); % create a Plot object and grab the current figure
% plt.BoxDim = [5.2, 4.3];
% plt.LineWidth = 2;
% plt.FontSize = 15;
% plt.LegendLoc = 'NorthEast';
% plt.FontName = 'Arial';
% plt.export(['figures/flux_probability_',target,'.tiff']);
% %% data inspection
% % ppp RC02035 -0.1:0.025:1.1;
% % bcaa RM01214  -0.1:0.025:0.5;
% % shunt RM04432 0:0.025:0.8;
% 
% % RC01600 looks a biased sampling 
% tail = 0.01;
% target = 'RC01600';
% figure(1)
% histogram(cutextreme(mysample_ori.samples(:,find(strcmp(model.rxns,target))),tail),'Normalization','probability')
% title('intact')
% xlabel(['zero perct = ', num2str(sum(mysample_ori.samples(:,find(strcmp(model.rxns,target)))==0)/size(mysample_ori.samples,1))]);
% figure(2)
% histogram(cutextreme(mysample_ori_bio.samples(:,find(strcmp(model.rxns,target))),tail),'Normalization','probability')
% title('intact+biomass')
% xlabel(['zero perct = ', num2str(sum(mysample_ori_bio.samples(:,find(strcmp(model.rxns,target)))==0)/size(mysample_ori.samples,1))]);
% figure(3)
% histogram(cutextreme(mysample_exp.samples(:,find(strcmp(model.rxns,target))),tail),'Normalization','probability')
% title('expression')
% xlabel(['zero perct = ', num2str(sum(mysample_exp.samples(:,find(strcmp(model.rxns,target)))==0)/size(mysample_ori.samples,1))]);
% figure(4)
% histogram(cutextreme(mysample_merged.samples(:,find(strcmp(model.rxns,target))),tail),'Normalization','probability')
% title('merged')
% xlabel(['zero perct = ', num2str(sum(mysample_merged.samples(:,find(strcmp(model.rxns,target)))==0)/size(mysample_ori.samples,1))]);
% 
% %plt.export(['figures/flux_probability_',target,'.tiff']);
% %%
% met = 'leu-L[m]'; % 'accoa[c]'
% tbl1 = listRxn(model,myCSM_exp.OFD,met);
% tbl2 = listRxn(model_adj,myCSM_merged.OFD,met);
% 
% 
% %%
% met = 'g6p-B[c]'; % 'accoa[c]'
% tbl3 = listRxn(model,myCSM.OFD,met);
% %% 
% model = changeRxnBounds(model,'EXC0050',-1,'l');% free bacteria uptake for integration
% tmp = optimizeCbModel(model);
% model.lb(strcmp(model.rxns,'BIO0010')) = tmp.obj * 0.9;
% options.nFiles = 50;
% options.nPointsReturned = 10000;
% [modelSampling,samples,volume] = sampleCbModel(model,'tmp','ACHR');% ,'secondaryFrac',0,'vars',1:length(model.rxns),'debug',1
% save('achr_result.mat','modelSampling','samples','volume');
% %%
% % BIASED reaction RM01090
% tail = 0;
% target = 'BIO0010';
% tmp = [cutextreme(mysample_ori_bio.samples(:,find(strcmp(model.rxns,target))),tail);...
%     cutextreme(samples(find(strcmp(modelSampling.rxns,target)),:),tail)'];
% 
% h0 = histogram(tmp);
% figure(1)
% histogram(cutextreme(mysample_ori_bio.samples(:,find(strcmp(model.rxns,target))),tail),'Normalization','probability','BinEdges',h0.BinEdges)
% title('intact+biomass-Gapsplit')
% xlabel(['zero perct = ', num2str(sum(abs(mysample_ori_bio.samples(:,find(strcmp(model.rxns,target)))-0)<1e-5)/size(mysample_ori_bio.samples,1))]);
% figure(2)
% histogram(cutextreme(samples(find(strcmp(modelSampling.rxns,target)),:),tail),'Normalization','probability','BinEdges',h0.BinEdges)
% title('intact+biomas-ACHR')
% xlabel(['zero perct = ', num2str(sum(abs(samples(find(strcmp(modelSampling.rxns,target)),:)-0)<1e-5)/size(samples,2))]);

%% old notes
% modify the DE evidence to reconcile conflicts at GPR level
% (1) F58H1.3 methionine salvage pathway. 
% related to RCC0010 that is a lumped reaction. the original reactions have
% two step, first step by F58H1.3. This is a hard conflict where F58H1.3 is
% nonresponsive but F42F12.4 is responsive. It may related to incomplete
% network reconstruction (like missing a branch). We simply bypass it by
% allow this reaction to be on (put F58H1.3 to unknown).
% (2) polg-1 mito pol gamma: the gene is moderately expressed and unknown
% about the RNAi efficiency (target not down based on DE). We consider it
% as unknown 
% (3) fntb-1: it is known to be an essential subunit of the fnt complex,
% the RNAi is not responsive regardless of down reg of target for 2Fold. It
% may indicate some biology or the enzyme level is still sufficient.
% Any case, we know this reaction should be on. so we put fntb-1 to unknown
% ExpCateg.nonresponsive = setdiff(ExpCateg.nonresponsive, {'F58H1.3','polg-1','fntb-1'});


% we start to have a lot of GPR level conflict due to non-responsive
% conditions, such as pdha-1. we need to consider either extensively refine
% or soften the use (not assigning to zero category) of nonresponsive
% maybe we can consider force responsive genes to high gene (instead of
% requiring a high gene to be associated with at least one high reaction)
% -- this can be done by revising iMAT line 183-187 (meaning allow
% competition of high gene and low reaction fitting on the same reaction) -
% likely high gene will never win - so think about it
% ==> we focus on the conflict between high and zero, that are basically
% generated by conflicts between responsive and non-responsive genes. we
% impose responsive genes in high gene set anyways and check the conflicts
% in fitting for details 

% ==> using new criteria (conflict with zero rxns), there are not too many
% conflicts. often one zero gene caused multiple high problem. we just fix
% by curation. 
