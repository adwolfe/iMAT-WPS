function [OFD, PFD, N_highFit, N_zeroFit, minLow, minMetBalanceLoss, minTotal, minTotal_OFD, MILP, MILP_PFD, HGenes, RLNames, OpenGene, latentRxn, Nfit_latent, wasteDW, RHNames, branchMets, minimizedLowRxns] = IMATplusplus_wiring_triple_inetgration_final(model,epsilon_f,epsilon_r,ExpCateg, branchTbl, modelType, speedMode,minLowTol,minCouplingtol,doMinPFD, doLatent,latentCAP,storeProp,SideProp,ATPm,verbose) % bacCoef
% Uses the iMAT-WPS algorithm (`REF`) to find the optimal flux distribution 
% that best fit triple data including, absolute gene expression, WPS 
% responsiveness and WPS similarity. This function can be applied to other
% similar dataset and general models. 
%
% iMAT-WPS algorithm performs multi-step fitting to find a flux distribution
% that best agrees with categorized gene expression data based on expression 
% levels and WPS responsiveness, and simontanously the WPS similarity data. 
% It builds upon the previous iMAT++ algorithm that integrates gene
% expression data.
%
% USAGE:
%
%    OFD = IMATplusplus_wiring_triple_inetgration_final(model, epsilon_f, epsilon_r, ExpCatag, branchTbl, modelType)
%
% INPUTS:
%    model:             input model (COBRA model structure)
%    epsilon_f:         the epsilon sequence for the forward direction of all reactions in the model. Non-applicable reactions (i.e., for a irreversible,
%                       reverse-only reaction) should have a non-zero default to avoid numeric error.
%    epsilon_r:         the epsilon sequence for the reverse direction of all reactions in the model. Non-applicable reaction (i.e., for a irreversible,
%                       forward-only reaction) should have a non-zero default to avoid numeric error.
%    ExpCateg:          expression categories used for iMAT++ fitting. 
%                       The `ExpCateg` should be a structure variable with "high", "low", "dynamic" and
%                       "zero" four fields. See the walkthrough scripts (and our Github) on how to generate
%                       them from raw expression quantification data (i.e., TPM)
%    branchTbl:         table listing all producing/consuming (P/C)
%                       reaction pairs to be integrated with iMAT-WPS. Must
%                       have must have 'mets', 'rxn1','rxn2', and
%                       'maxCosine' columns. 
% OPTIONAL INPUTS:
%    modelType:         an integer to specify the input COBRA model. 
%                       1 == dual C. elegans tissue model
%                       2 == generic C. elegans model
%                       3 == other COBRA models. 
%                       iMAT++ has some custom actions (such as calculating
%                       bacteria waste) that are specific to C. elegans
%                       models. So, user should specify the correct model
%                       type to avoid errors from the custom actions.
%    speedMode:         (1, 2, or 3) to indicate which speed mode is used.
%                       We offer three modes of iMAT++ for the different 
%                       balance between computational speed and optimization
%                       stringency. The speed mode level 1 is the original
%                       iMAT++ configuration used in our tissue modeling.
%                       The level 3 gives fastest speed, so is recommanded 
%                       for all complex models, such as RECON models. 
%                       In level 3, the low reactions (dependent on rarely 
%                       and lowly expressed genes) will be constrained by 
%                       rigid boundaries after flux minimization. Additionally, 
%                       we release the MILP stringency to gain computational 
%                       speed. In level 2, we only release MILP strigency.
%                       But in general, the three modes give similar flux
%                       predictions.
%    minLowTol:         the relative tolerance for total flux of the 
%                       reactions dependent on lowly and rarely expressed 
%                       genes. This parameter tunes the strigency of 
%                       constraining low reaction fluxes to the minimal level.
%                       Default value (0.05, 5%) provides a balance between
%                       fitting of expression+responsiveness data and that
%                       of similarity data. 
%    minCouplingtol:    the relative tolerance for the objective function
%                       fitting of the similarity data. Default is 0.05 (5%)
%                       to offer a good balance between the fitting of
%                       minLow and that of flux coupling (similarity data).
%    doMinPFD:          whether to perform flux minimization of PFD 
%                       (after minimizing total flux for lowly and rarely expressed genes). 
%                       PFD flux minimization is required to obtain OFD, but  
%                       could be omitted if one only wants to get MILP fitting 
%                       result.
%    doLatent:          whether to perform recursive fitting of latent reactions
%    latentCAP:         the total flux cap for recursive fitting of latent
%                       reactions. The total flux will be capped at (1 +
%                       latentCAP)*OriginalTotalFlux; The default cap is
%                       0.05 (5%)
%    storeProp:         (dual-C. elegans tissue network only, use empty ('[]') for other models) the maximum allowed storage molecule uptake, as w/w percentage of bacteria uptake                  
%    SideProp:          (dual-C. elegans tissue network only, use empty ('[]') for other models) the maximum allowed side metabolites uptake, as w/w percentage of bacteria uptake                  
%    ATPm:              (C. elegans network only, use empty ('[]') for other models) the Non Growth Associated Maintenance (NGAM) used in fitting               
%    verbose:           (0 or 1) to show the MILP log or not
%
%
% OUTPUT:
%   OFD:                the Optimal Flux Distribution (OFD)
% OPTIONAL OUTPUTS:
%   PFD:                the primary flux distribution (PFD)
%   N_highFit:          the number of highly expressed genes fitted
%   N_zeroFit:          the number of rarely expressed genes fitted
%   minLow:             the minimal total flux of reactions dependent on
%                       rarely and lowly expressed genes
%   minTotal:           the minimal total flux of PFD (primary flux distribution)
%   minTotal_OFD:       the minimal total flux of OFD
%   MILP:               the MILP problem (in COBRA format) in the final
%                       flux minimization of OFD. This MILP can serve as a
%                       startpoint for any further analysis such FVA. This
%                       is referred to as Optimal Flux Model (OFM) in the
%                       paper.
%   MILP_PFD:           the MILP problem (in COBRA format) right before the
%                       first Flux minimization. This MILP is readily
%                       constrained for high/zero genes' fitting , and flux
%                       minimization of low reactions. This includes all
%                       data-driven constraints, so it is a conservative 
%                       but good startpoint for further custom analysis.
%                       This is referred to as Primary Flux Model (PFM) in
%                       the paper.
%   HGenes:             the list of highly expressed genes to be fitted
%   RLNames:            the list of reactions dependent on rarely expressed genes
%   OpenGene:           the list of fitted (carrying flux) highly expressed genes
%   latentRxn:          the list of all identified latent reactions to be fitted 
%   Nfit_latent:        the (total) number of latent reactions fitted
%   wasteDW:            (C. elegans network only) the percentage of wasted bacterial biomass (DW/DW)
%   RHNames:            List of "on" reactions in the fitting.
%   branchMets:         List of metabolites whose flux splits are fitted
%   minimizedLowRxns    List of low reactions whose flux are minimized.
%
% `Author. (202x). Final Tittle and journal.
%
% .. Author: - (COBRA implementation) Xuhang Li, Mar 202x
%% this is the main integration function
% apply default parameters
if (nargin < 6) % type of the input model 
    modelType = 1; % 1==>dual C.elegans; 2==>generic C. elegans; 3==>non-c.elegans
end
if (nargin < 7) 
    speedMode = 1; % by default, do normal IMAT++ (speed mode = 1)
end
if (nargin < 8)
    if speedMode ~= 3
        minLowTol = 1e-5; % default value of low reaction flux tolerance
    else
        minLowTol = 0; % not applicable to speed mode 3
    end
end
if (nargin < 9)
    minCouplingtol = 1e-5; % default value of flux coupling tolerance
end
if (nargin < 10)
    doMinPFD = true;
end
if (nargin < 11)
    doLatent = true;
end
if (nargin < 12) % the flux tolerance cap of total flux in the latent fitting. By default, we use 5%
    latentCAP = 0.05;
end
if (nargin < 13) % for C. elegans' dual-tissue model only
    if modelType == 1
        storeProp = 0.01;
    else
        storeProp = [];
    end
end
if (nargin < 14) % for C. elegans' dual-tissue model only
    if modelType == 1
        SideProp = 0.02;
    else
        SideProp = [];
    end
end
if (nargin < 15) % for C. elegans' model only
    if modelType ~= 3
        ATPm = 10;
    else
        ATPm = [];
    end
end
% if (nargin < 14) 
%     bacCoef = 1.01; % by default, don't show the MILP details
% end
if (nargin < 16) 
    verbose = 0; % by default, don't show the MILP details
end

% set global constant 
bacMW=966.28583751; % only will be used for C. elegans model
if speedMode > 1
    relMipGapTol = 0.001; % we release the MIPgap to 0.1% 
    % this released MipGap only applies to latent flux fitting step.
    % The optimization stringency of PFD is still kept.
    % users can release the MipGap for PFD manually if needed
else
    relMipGapTol = 1e-12;
end
% Check if is running on gurobi solver
solverOK = changeCobraSolver('gurobi', 'MILP',0);
if ~solverOK
    fprintf('The solver parameter auto-tuning is not supported for the current solver! Please use Gurobi for best performance!\n')
end
%% mapping the gene categories to reactions
fprintf('Start flux fitting... \n');
tic()
fprintf('Processing expression data... \n');
worm = model;
% process gene expression data
expression_gene=struct;
expression_gene.value = [3*ones(length(ExpCateg.high),1);2*ones(length(ExpCateg.dynamic),1);1*ones(length(ExpCateg.low),1);zeros(length(ExpCateg.zero),1)];
expression_gene.gene = [ExpCateg.high;ExpCateg.dynamic;ExpCateg.low;ExpCateg.zero];
[expressionRxns] = mapExpressionToReactions_DE_xl(worm, expression_gene, ExpCateg.responsive,ExpCateg.nonresponsive);
RHNames = worm.rxns(expressionRxns == 3);
RLNames = worm.rxns(expressionRxns == 0); % for the first step integration

% any non-zero reaction associated with the hard high genes should be high
% reaction (if it is originally moderate, it should be already high; if it
% is originally low bc of either low expression or nonresponsive, it should
% be overide; if it is originally zero, it should be notified by warning in
% the level mapping)
for i = 1: length(ExpCateg.responsive)
    myrxns = worm.rxns(any(worm.rxnGeneMat(:,strcmp(ExpCateg.responsive(i),worm.genes)),2),1);
    if any(expressionRxns(ismember(worm.rxns,myrxns)) > 0)
        nonZeroRxns = intersect(myrxns, worm.rxns(expressionRxns > 0));
        RHNames = union(RHNames, nonZeroRxns);
    end
end
% update the reaction labels - bug fixed 7/12/23
% these high reaction overides the original levels (rxn levels to to 3)
% this will remove any originally low reaction that has a responsive gene
expressionRxns(ismember(worm.rxns, RHNames)) = 3;


% for high gene (based on absolute expression), only genes 
% associated with high reactions are high genes! (not the input
% high category genes!)
HGenes = {};
maybeHigh = setdiff(union(ExpCateg.high, ExpCateg.responsive),ExpCateg.nonresponsive); % we exclude the nonresponsive genes from high set (regardless of their absolute expression)
for i = 1: length(maybeHigh)
    myrxns = worm.rxns(any(worm.rxnGeneMat(:,strcmp(maybeHigh(i),worm.genes)),2),1);
    if any(ismember(myrxns,RHNames))
        HGenes = [HGenes;maybeHigh(i)];
    end
end

% now, the responsive genes should be all in high gene set, QC here
if (any(~ismember(ExpCateg.responsive, HGenes)))
    error('the responsive genes are not correctly handled! check!')
end


%% begin to do iMAT++ integration
%% Step1: do MILP integration of highly expressed genes and rarely expressed genes
toc()
fprintf('fitting high genes with MILP... \n');
tic()
% setup model specifc constraints for C. elegans mdoel
if modelType == 1
    if storeProp(2) >=0 && SideProp(2) >=0 && ATPm>=0
        % adjust the ATPm according to the bacterial uptake. This requires 
        % manual tuning
        worm = changeRxnBounds(worm,'RCC0005_I',ATPm,'l');
        worm = changeRxnBounds(worm,'RCC0005_X',ATPm,'l');
        % change the side proportion
        worm.S(end-1, strcmp('EXC0050_L',worm.rxns)) = storeProp*bacMW*0.01;
        worm.S(end, strcmp('EXC0050_L',worm.rxns)) = SideProp*bacMW*0.01;
    else
        error('Please check your input parameters for side proportion, storage proportion, and ATPm!');
    end
elseif modelType == 2
    worm = changeRxnBounds(worm,'RCC0005',ATPm,'l');
    if storeProp(2) >=0 && SideProp(2) >=0 && ATPm>=0
        % adjust the ATPm according to the bacterial uptake. This requires 
        % manual tuning
        worm = changeRxnBounds(worm,'RCC0005',ATPm,'l');
        % change the side proportion
        worm.S(storeProp(1), strcmp('EXC0050',worm.rxns)) = storeProp(2)*bacMW*0.01;
        worm.S(SideProp(1), strcmp('EXC0050',worm.rxns)) = SideProp(2)*bacMW*0.01;
    else
        error('Please check your input parameters for side proportion, storage proportion, and ATPm!');
    end
elseif modelType == 3
    fprintf('Doing integration for a non-C. elegans model. Make sure you constrained the ATPm before running iMAT++! \n');
end
% perform the gene-centric MILP fitting of highly expressed genes and reactions dependent on rarely expressed genes
[~, ~,~,solution, MILProblem] = ExpressionIntegrationByMILP(worm, RHNames, RLNames,HGenes, epsilon_f, epsilon_r,[],[],verbose);
toc()
%% Step2: do minimization of low reactions
fprintf('Minimizing low flux... \n');
tic()
% convert the objective value in a solution to the model constraints
MILProblem =  solution2constraint(MILProblem,solution);
NfitInd = length(MILProblem.b);%save the index of this constraint for later use

% add new variables (absolute flux proxy) for flux minimization
MILProblem = addAbsFluxVariables(MILProblem, worm);

% set the initial solution (as the last step solution)
% only recommend to use when numerical issue is severe
% MILProblem.x0 = [solution.full;abs(solution.full(1:length(worm.rxns)))];

% minimize low flux (reactions dependent on lowly and rarely expressed
% genes)
% please note that, to have minimal number of variables in the MILP problem, 
% we do the flux minimization by the absolute flux proxy variables we just
% generated for total flux minimization

% create a new objective function
% Creating c vector (objective function)
c_v = zeros(size(worm.S,2),1);
c_yh1 = zeros(length(RHNames),1);
c_yl = zeros(length(RLNames),1);
c_yh2 = zeros(length(RHNames),1);
c_minFluxLowRxns = zeros(length(worm.rxns),1);
c_minFluxLowRxns(expressionRxns == 1 | expressionRxns == 0) = 1; % set the coeffi of "low" and "rare" reactions as 1
minimizedLowRxns = worm.rxns(expressionRxns == 1 | expressionRxns == 0);
c_gene = zeros(length(HGenes),1);
c = [c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns];
MILProblem.c = c;
% set the objective sense as minimize 
MILProblem.osense = 1;
% sovle the MILP problem
solution = solveCobraMILP_XL(MILProblem, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', verbose);
% we use customized solver interface to fine-tune the solver parameters for gurobi.
% However, user is able to replace `solveCobraMILP_XL` with original
% `solveCobraMILP` function in COBRAtoolbox. If using the original
% function, please specify the optimal, feasible and integer tolerance!

if solution.stat == 0 && solverOK % when failed to solve (infeasible), we start to tune solver parameter 
    % NOTE: SPECIFIC TO GUROBI SOLVER!
    fprintf('initial solving failed! Start the auto-tune...#1\n')
    gurobiParameters = struct();
    gurobiParameters.Presolve = 0;
    solution = solveCobraMILP_XL(MILProblem,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
    if solution.stat == 0
        fprintf('initial solving failed! Start the auto-tune...#2\n')
        gurobiParameters.NumericFocus = 3;
        solution = solveCobraMILP_XL(MILProblem,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
        if solution.stat == 0
            error('MILP solving failed! Please inspect the reason!');
        end
    end
end        
if solution.stat ~= 1
    if solution.stat == 3
        % pass with warning
        warning('Low fluxes are not fully eliminated! (solver time out! Change timeLimit if needed!)');
    else % unknown error
        error('MILP solving failed! Please inspect the reason!');
    end
end
toc()
fprintf('Minimizing low flux completed! \n');
minLow = solution.obj;

% add the minlow constraint
if speedMode == 3
    numTol = 1e-8; % too small tolerance causes numerical instability; too large lose the flux minimization's effect
    % constrain by current flux plus a numeric tolerance 
    MILProblem.ub(expressionRxns == 1 | expressionRxns == 0) = abs(solution.full(expressionRxns == 1 | expressionRxns == 0)) + numTol;
    MILProblem.lb(expressionRxns == 1 | expressionRxns == 0) = -abs(solution.full(expressionRxns == 1 | expressionRxns == 0)) - numTol;   
    % all absolute zero fluxes stay zero
    MILProblem.ub((expressionRxns == 1 | expressionRxns == 0) & (solution.full(1:length(worm.rxns)) == 0)) = 0;
    MILProblem.lb((expressionRxns == 1 | expressionRxns == 0) & (solution.full(1:length(worm.rxns)) == 0)) = 0;
    
    % then remove redundant variables (the zero reaction integer variables)
    c_v = zeros(size(worm.S,2),1);
    c_yh1 = zeros(length(RHNames),1);
    c_yl = ones(length(RLNames),1);
    c_yh2 = zeros(length(RHNames),1);
    c_minFluxLowRxns = zeros(length(worm.rxns),1);
    c_gene = zeros(length(HGenes),1);
    lowInd = logical([c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns]);
    MILProblem.A(:,lowInd) = [];% remove variables
    NfitZero = sum(abs(solution.full(expressionRxns == 0)) < 1e-8); % default int tol of iMAT++
    MILProblem.b(NfitInd) = MILProblem.b(NfitInd) - NfitZero; % update the Nfit constraints
    MILProblem.lb(lowInd) = [];% update boundary
    MILProblem.ub(lowInd) = [];% update boundary 
    MILProblem.vartype(lowInd) = [];% update variable type
    RLNames_ori = RLNames;
    RLNames = [];% update RLNames (should be empty because variables are deleted)
else
    % set a minFlux tolerance; 
    % to facilitate tuning, the tolerance now is a relative tolerence
    tol = minLowTol .* solution.obj; % it could be solver tolarence or some larger number to increase the numeric stability and allow flexible fitting of zero and low reactions
    solution.obj = solution.obj + tol;
    MILProblem = solution2constraint(MILProblem,solution);
end
minLowInd = length(MILProblem.b);
% MILProblem.x0 = solution.full;
%% Step 3: integrate the DE similarity constraint
% this is equivalent to doing optimization on top of the original PFD; we
% integrate the calculations in a single function to make it more
% streamlined; this will produce a new PFD that has DE similarity
% constraints

% algorithm:
% idea: The core idea of the algorithm is that the flux wiring can be
% mathmetically represented as the flux balance around the branching
% metabolite. The wired (ie DE similarity coupled) reactions should have
% balanced flux for this metabolite. To force the flux wiring around the 
% coupled reactions by cosine similarity, we simply extract the mass
% balance for corresponding metabolite around coupled reactions, use cosine
% to weight the importance of this balance, and optimize the best balance.

% mathmatically,

% for each metabolite i in the branchTble,
%   construct the mass-balance variable m(i)
%   m(i) = |C1(i) * V1(i) + C2(i) * V2(i) + ...|, V1, V2, V3,..., all
%                                                               reactions in the flux coupling module of this metabolite based on
%                                                               significant cosine (usually only two reactions). 
%   represented as the following in the minimization problem;
%   m(i) >= C1(i) * V1(i) + C2(i) * V2(i) + ...
%   m(i) >= -(C1(i) * V1(i) + C2(i) * V2(i) + ...)
%
%   minimize the total m(i),
%   min( w(1)*m(1) + w(2)*m(2) + w(3)*m(3) ...),
%       w(i) = max(maxCosine(pair1), maxCosine(pair2),...), all reaction
%                                                            pairs in the flux coupling module of m(i)
%b
%   the objective values represent the loss function of the unfitted flux
%   coupling in the branchTbl.

branchMets = unique(branchTbl.mets);

% add the met-mass-balance variables
% attach a 2*n(mets)-by-n(mets) sparse matrix to the original A at its
% lower right corner
A = sparse(size(MILProblem.A,1)+2*(length(branchMets)),size(MILProblem.A,2)+length(branchMets)); 
[m,n,s] = find(MILProblem.A);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end
% add variables 
%   m(i) >= C1(i) * V1(i) + C2(i) * V2(i) + ...
% => -C1(i) * V1(i) - C2(i) * V2(i) - ... + m(i) >= 0
for i = 1:length(branchMets)
    % first locate the reactions in the module 
    myrxns = unique([branchTbl.rxn1(strcmp(branchTbl.mets, branchMets{i})); branchTbl.rxn2(strcmp(branchTbl.mets, branchMets{i}))]);
    % set up the coefficient
    A(size(MILProblem.A,1)+i,ismember(worm.rxns, myrxns)) = -MILProblem.A(strcmp(worm.mets, branchMets{i}), ismember(worm.rxns, myrxns));
    A(size(MILProblem.A,1)+i,size(MILProblem.A,2)+i) = 1;
end
%   m(i) >= -(C1(i) * V1(i) + C2(i) * V2(i) + ...)
% => C1(i) * V1(i) + C2(i) * V2(i) + ... + m(i) >= 0
for i = 1:length(branchMets)
    % first locate the reactions in the module 
    myrxns = unique([branchTbl.rxn1(strcmp(branchTbl.mets, branchMets{i})); branchTbl.rxn2(strcmp(branchTbl.mets, branchMets{i}))]);
    % set up the coefficient
    A(size(MILProblem.A,1)+length(branchMets)+i,ismember(worm.rxns, myrxns)) = MILProblem.A(strcmp(worm.mets, branchMets{i}), ismember(worm.rxns, myrxns));
    A(size(MILProblem.A,1)+length(branchMets)+i,size(MILProblem.A,2)+i) = 1;
end

% update other fields
lb = [MILProblem.lb;zeros(length(branchMets),1)];
ub = [MILProblem.ub;1000*ones(length(branchMets),1)];
b = [MILProblem.b;zeros(2*length(branchMets),1)];
if size(MILProblem.csense,1) == 1 %for some model, the csense is a string instead of vector
    csense1(1:(2*length(branchMets))) = 'G';
    csense = [MILProblem.csense,csense1];
else %is a vector
    csense1(1:(2*length(branchMets)),1) = 'G';
    csense = [MILProblem.csense; csense1];
end
% modify vartype if it is an MILP
if isfield(MILProblem,'vartype')
    vartype1(1:length(branchMets),1) = 'C';
    vartype = [MILProblem.vartype;vartype1];
    MILProblem.vartype = vartype;
end
% generate function output
MILProblem.A = A;
MILProblem.b = b;
MILProblem.lb = lb;
MILProblem.ub = ub;
MILProblem.csense = csense;

% construct the objective function (we still expands all variable by
% groups)
c_v = zeros(size(worm.S,2),1);
c_yh1 = zeros(length(RHNames),1);
c_yl = zeros(length(RLNames),1);
c_yh2 = zeros(length(RHNames),1);
c_gene = zeros(length(HGenes),1);
c_minFluxLowRxns = zeros(length(worm.rxns),1);
c_fluxCoupling = zeros(length(branchMets),1);
for i = 1:length(branchMets)
    c_fluxCoupling(i) = max(branchTbl.maxCosine(strcmp(branchTbl.mets, branchMets{i})));
end
c = [c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns;c_fluxCoupling];
MILProblem.c = c;

% optimiaze the objective
fprintf('Best fitting the DE similarity... \n');
tic()
solution = solveCobraMILP_XL(MILProblem, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', verbose);
if solution.stat == 0 && solverOK% when failed to solve, we start to tune solver parameter #NOTE: SPECIFIC TO GUROBI SOLVER!%
    fprintf('initial solving failed! Start the auto-tune...#1\n')
    gurobiParameters = struct();
    gurobiParameters.Presolve = 0;
    solution = solveCobraMILP_XL(MILProblem,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
    if solution.stat == 0
        fprintf('initial solving failed! Start the auto-tune...#2\n')
        gurobiParameters.NumericFocus = 3;
        solution = solveCobraMILP_XL(MILProblem,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
        if solution.stat == 0
            error('MILP solving failed! Please inspect the reason!');
        end
    end
end  
if solution.stat ~= 1
    if solution.stat == 3
        % pass with warning
        warning('Total loss of flux coupling is not fully minimized! (solver time out! Change timeLimit if needed!)');
    else % unknown error
        error('MILP solving failed! Please inspect the reason!');
    end
end
minMetBalanceLoss = solution.obj;
tol = minCouplingtol .* solution.obj; % it could be solver tolarence or some larger number to increase the numeric stability and allow flexible fitting of DE similarity
solution.obj = solution.obj + tol;
MILProblem = solution2constraint(MILProblem,solution);
MILProblem.x0 = solution.full;
minMetBalanceLossInd = length(MILProblem.b);
fprintf('Best fitting the DE similarity completed! \n');
toc()
%% save check point
MILP_PFD = MILProblem;
MILP = MILProblem; % write the MILP output. If latent is done, this MILP will be updated.

%% Step 4: minimize total flux as needed

if doMinPFD
    % minimize total flux
    % again, to reduce computational burden, we use the same proxy
    % variables used in the low flux minimization
    % create a new objective function
    % Creating c vector (objective function)
    c_v = zeros(size(worm.S,2),1);
    c_yh1 = zeros(length(RHNames),1);
    c_yl = zeros(length(RLNames),1);
    c_yh2 = zeros(length(RHNames),1);
    c_minFluxLowRxns = ones(length(worm.rxns),1);
    c_gene = zeros(length(HGenes),1);
    c_fluxCoupling = zeros(length(branchMets),1);
    c = [c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns;c_fluxCoupling];
    MILProblem.c = c;
    fprintf('Minimizing total flux... \n');
    tic()
    solution = solveCobraMILP_XL(MILProblem, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', verbose);
    if solution.stat == 0 && solverOK% when failed to solve, we start to tune solver parameter #NOTE: SPECIFIC TO GUROBI SOLVER!%
        fprintf('initial solving failed! Start the auto-tune...#1\n')
        gurobiParameters = struct();
        gurobiParameters.Presolve = 0;
        solution = solveCobraMILP_XL(MILProblem,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
        if solution.stat == 0
            fprintf('initial solving failed! Start the auto-tune...#2\n')
            gurobiParameters.NumericFocus = 3;
            solution = solveCobraMILP_XL(MILProblem,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
            if solution.stat == 0
                error('MILP solving failed! Please inspect the reason!');
            end
        end
    end  
    if solution.stat ~= 1
        if solution.stat == 3
            % pass with warning
            warning('Total flux is not fully minimized! (solver time out! Change timeLimit if needed!)');
        else % unknown error
            error('MILP solving failed! Please inspect the reason!');
        end
    end
    minTotal = solution.obj;
    solution.obj = solution.obj*(1+latentCAP); % 1% minFlux constraint!
    MILProblem2 = solution2constraint(MILProblem,solution);
    MILProblem2.x0 = solution.full;
    MILP2 = MILProblem2; 
    MILP2.minTotalObj = c;
    minTotalInd = length(MILP2.b);


    % this is a special error-handling section. In rare cases, the highly
    % expressed genes are strongly conflicting with rarely expressed genes.
    % so, our minimization of low reaction fluxes will essentially kills
    % the fitting of most if not all high genes. Finally the minimal total flux
    % could be a very small number, even zero. We raise an error for such 
    % invalid fitting.
    if solution.obj <= minLowTol % minLowTol represents a small flux tolerance 
        error('All the highly expressed genes conflict with at least one rarely expressed gene. The flux fitting is not valid!')
    end
else
    minTotal = 0;
end
fprintf('Primary Flux Fitting completed! \n');
toc()
% make output of primary flux distribution and some optional outputs
FluxDistribution = solution.full(1:length(worm.rxns));
PFD = solution.full(1:length(worm.rxns));
if speedMode == 3 % the low fitting is not controlled by binary variable
    nameL = worm.rxns(ismember(worm.rxns,RLNames_ori));
    fluxL = solution.full(ismember(worm.rxns,RLNames_ori));
    ClosedLReaction = nameL(abs(fluxL)<1e-9);
else
    nameL = worm.rxns(ismember(worm.rxns,RLNames));
    yL = logical(solution.int((length(RHNames)+1):(length(RHNames)+length(RLNames))));
    ClosedLReaction = nameL(yL);
end
OpenGene = HGenes(logical(solution.int((2*length(RHNames)+length(RLNames))+1:end)));
N_highFit = length(OpenGene);
N_zeroFit = length(ClosedLReaction);
if doLatent
    %% step4. make the latent rxns fitting
    [FluxDistribution,latentRxn,Nfit_latent,minTotal_OFD,MILP] = fitLatentFluxes_triple_integration(MILP2, worm,PFD, HGenes,epsilon_f,epsilon_r,latentCAP,verbose,relMipGapTol);
    OFD = FluxDistribution;
else
    OFD = [];
    latentRxn = [];
    Nfit_latent = [];
    minTotal_OFD = [];
end
%% provide some major evaluation of the flux distribution (only for C. elegans model)
if modelType == 1
    fprintf('the bacterial uptake is: %f\n',FluxDistribution(strcmp(worm.rxns,'EXC0050_L'))); % bac
    bacWaste = ...
        FluxDistribution(strcmp(worm.rxns,{'EXC0051_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'protein_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0052_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'rna_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0053_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'dna_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0054_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'peptido_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0055_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'lps_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0056_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'lipa_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0057_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'outerlps_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0058_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'pe_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0059_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'pg_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0060_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'ps_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0063_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'clpn_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0065_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'glycogen_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0072_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'soluble_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0142_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'phospholipid_BAC_L[e]')});
    wasteDW = bacWaste / (-FluxDistribution(strcmp(worm.rxns,'EXC0050_L')) * MolMass(model.metFormulas{strcmp(model.mets,'BAC_L[e]')}));
    fprintf('the total waste of bulk bacteria is: %f\n',wasteDW);
elseif modelType == 2
    fprintf('the bacterial uptake is: %f\n',FluxDistribution(strcmp(worm.rxns,'EXC0050'))); % bac
    bacWaste = ...
        FluxDistribution(strcmp(worm.rxns,{'EXC0051'})) * MolMass(model.metFormulas{strcmp(model.mets,'protein_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0052'})) * MolMass(model.metFormulas{strcmp(model.mets,'rna_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0053'})) * MolMass(model.metFormulas{strcmp(model.mets,'dna_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0054'})) * MolMass(model.metFormulas{strcmp(model.mets,'peptido_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0055'})) * MolMass(model.metFormulas{strcmp(model.mets,'lps_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0056'})) * MolMass(model.metFormulas{strcmp(model.mets,'lipa_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0057'})) * MolMass(model.metFormulas{strcmp(model.mets,'outerlps_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0058'})) * MolMass(model.metFormulas{strcmp(model.mets,'pe_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0059'})) * MolMass(model.metFormulas{strcmp(model.mets,'pg_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0060'})) * MolMass(model.metFormulas{strcmp(model.mets,'ps_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0063'})) * MolMass(model.metFormulas{strcmp(model.mets,'clpn_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0065'})) * MolMass(model.metFormulas{strcmp(model.mets,'glycogen_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0072'})) * MolMass(model.metFormulas{strcmp(model.mets,'soluble_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0142'})) * MolMass(model.metFormulas{strcmp(model.mets,'phospholipid_BAC[e]')});
    wasteDW = bacWaste / (-FluxDistribution(strcmp(worm.rxns,'EXC0050')) * MolMass(model.metFormulas{strcmp(model.mets,'BAC[e]')}));
    fprintf('the total waste of bulk bacteria is: %f\n',wasteDW);
else
    wasteDW = [];
    fprintf('Please inspect the flux distribution manually for a non-C. elegans model\n');
end
%% fix the output for speed mode level 3
if speedMode == 3
    RLNames = RLNames_ori;% recover RLNames list for output
end
%% some QC of flux coupling fitting
MetFlux = 0;
for i = 1:length(branchMets)
    % first locate the reactions in the module 
    myrxns = unique([branchTbl.rxn1(strcmp(branchTbl.mets, branchMets{i})); branchTbl.rxn2(strcmp(branchTbl.mets, branchMets{i}))]);
    MetFlux = MetFlux + sum(abs(worm.S(strcmp(worm.mets, branchMets{i}), ismember(worm.rxns, myrxns)) .* FluxDistribution(ismember(worm.rxns, myrxns))'))/ 2;
end
fprintf('Flux coupling unbalanced mass %.2f/%.2f = %.2f%%! \n',minMetBalanceLoss, MetFlux, minMetBalanceLoss/MetFlux*100);
%% add some useful index 
MILP.minLowInd = minLowInd;
if doMinPFD
    MILP.minTotalInd = minTotalInd;
else
    MILP.minTotalInd = nan;
end
MILP.minMetBalanceLossInd = minMetBalanceLossInd;
end
