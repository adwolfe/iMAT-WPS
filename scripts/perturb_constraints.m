function [FVA_lb, FVA_ub] = perturb_constraints(MILPproblem_minFlux, model, targetRxn,parforFlag,relMipGapTol,verbose)
% Perform Flux Variability Analysis (FVA) given a MILProblem and target reactions. 
% The first nRxn variables in the MILProblem must be the fluxes of reactions in the input
% model.
%
% USAGE:
%
%    [lb, ub] = FVA_MILP(MILProblem, model, targetRxns, parforFlag)
%
% INPUTS:
%    MILPproblem_minFlux: the input MILP problem (COBRA MILP structure). 
%                       The MILP should be readily constrained for FVA
%                       calculation. For example, the total flux cap should
%                       be already set.
%    model:             input model (COBRA model structure)
%    targetRxns:        cell of target reactions to perform FVA on
%    parforFlag:        (0 or 1) whether to use parallel computing (parfor)
%    relMipGapTol:      the relative MIP gap to use in the FVA calculation
%
% OUTPUT:
%   ub:                 a vector of upper boundaries of queried reactions
%   lb:                 a vector of lower boundaries of queried reactions
%
% Additional Notice:    Please make sure the S matrix of the input MILP follows the structure of iMAT++ MILP. Some variables such as absolute flux proxy will be assumed to be at specifc positions, so errors will occur if the S matrix is not formed as standard iMAT++. 
%
% `Yilmaz et al. (2020). Final Tittle and journal.
% .. Author: - Xuhang Li, Mar 2020

if nargin < 4 || isempty(parforFlag)
    parforFlag = true;
end
if (nargin < 5) 
    relMipGapTol = 1e-12;
end
if (nargin < 6) 
    verbose = true;
end

fprintf('Start to perform the FVA...\n');
% Check if is running on gurobi solver
solverOK = changeCobraSolver('gurobi', 'MILP',0);
if ~solverOK
    fprintf('The solver parameter auto-tuning is not supported for current solver! Please use Gurobi for best performance!\n')
end

FluxObj = find(ismember(model.rxns,targetRxn)); 
MILPproblem_minFlux_ori = MILPproblem_minFlux;
%% analyze FVA
if parforFlag
    environment = getEnvironment();

    WaitMessage = parfor_wait(length(MILPproblem_minFlux.c), 'Waitbar', true); 

    parfor i = 1:length(MILPproblem_minFlux_ori.b)
        restoreEnvironment(environment);
        MILPproblem_minFlux = MILPproblem_minFlux_ori;

        % perturb one constraint
        MILPproblem_minFlux.A(i,:) = [];
        MILPproblem_minFlux.b(i) = [];
        MILPproblem_minFlux.csense(i) = [];

        % create a new objective function
        c = zeros(size(MILPproblem_minFlux.A,2),1);
        c(FluxObj) = 1;
        MILPproblem_minFlux.c = c;
        % reverse direction (lb)
        MILPproblem_minFlux.osense = 1;
        % solve the MILP
        solution = autoTuneSolveMILP(MILPproblem_minFlux,solverOK,relMipGapTol,targetRxn);       
        FVA_lb(i) = solution.obj;
        if verbose
            mylb = solution.obj;
        end
        % forward direction (ub)
        MILPproblem_minFlux.osense = -1;
        % solve the MILP
        solution = autoTuneSolveMILP(MILPproblem_minFlux,solverOK,relMipGapTol,targetRxn); 
        FVA_ub(i) = solution.obj;
        if verbose
            fprintf('Const.#%d: [%f, %f] \n',i,mylb, solution.obj);
        end
        WaitMessage.Send; 
    end
    WaitMessage.Destroy
else % same thing but in a for loop
    for i = 1:length(MILPproblem_minFlux_ori.b)
        
        MILPproblem_minFlux = MILPproblem_minFlux_ori;

        % perturb one constraint
        MILPproblem_minFlux.A(i,:) = [];
        MILPproblem_minFlux.b(i,:) = [];
        MILPproblem_minFlux.csense(i) = [];

        % create a new objective function
        c = zeros(size(MILPproblem_minFlux.A,2),1);
        c(FluxObj) = 1;
        MILPproblem_minFlux.c = c;
        MILPproblem_minFlux.osense = 1;
        % when parfor is not used, go with default (use up cores)
        % solve the MILP
        solution = autoTuneSolveMILP(MILPproblem_minFlux,solverOK,relMipGapTol,targetRxn,0);       
        FVA_lb(i) = solution.obj;
        if verbose
            mylb = solution.obj;
        end
        MILPproblem_minFlux.osense = -1;
        % solve the MILP
        solution = autoTuneSolveMILP(MILPproblem_minFlux,solverOK,relMipGapTol,targetRxn,0); 
        FVA_ub(i) = solution.obj;
        if verbose
            fprintf('Const.#%d: [%f, %f] \n',i,mylb, solution.obj);
        end
    end
end
end