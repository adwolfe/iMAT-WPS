function FinalBranchTbl = extract_DEsimilarity_constraints(conditionInfo)


%% parameters 

% set the significant cosine threshod; gene pairs with similarity higher
% than the threshold will be included in the analysis. This threshold is a
% rough heuristic estiamte based on the distributions of maximum cosines
% (see REWIRING paper for details).
sigCosine = 0.2;

%% load the tables
% load DE similarity matrix 
DEsim = readtable('./input/WPS/cosineSimilarity_FC_denoised_stationery_metabolic.csv','ReadRowNames',true);

% subset the DE sim matrix for only genes in the sample 
DEsim = DEsim(ismember(DEsim.Properties.RowNames, regexprep(conditionInfo.RNAiID,' ','_')),...
                ismember(DEsim.Properties.RowNames, regexprep(conditionInfo.RNAiID,' ','_')));

% load model
load('input/model/makeWormModel/iCEL1314_withUptakes.mat');
% load some annotation tables
lookupTbl = readtable('./input/model/IDtbl.csv');
% load the flux capacity defined by running FVA (this FVA code is provided
% in a2_1_run_integration.m
load('input/model/capacity_generic_withUptakes.mat');
% update the boundary -- so indirectly inreversible reactions are now noted
% in its boundary 
model.ub(~capacity_f) = min(model.ub(~capacity_f),zeros(sum(~capacity_f),1));
model.lb(~capacity_r) = max(model.lb(~capacity_r),zeros(sum(~capacity_r),1));

%% analyze branching points
% loop through all metabolites and look at all metabolites with DE similarity
% data. This is to identify all producing/consuming reaction pairs (P/C pair)
% in the network that are associated with WPS DE similarity data

% gether measured genes (responsive genes in the WPS dataset)
measuredGenes = conditionInfo.RNAi_WBID(ismember(regexprep(regexprep(conditionInfo.RNAiID,' ','_'),'\.','_'), DEsim.Properties.VariableNames));
measuredGenes = lookupTbl.ICELgene(ismember(lookupTbl.WormBase_Gene_ID, measuredGenes));
branchTbl = table(); % first inspect all possible pairs in all branches with data
mets = [];
rxn1 = [];
rxn2 = [];
type = [];
cosines = {};
maxCosine = [];

% we exclude the hub metabolites as they are by-products in many reactions
% or they are usually not associated with major flux wiring (too many
% reactions is associated with each other by hub metabolites which is not
% likely a realistic coupling mediated by the metabolite) (and they are
% hard to model with FBA as the constriant/obj becomes ambigous)

% hub metabolites are also often associated with a few genes that are
% associated with many reactions (such as hacd-1). in this case, if we
% equalize the production and consumption of too many reactions, it
% essentially does not inform about useful flux wiring. so the hubs should
% be removed 

% some inorganic metaboites that are not likely couple fluxes were also
% removed
metExcl = {'h2o[c]', 'h2o[m]','h2o[e]','h[c]','h[e]','h[m]','o2[c]', 'o2[m]', 'o2[e]'...
            'co2[c]', 'co2[m]', 'co2[e]','ppi[c]', 'ppi[m]', 'ppi[e]', 'pi[c]', 'pi[m]', 'pi[e]'...
            'nad[m]','nad[c]','nad[e]', 'nadh[m]', 'nadh[c]', 'nadh[e]',...
            'fad[m]', 'fad[c]', 'fad[e]', 'fadh2[m]', 'fadh2[c]', 'fadh2[e]',...
            'etfox[m]','etfrd[m]',...
            'nadp[m]', 'nadp[c]', 'nadp[e]', 'nadph[m]', 'nadph[c]', 'nadph[e]',...
            'atp[c]', 'atp[m]', 'atp[e]', 'adp[c]', 'adp[m]', 'adp[e]',...
            'coa[c]','coa[m]','coa[e]'};
% accoa was not excluded bc it has only a few P/C type correlation - likely
% meaningful

% NOTE: based on the criteria, we had to remove nadph (it indeed introduces stupid
% coupling for nadp[m]; so the coupling between PPP and FA syn is not 
% directly used as an input data for flux prediction

% generate the metabolite split table
RNAiGeneNameList = repmat({'NA'},size(DEsim,1),1);
[A B] = ismember(DEsim.Properties.VariableNames, regexprep(regexprep(conditionInfo.RNAiID,' ','_'),'\.','_'));
RNAiGeneNameList(A) = conditionInfo.RNAi_WBID(B(A));
[A B] = ismember(RNAiGeneNameList, lookupTbl.WormBase_Gene_ID);
RNAiGeneNameList(A) = lookupTbl.ICELgene(B(A));

for i = 1:length(model.mets)
    if (~strcmp(model.mets(i), metExcl))
        % all associated rxns
        myrxns = model.rxns(any(model.S(strcmp(model.mets,model.mets{i}),:),1));
        myInds = any(model.S(strcmp(model.mets,model.mets{i}),:),1);
        % keep reactions with data
        myRxnGeneMat = model.rxnGeneMat(myInds,ismember(model.genes,measuredGenes));
        myrxns = myrxns(any(myRxnGeneMat,2));
        % fprintf('calculating #%d: %s ... %d rxns\n',i, model.mets{i}, length(myrxns));
        % produce a table with all pairwise combination of these rxns
        for ii = 1:(length(myrxns)-1)
            for jj = (1+ii): length(myrxns)
                mets = [mets;model.mets(i)];
                rxn1 = [rxn1; myrxns(ii)];
                rxn2 = [rxn2; myrxns(jj)];
                if (  (model.lb(strcmp(model.rxns,myrxns(ii))) < 0 && model.ub(strcmp(model.rxns,myrxns(ii))) > 0) || ... 
                      (model.lb(strcmp(model.rxns,myrxns(jj))) < 0 && model.ub(strcmp(model.rxns,myrxns(jj))) > 0) ) % any reaction runs in both direction
                    type = [type; {'P/C'}]; % production - consumption connected
                elseif (  ( model.S(i,strcmp(model.rxns,myrxns(ii))) * model.lb(strcmp(model.rxns,myrxns(ii))) > 0 || ... 
                            model.S(i,strcmp(model.rxns,myrxns(ii))) * model.ub(strcmp(model.rxns,myrxns(ii))) > 0 ...
                            ) && ... % ii producing
                          ( model.S(i,strcmp(model.rxns,myrxns(jj))) * model.lb(strcmp(model.rxns,myrxns(jj))) < 0 || ... 
                            model.S(i,strcmp(model.rxns,myrxns(jj))) * model.ub(strcmp(model.rxns,myrxns(jj))) < 0 ...
                            )  ... % jj consuming;;
                               ) || ... 
                       (  ( model.S(i,strcmp(model.rxns,myrxns(jj))) * model.lb(strcmp(model.rxns,myrxns(jj))) > 0 || ... 
                            model.S(i,strcmp(model.rxns,myrxns(jj))) * model.ub(strcmp(model.rxns,myrxns(jj))) > 0 ...
                            ) && ... % jj producing
                          ( model.S(i,strcmp(model.rxns,myrxns(ii))) * model.lb(strcmp(model.rxns,myrxns(ii))) < 0 || ... 
                            model.S(i,strcmp(model.rxns,myrxns(ii))) * model.ub(strcmp(model.rxns,myrxns(ii))) < 0 ...
                            ) ... % ii consuming;;
                               )
                    % one consume and one produce 
                    type = [type; {'P/C'}];
                elseif (  ( model.S(i,strcmp(model.rxns,myrxns(ii))) * model.lb(strcmp(model.rxns,myrxns(ii))) > 0 || ... 
                            model.S(i,strcmp(model.rxns,myrxns(ii))) * model.ub(strcmp(model.rxns,myrxns(ii))) > 0 ...
                            ) && ... % ii producing
                          ( model.S(i,strcmp(model.rxns,myrxns(jj))) * model.lb(strcmp(model.rxns,myrxns(jj))) > 0 || ... 
                            model.S(i,strcmp(model.rxns,myrxns(jj))) * model.ub(strcmp(model.rxns,myrxns(jj))) > 0 ...
                            ) ... % jj producing;;
                                )
                    type = [type; {'P/P'}];
                elseif (  ( model.S(i,strcmp(model.rxns,myrxns(ii))) * model.lb(strcmp(model.rxns,myrxns(ii))) < 0 || ... 
                            model.S(i,strcmp(model.rxns,myrxns(ii))) * model.ub(strcmp(model.rxns,myrxns(ii))) < 0 ...
                            ) && ... % ii comsuming
                          ( model.S(i,strcmp(model.rxns,myrxns(jj))) * model.lb(strcmp(model.rxns,myrxns(jj))) < 0 || ... 
                            model.S(i,strcmp(model.rxns,myrxns(jj))) * model.ub(strcmp(model.rxns,myrxns(jj))) < 0 ...
                            ) ... % jj comsuming;;
                                )
                    type = [type; {'C/C'}];
                elseif model.lb(strcmp(model.rxns,myrxns(ii))) >= model.ub(strcmp(model.rxns,myrxns(ii))) ...
                        || ...
                        model.lb(strcmp(model.rxns,myrxns(jj))) >= model.ub(strcmp(model.rxns,myrxns(jj))) % one reaction is blocked
                    type = [type; {'B/X'}];
                else
                    error('there is a bug in the logic')
                end
                % make the reaction-reaction DE similarity metric by
                % taking the maximum 
                geneset1 = model.genes(model.rxnGeneMat(strcmp(model.rxns,myrxns(ii)),:) == 1);
                geneset2 = model.genes(model.rxnGeneMat(strcmp(model.rxns,myrxns(jj)),:) == 1);
                subMat = DEsim(ismember(RNAiGeneNameList, geneset1) , ismember(RNAiGeneNameList, geneset2));
                % mask the correlation mediated by the same gene/GPR (bc then
                % the cosine is not driven by common metabolic consequence,
                % instead, by simply the same gene) or overlaped GPR (then
                % the cosince could be driven by effects on the same
                % reaction perturbation)
                if any(ismember(subMat.Properties.VariableNames, regexprep(subMat.Properties.RowNames,'\.','_')))
                    % 1/2/2024: updated the code to filter out all
                    % similarities involving overlaped genes are masked
                    % with nan to be avoid confounded interpretation (i.e.,
                    % the similarity can be indicating the effect caused by
                    % perturbing the same reaction instead of a coupled
                    % pair of reactions).
                    subMat{:,ismember(subMat.Properties.VariableNames, regexprep(subMat.Properties.RowNames,'\.','_'))} = nan;
                    subMat{ismember(regexprep(subMat.Properties.RowNames,'\.','_'), subMat.Properties.VariableNames),:} = nan;
                end
                % take maximum and save data
                cosines = [cosines; {subMat}];
                maxCosine = [maxCosine; max(max(subMat{:,:}))]; % to enrich for signals, we consider use the maximum to represent any possible coupling
            end
        end
    end
end

% 8/9/23: fixed a bug in the code but the output table is unchanged 

branchTbl = table(mets, rxn1, rxn2, type, maxCosine, cosines);
if size(branchTbl,1) > 0
    branchTbl.formula1 = printRxnFormula(model,branchTbl.rxn1, 0);
    branchTbl.formula2 = printRxnFormula(model,branchTbl.rxn2, 0);
else
    branchTbl.formula1 = {};
    branchTbl.formula2 = {};
end


% next, produce the table of significant couplings for integration use

% when the two reactions are both producing or comsuming, it usually links
% to the common influence in downstream or upstream (like both cause energy
% deficiency). It is not related to flux wiring in this case; so we focus
% on P/C type 
branchTbl2 = branchTbl(branchTbl.maxCosine > sigCosine & strcmp(branchTbl.type,'P/C'),:);

% annotate the degree (so we know how many branches are there)
for i = 1:size(branchTbl2,1)
    branchTbl2.degree(i) = sum(model.S(strcmp(model.mets, branchTbl2.mets{i}),:) ~= 0);
end

% Quality check: check for metabolites with more than one pairs 
% if it is more than a pair of reaction, we should check if it is a big
% inter-correlated group; if not, it may not indicate flux wiring

% to check if they for inter-correlated group, we used a graph analysis
% method that checks the connectivity. 
if size(branchTbl2,1)>0
    npair = tabulate(branchTbl2.mets);
    complexMets = npair([npair{:,2}] > 1,1);
    noPassMet = {};
    for i = 1:length(complexMets)
        allPairs = branchTbl2(strcmp(branchTbl2.mets, complexMets{i}),:);
        % we consider it as a graph and check for connected graph components;
        % if there is only one, it is good
        s = allPairs.rxn1;
        t = allPairs.rxn2;
        G = graph(s,t);
        bins = conncomp(G);
        if length(unique(bins)) > 1
            noPassMet = [noPassMet,complexMets(i)];
        end
    end
    FinalBranchTbl = branchTbl2(ismember(branchTbl2.mets, [npair([npair{:,2}] == 1,1); ...
                                                           setdiff(complexMets, noPassMet)]),:);

    
    % the majority 90% of the complex metabolites were all connected (meaning
    % only one big inter-correlated group). it is a bit surprising (may be due
    % to the fact that many of these were driven by a single genes used in many
    % reactions).
    
    % Three metabolites were manually cleaned up in the full data, we keep the
    % same if the metabolite is in the sample set 
    if any(strcmp(branchTbl2.mets, 'akg[m]') & strcmp(branchTbl2.rxn1,'RMC0003'))
        FinalBranchTbl = [FinalBranchTbl; branchTbl2(strcmp(branchTbl2.mets, 'akg[m]') & strcmp(branchTbl2.rxn1,'RMC0003'),:)];
    end
    if any(strcmp(branchTbl2.mets, 'adn[c]'))
        FinalBranchTbl = [FinalBranchTbl; branchTbl2(strcmp(branchTbl2.mets, 'adn[c]'),:)];
    end
else
    FinalBranchTbl = branchTbl2;
end


