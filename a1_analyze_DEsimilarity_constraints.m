%% parameters 
sigCosine = 0.2;

%% load the tables
DEsim = readtable('./../../MetabolicLibrary/2_DE/output/cosineSimilarity_FC_denoised_stationery_metabolic.csv','ReadRowNames',true);
load('input/model/makeWormModel/iCEL1314_withUptakes.mat');
lookupTbl = readtable('./input/model/IDtbl.csv');
condInfo = readtable('./../../MetabolicLibrary/2_DE/output/RNAi_condition_metaInfo.csv');

load('input/model/capacity_generic_withUptakes.mat');
% update the boundary -- so indirectly inreversible reactions are now noted
% in its boundary 
model.ub(~capacity_f) = min(model.ub(~capacity_f),zeros(sum(~capacity_f),1));
model.lb(~capacity_r) = max(model.lb(~capacity_r),zeros(sum(~capacity_r),1));

%% analyze branching points
% loop through all metabolites and look at all metabolites with DE similarity
% data. how much DE similarity is related to a branching point?

measuredGenes = condInfo.RNAi_WBID(ismember(regexprep(regexprep(condInfo.RNAiID,' ','_'),'\.','_'), DEsim.Properties.VariableNames));
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
[A B] = ismember(DEsim.Properties.VariableNames, regexprep(regexprep(condInfo.RNAiID,' ','_'),'\.','_'));
RNAiGeneNameList(A) = condInfo.RNAi_WBID(B(A));
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
        fprintf('calculating #%d: %s ... %d rxns\n',i, model.mets{i}, length(myrxns));
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
                
                cosines = [cosines; {subMat}];
                maxCosine = [maxCosine; max(max(subMat{:,:}))]; % to enrich for signals, we consider use the maximum to represent any possible coupling
            end
        end
    end
end

% 8/9/23: fixed a bug in the code but the output table is unchanged 

branchTbl = table(mets, rxn1, rxn2, type, maxCosine, cosines);
branchTbl.formula1 = printRxnFormula(model,branchTbl.rxn1, 0);
branchTbl.formula2 = printRxnFormula(model,branchTbl.rxn2, 0);

writetable(branchTbl,'input/WPS/all_branchPoint_table.csv');

% next, produce the table of significant couplings for integration

% when the two reactions are both producing or comsuming, it usually links
% to the common influence in downstream or upstream (like both cause energy
% deficiency). It is not related to flux wiring in this case; so we focus
% on P/C type 
branchTbl2 = branchTbl(branchTbl.maxCosine > sigCosine & strcmp(branchTbl.type,'P/C'),:);

% annotate the degree (so we know how many branches are there)
for i = 1:size(branchTbl2,1)
    branchTbl2.degree(i) = sum(model.S(strcmp(model.mets, branchTbl2.mets{i}),:) ~= 0);
end

% check for metabolites with more than one pairs 
% if it is more than a pair of reaction, we should check if it is a big
% inter-correlated group; if not, it may not indicate flux wiring

% to check if they for inter-correlated group, we used a graph analysis
% method that checks the connectivity. 
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
%% we manully inspected these and clean up the suspicous ones
i = 3;
allPairs = branchTbl2(strcmp(branchTbl2.mets, noPassMet{i}),:);
% we consider it as a graph and check for Connected graph components;
% if there is only one, it is good
s = allPairs.rxn1;
t = allPairs.rxn2;
G = graph(s,t);
bins = conncomp(G);
plot(G)
% inspection result:
% 'udp[c] --> skip; low corr and doesnt seem reasonable
% 'akg[m]' --> keep the high correlation and likely meaningful one (glu -->
% akg)
% 'adn[c]' --> both are high corr and reasonable; likely to be in two
% tissues; we just model them together 

% we keep a set of manually curated good pairs
GoodPairs = [branchTbl2(strcmp(branchTbl2.mets, 'akg[m]') & strcmp(branchTbl2.rxn1,'RMC0003'),:);
             branchTbl2(strcmp(branchTbl2.mets, 'adn[c]'),:);];
FinalBranchTbl = [FinalBranchTbl; GoodPairs];  

% save table
writetable(FinalBranchTbl,'input/WPS/final_branchPoint_table.csv');

% final table contains information about 90 unique metabolite (95 with
% compartment) and 174 reaction pairs and 181 unique reaction. they are 
% related to 198 unique RNAi conditions and 182 unique genes
%% save figures for this
% the pie chart to show number of responsive gene covered 
% also some numbers 
N_met = length(unique(FinalBranchTbl.mets));
N_met_noComp = length(unique(regexprep(FinalBranchTbl.mets,'\[(m|c|e)\]$','')));
N_reaction = length(unique([FinalBranchTbl.rxn1; FinalBranchTbl.rxn2]));
N_reaction_covered = sum(any(model.S(ismember(model.mets, unique(FinalBranchTbl.mets)),:) ~= 0,1));

UsedCond = {};
for i = 1:size(FinalBranchTbl)
    UsedCond = [UsedCond; FinalBranchTbl.cosines{i}.Properties.RowNames; FinalBranchTbl.cosines{i}.Properties.VariableNames'];
end
UsedCond = unique(regexprep(UsedCond,'\.','_','all'));
UsedGenes = unique(RNAiGeneNameList(ismember(regexprep(DEsim.Properties.RowNames,'\.','_','all'), UsedCond)));
% pie chart 
figure;
text = {'contributing gene','not contributing gene'};
% create pie chart
pie([length(UsedGenes), length(setdiff(measuredGenes, UsedGenes))], false(2,1), text);
% save figure
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [3, 3];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.LegendLoc = 'NorthWest';
plt.Title = {['covered metabolites = ', num2str(N_met)];['covered rxns = ', num2str(N_reaction_covered)]};
plt.export(['figures/DEsimilarity_piechart.pdf']);



% label used pairs in the DE similarity heatmap 
% do this in R
% let's save the used gene pairs to the table (i.e. that gives the max cos)
for i = 1:size(FinalBranchTbl)
    mat = FinalBranchTbl.cosines{i};
    [maxVal, indx] = max(table2array(mat),[],'all','linear');
    [row_index, col_index] = ind2sub(size(mat), indx);
    if maxVal ~= mat{row_index, col_index}
        error('??')
    end

    FinalBranchTbl.pair_gene1(i) = mat.Properties.RowNames(row_index);
    FinalBranchTbl.pair_gene2(i) = mat.Properties.VariableNames(col_index);
end
% unique(regexprep([FinalBranchTbl.pair_gene1; FinalBranchTbl.pair_gene2],'\.','_','all'))
writetable(FinalBranchTbl, 'input/WPS/DEsim_table_integration_summary.csv');

%% notes about integration beyond simple splits
% one reasonable extension of the flux coupling around a bridging
% metabolite is to allow considering information in non-adjacent reactions
% in a linearly connected pathway to a split. For instance, we can consider
% for each metabolite, allow extension if the two adjacent reactions
% are fully coupled (this is given by weighted distance <= 1 or flux 
% coupling analysis [all rxn vs all rxn, we have the pipeline set up for 
% rewiring analysis]). To do this, we need to first compare how many new 
% constraints can be estabolished with the extension. This impelmentation
% is not technically simply and it may introduce a lot noises in fitting
% since this extension assumes indirect DE similarly can be explained by
% the network mechanism, whcih can be not true consideirng the regulation
% and incompleteness of the network. So, we opted to not complicate the
% algorithm by including more. (remember we already have constraints for
% 100 metabolites). 
%% notes in manual inspection of the table
% it may need curation before applied. eg focytC[m] couplings OXPHOS with focytC
% synthesis, which is obviously a bad coupling for flux wiring
% the followings are the apparently unreasonable couplings 
% (low cosine and unintuitive)
% check if they can be buffered or they had a negative impact on the flux
% integration
% (VERY BAD) 'ala-L[c]'	'RC03038'	'RCC0142'	'P/C'	0.306444282724250	1x1 table	'atp[c] + ala-L[c] + trnaala[c]  -> amp[c] + ppi[c] + alatrna[c] '	'2 atp[c] + h2o[c] + 2 cys-L[c] + 2 trdrd[c] + cpmp[c]  -> 4 h[c] + 2 amp[c] + 2 ppi[c] + 2 ala-L[c] + 2 trdox[c] + molybd[c] '	15
% (likely both producing) 'msa[m]'	'RM00908'	'RM01608'	'P/C'	0.217099278516508	1x1 table	'akg[m] + ala-B[m]  <=> glu-L[m] + msa[m] '	'nad[m] + 3hpp[m]  -> nadh[m] + h[m] + msa[m] '	4
% 'utp[c]'	'RC00416'	'BIO0003'	'P/C'	0.230503569066053	1x14 table	'h[c] + utp[c] + acgam1p[c]  <=> ppi[c] + uacgam[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	8
% 'glu-L[c]'	'RC07396'	'RC00573'	'P/C'	0.243073574138845	1x1 table	'glu-L[c] + 2kmb[c]  -> akg[c] + met-L[c] '	'atp[c] + h2o[c] + gln-L[c] + utp[c]  -> adp[c] + 2 h[c] + pi[c] + glu-L[c] + ctp[c] '	33
% 'glu-L[c]'	'RC00734'	'RC00573'	'P/C'	0.243073574138845	2x1 table	'akg[c] + tyr-L[c]  <=> glu-L[c] + 34hpp[c] '	'atp[c] + h2o[c] + gln-L[c] + utp[c]  -> adp[c] + 2 h[c] + pi[c] + glu-L[c] + ctp[c] '	33
% 'glu-L[c]'	'RC00694'	'RC00573'	'P/C'	0.243073574138845	2x1 table	'akg[c] + phe-L[c]  <=> glu-L[c] + phpyr[c] '	'atp[c] + h2o[c] + gln-L[c] + utp[c]  -> adp[c] + 2 h[c] + pi[c] + glu-L[c] + ctp[c] '	33
% 'glu-L[m]'	'RM01648'	'RM00258'	'P/C'	0.231206999741968	1x1 table	'akg[m] + 4abut[m]  -> sucsal[m] + glu-L[m] '	'akg[m] + ala-L[m]  <=> pyr[m] + glu-L[m] '	27
% 'glu-L[m]'	'RM04188'	'RM00258'	'P/C'	0.231206999741968	1x1 table	'glu-L[m] + 2mop[m]  -> akg[m] + 3aib[m] '	'akg[m] + ala-L[m]  <=> pyr[m] + glu-L[m] '	27
% 'glu-L[m]'	'RM00908'	'RM00258'	'P/C'	0.231206999741968	1x1 table	'akg[m] + ala-B[m]  <=> glu-L[m] + msa[m] '	'akg[m] + ala-L[m]  <=> pyr[m] + glu-L[m] '	27

% likely flux coupled but not major flux wiring (check the fitting result
% on these) (big flux --> small flux)
% 'amp[c]'	'RC00127'	'RC00185'	'P/C'	0.277577285608170	1x2 table	'atp[c] + amp[c]  <=> 2 adp[c] '	'atp[c] + adn[c]  -> adp[c] + h[c] + amp[c] '	69
% 'amp[c]'	'RC00127'	'RC05578'	'P/C'	0.255730102681089	1x1 table	'atp[c] + amp[c]  <=> 2 adp[c] '	'atp[c] + glu-L[c] + trnaglu[c]  -> amp[c] + ppi[c] + glutrna[c] '	69
% 'amp[c]'	'RC00127'	'RC03652'	'P/C'	0.227106926045622	1x1 table	'atp[c] + amp[c]  <=> 2 adp[c] '	'atp[c] + gln-L[c] + trnagln[c]  -> amp[c] + ppi[c] + glntrna[c] '	69
% 'amp[c]'	'RC00127'	'RC03658'	'P/C'	0.305276841143129	1x1 table	'atp[c] + amp[c]  <=> 2 adp[c] '	'atp[c] + lys-L[c] + trnalys[c]  -> amp[c] + ppi[c] + lystrna[c] '	69
% 'amp[c]'	'RC00127'	'RC00181'	'P/C'	0.252367952983668	1x1 table	'atp[c] + amp[c]  <=> 2 adp[c] '	'h[c] + amp[c] + h2o[c]  -> nh4[c] + imp[c] '	69
% 'amp[c]'	'RC05578'	'RC00181'	'P/C'	0.221284644491998	1x1 table	'atp[c] + glu-L[c] + trnaglu[c]  -> amp[c] + ppi[c] + glutrna[c] '	'h[c] + amp[c] + h2o[c]  -> nh4[c] + imp[c] '	69
% 'amp[c]'	'RC03658'	'RC00181'	'P/C'	0.287638620806091	1x1 table	'atp[c] + lys-L[c] + trnalys[c]  -> amp[c] + ppi[c] + lystrna[c] '	'h[c] + amp[c] + h2o[c]  -> nh4[c] + imp[c] '	69
% 'udpg[c]'	'RC01005'	'RC00289'	'P/C'	0.811145061649767	1x1 table	'dolp[c] + udpg[c]  -> udp[c] + dolglcp[c] '	'g1p[c] + h[c] + utp[c]  <=> ppi[c] + udpg[c] '	10
% 'amp[m]'	'RM03658'	'RM00127'	'P/C'	0.305276841143129	1x1 table	'atp[m] + lys-L[m] + trnalys[m]  -> ppi[m] + amp[m] + lystrna[m] '	'atp[m] + amp[m]  -> 2 adp[m] '	23

% interesting to check:
% (1) whether the g3p constraint eventually lead to the prediction of large
% flux in gluconeogeneiss 
% 'g3p[c]'	'RC01070'	'RC01066'	'P/C'	0.282639515581195	1x2 table	'fdp[c]  <=> dhap[c] + g3p[c] '	'2dr5p[c]  <=> g3p[c] + acald[c] '
% 'g3p[c]'	'RC01641'	'RC01827'	'P/C'	0.597500826723647	1x1 table	'g3p[c] + s7p[c]  <=> xu5p-D[c] + r5p[c] '	'g3p[c] + s7p[c]  <=> f6p[c] + e4p[c] '
% 'g3p[c]'	'RC01641'	'RC01066'	'P/C'	0.223896084077802	1x2 table	'g3p[c] + s7p[c]  <=> xu5p-D[c] + r5p[c] '	'2dr5p[c]  <=> g3p[c] + acald[c] '
% 'g3p[c]'	'RC01830'	'RC01827'	'P/C'	0.597500826723647	1x1 table	'f6p[c] + g3p[c]  <=> xu5p-D[c] + e4p[c] '	'g3p[c] + s7p[c]  <=> f6p[c] + e4p[c] '
% 'g3p[c]'	'RC01830'	'RC01066'	'P/C'	0.223896084077802	1x2 table	'f6p[c] + g3p[c]  <=> xu5p-D[c] + e4p[c] '	'2dr5p[c]  <=> g3p[c] + acald[c] '

% (2) if gluconeogenesis flux is correctly maintained and if RNA syn is
% disrupted 
% 'gtp[c]'	'RC00431'	'BIO0003'	'P/C'	0.429026643869237	2x14 table	'gtp[c] + oaa[c]  <=> pep[c] + gdp[c] + co2[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	15
% 'gtp[c]'	'RC00330'	'BIO0003'	'P/C'	0.322909515508412	1x14 table	'atp[c] + gdp[c]  <=> adp[c] + gtp[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	15
% 'oaa[c]'	'RC00431'	'RC00352'	'P/C'	0.211766335917045	2x2 table	'gtp[c] + oaa[c]  <=> pep[c] + gdp[c] + co2[c] '	'atp[c] + coa[c] + cit[c]  -> adp[c] + pi[c] + oaa[c] + accoa[c] '
% 'oaa[c]'	'RC00431'	'RC00355'	'P/C'	0.229297644288476	2x1 table	'gtp[c] + oaa[c]  <=> pep[c] + gdp[c] + co2[c] '	'akg[c] + asp-L[c]  <=> oaa[c] + glu-L[c] '
% 'gdp[c]'	'RC00431'	'RC00328'	'P/C'	0.348251665154161	2x1 table	'gtp[c] + oaa[c]  <=> pep[c] + gdp[c] + co2[c] '	'h2o[c] + gdp[c]  -> h[c] + pi[c] + gmp[c] '	18

% (3) pyruvate flux is coupled to alanine (although very weak) how would
% this influence the flux preditcion
% 'pyr[m]'	'RM00209'	'RM00258'	'P/C'	0.207368207147373	4x1 table	'pyr[m] + coa[m] + nad[m]  -> accoa[m] + co2[m] + nadh[m] '	'akg[m] + ala-L[m]  <=> pyr[m] + glu-L[m] '	12

% (4) accoa largely coupled  to mito beta and PDH; check how this influence
% the aa/shunt sourced accoa
% 'accoa[m]'	'RM00209'	'RM00351'	'P/C'	0.340519083595258	4x1 table	'pyr[m] + coa[m] + nad[m]  -> accoa[m] + co2[m] + nadh[m] '	'accoa[m] + h2o[m] + oaa[m]  -> coa[m] + h[m] + cit[m] '	35
% 'accoa[m]'	'RM00351'	'RMC0088'	'P/C'	0.257813234766866	1x4 table	'accoa[m] + h2o[m] + oaa[m]  -> coa[m] + h[m] + cit[m] '	'4 coa[m] + 4 nad[m] + 4 h2o[m] + 3 fad[m] + fa16p1n7coa[m]  -> 4 accoa[m] + 4 nadh[m] + 4 h[m] + 3 fadh2[m] + occoa[m] '	35
% and many others
% the branching points will be integrated in the iMAT as one more step of
% optimization;

% (5) cytosolic accoa largely coupled to histone acetylation; this is very
% surprisingly predicted by iMAT++ itself; double check if it changes 
% 'accoa[c]'	'RC01978'	'RC03552'	'P/C'	0.363271878610942	1x1 table	'h2o[c] + accoa[c] + aacoa[c]  -> h[c] + coa[c] + hmgcoa[c] '	'accoa[c] + Nlys[c]  <=> h[c] + coa[c] + Naclys[c] '	40
% 'accoa[c]'	'RCC0019'	'RC03552'	'P/C'	0.363704525616254	1x1 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'accoa[c] + Nlys[c]  <=> h[c] + coa[c] + Naclys[c] '	40
% 'accoa[c]'	'RC02058'	'RC03552'	'P/C'	0.270926143382104	2x1 table	'accoa[c] + gam6p[c]  -> h[c] + coa[c] + acgam6p[c] '	'accoa[c] + Nlys[c]  <=> h[c] + coa[c] + Naclys[c] '	40

% (6) cytC syn is coupled with ETC reduction of cytc; see if it has a negative
% impact on flux (should be correctly buffered)
% 'focytC[m]'	'RM00197'	'RM00081'	'P/C'	0.221484378466578	1x4 table	'lac-D[m] + 2 ficytC[m]  -> pyr[m] + 2 h[m] + 2 focytC[m] '	'8 h[m] + 4 focytC[m] + o2[m]  -> 4 h[c] + 2 h2o[m] + 4 ficytC[m] '	7
% 'focytC[m]'	'RM02161'	'RM00081'	'P/C'	0.701882771479615	8x4 table	'2 h[m] + qh2[m] + 2 ficytC[m]  -> 4 h[c] + q[m] + 2 focytC[m] '	'8 h[m] + 4 focytC[m] + o2[m]  -> 4 h[c] + 2 h2o[m] + 4 ficytC[m] '	7
% 'focytC[m]'	'RM00081'	'RM02480'	'P/C'	0.315962093857402	4x2 table	'8 h[m] + 4 focytC[m] + o2[m]  -> 4 h[c] + 2 h2o[m] + 4 ficytC[m] '	'apocytc[m] + pheme[m]  -> focytC[m] '	7

% (7) methionine influx: the weak corr between mtrr-1 and sams-1 adds
% coupling between folate-met influx and met/sam cycle. This may strongly
% change the flux to cyclic pattern that will be problemmatic; see how it
% will go! (to fix it, we may have to identify all no-measure pairs for all
% metabolites considered here, and always include these no-measure pairs in
% the total metabolite flux optimalziation; this will add the dietary met
% flux back into the influxes considered for met-L) But first see if it can
% be buffered.
% 'met-L[c]'	'RC00946'	'RC00177'	'P/C'	0.293263547288988	3x4 table	'hcys-L[c] + 5mthf[c]  -> h[c] + met-L[c] + thf[c] '	'atp[c] + h2o[c] + met-L[c]  -> pi[c] + ppi[c] + amet[c] '	8

% (8) ketone supply of hmgcoa (see if this will influence flux prediction)
% 'hmgcoa[c]'	'RC01978'	'RC02082'	'P/C'	0.203254558460800	1x1 table	'h2o[c] + accoa[c] + aacoa[c]  -> h[c] + coa[c] + hmgcoa[c] '	'2 h[c] + 2 nadph[c] + hmgcoa[c]  -> 2 nadp[c] + coa[c] + mev-R[c] '	3

% (9) are the many weak couplings around nucleic acid and DNA/RNA syn meaningful?
% 'ctp[c]'	'RC01799'	'RC00570'	'P/C'	0.259224683868232	1x1 table	'h[c] + pa_pl[c] + ctp[c]  -> ppi[c] + cdpdag[c] '	'atp[c] + cdp[c]  <=> adp[c] + ctp[c] '	10
% 'ctp[c]'	'RC00571'	'BIO0003'	'P/C'	0.240322150587729	1x14 table	'atp[c] + nh4[c] + utp[c]  -> adp[c] + 2 h[c] + pi[c] + ctp[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	10
% 'ctp[c]'	'RC00573'	'BIO0003'	'P/C'	0.240322150587729	1x14 table	'atp[c] + h2o[c] + gln-L[c] + utp[c]  -> adp[c] + 2 h[c] + pi[c] + glu-L[c] + ctp[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	10
% 'ctp[c]'	'RC00570'	'BIO0003'	'P/C'	0.322909515508412	1x14 table	'atp[c] + cdp[c]  <=> adp[c] + ctp[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	10
% 'utp[c]'	'RC00156'	'BIO0003'	'P/C'	0.322909515508412	1x14 table	'atp[c] + udp[c]  <=> adp[c] + utp[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	8
% 'utp[c]'	'RC00416'	'BIO0003'	'P/C'	0.230503569066053	1x14 table	'h[c] + utp[c] + acgam1p[c]  <=> ppi[c] + uacgam[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	8
% 'dctp[c]'	'RC02326'	'BIO0004'	'P/C'	0.309397287122737	1x3 table	'atp[c] + dcdp[c]  -> adp[c] + dctp[c] '	'0.368 atp[c] + 0.368 h2o[c] + 0.325 datp[c] + 0.175 dgtp[c] + 0.175 dctp[c] + 0.325 dttp[c]  -> 0.368 adp[c] + 0.368 h[c] + 0.368 pi[c] + ppi[c] + DNA[c] '	2
% 'dttp[c]'	'RC02093'	'BIO0004'	'P/C'	0.309397287122737	1x3 table	'atp[c] + dtdp[c]  <=> adp[c] + dttp[c] '	'0.368 atp[c] + 0.368 h2o[c] + 0.325 datp[c] + 0.175 dgtp[c] + 0.175 dctp[c] + 0.325 dttp[c]  -> 0.368 adp[c] + 0.368 h[c] + 0.368 pi[c] + ppi[c] + DNA[c] '	3

% (10) dag: not likely correct but see if it is buffered or actually
% meaningful
% 'dag_pl[c]'	'RC02239'	'RC02057'	'P/C'	0.221548350905133	2x1 table	'h2o[c] + pa_pl[c]  -> pi[c] + dag_pl[c] '	'dag_pl[c] + cdpea[c]  <=> h[c] + cmp[c] + pe[c] '	7
% 'dag_pl[c]'	'RC02057'	'RC03435'	'P/C'	0.219641996542157	1x3 table	'dag_pl[c] + cdpea[c]  <=> h[c] + cmp[c] + pe[c] '	'h2o[c] + pail45p[c]  -> h[c] + dag_pl[c] + mi145p[c] '	7
% 'dag_pl[c]'	'RC02057'	'RC03332'	'P/C'	0.219641996542157	1x3 table	'dag_pl[c] + cdpea[c]  <=> h[c] + cmp[c] + pe[c] '	'h2o[c] + pail[c]  -> h[c] + dag_pl[c] + mi1p-D[c] '	7

% (11) the source of ethamp; coupling is likely wrong but interesting to
% see how the model reacts to it
% 'ethamp[c]'	'RC02037'	'RC02464'	'P/C'	0.318065858870702	1x1 table	'amet[c] + ethamp[c]  -> h[c] + ahcys[c] + methamp[c] '	'sph1p_ce[c]  -> ethamp[c] + m15hxdcal[c] '	7

% (12) coupling between peroxi oxi and the FA syn is captured. see how
% model reacts
% 'pmtcoa[c]'	'RC07758'	'RCC0126'	'P/C'	0.399842366054898	2x4 table	'h[c] + pmtcoa[c] + malcoa[c]  -> co2[c] + coa[c] + 3oodcoa[c] '	'h2o[c] + nad[c] + coa[c] + o2[c] + stcoa[c]  -> h[c] + nadh[c] + accoa[c] + h2o2[c] + pmtcoa[c] '	11
% 'pmtcoa[c]'	'RCC0165'	'RCC0126'	'P/C'	0.692051466068343	1x4 table	'o2[c] + pmtcoa[c]  -> h2o2[c] + hdd2coa[c] '	'h2o[c] + nad[c] + coa[c] + o2[c] + stcoa[c]  -> h[c] + nadh[c] + accoa[c] + h2o2[c] + pmtcoa[c] '	11
% 'stcoa[c]'	'RCC0124'	'RCC0125'	'P/C'	0.412245389710732	4x4 table	'9 h[c] + 6 nadph[c] + 3 malcoa[c] + stcoa[c]  -> 3 h2o[c] + 3 co2[c] + 6 nadp[c] + 3 coa[c] + lgnccoa[c] '	'3 h2o[c] + 3 nad[c] + 3 coa[c] + 3 o2[c] + lgnccoa[c]  -> 3 h[c] + 3 nadh[c] + 3 accoa[c] + 3 h2o2[c] + stcoa[c] '	9
% 'eicostetcoa[c]'	'RCC0063'	'RCC0096'	'P/C'	0.420804806043524	5x4 table	'3 h[c] + 2 nadph[c] + malcoa[c] + strdnccoa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + eicostetcoa[c] '	'2 h2o[c] + 2 nad[c] + 2 coa[c] + 2 o2[c] + eicostetcoa[c]  -> 2 h[c] + 2 nadh[c] + 2 accoa[c] + 2 h2o2[c] + fa16p4n3coa[c] '	7
% 'dlnlcgcoa[c]'	'RCC0064'	'RCC0093'	'P/C'	0.420804806043524	5x4 table	'3 h[c] + 2 nadph[c] + malcoa[c] + lnlncgcoa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + dlnlcgcoa[c] '	'2 h2o[c] + 2 nad[c] + 2 coa[c] + 2 o2[c] + dlnlcgcoa[c]  -> 2 h[c] + 2 nadh[c] + 2 accoa[c] + 2 h2o2[c] + fa16p3n6coa[c] '	7
% 'fa16p1n7coa[c]'	'RCC0069'	'RCC0070'	'P/C'	0.493184449254954	1x13 table	'o2[c] + pmtcoa[c] + 2 focytb[c]  -> 2 h2o[c] + 2 ficytb[c] + fa16p1n7coa[c] '	'3 h[c] + 2 nadph[c] + malcoa[c] + fa16p1n7coa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + vacccoa[c] '	8
% 'fa16p1n7coa[c]'	'RCC0070'	'RCC0089'	'P/C'	0.420804806043524	13x4 table	'3 h[c] + 2 nadph[c] + malcoa[c] + fa16p1n7coa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + vacccoa[c] '	'h2o[c] + nad[c] + coa[c] + o2[c] + odecoa[c]  -> h[c] + nadh[c] + accoa[c] + h2o2[c] + fa16p1n7coa[c] '	8
% 'vacccoa[c]'	'RCC0070'	'RCC0090'	'P/C'	0.420804806043524	13x4 table	'3 h[c] + 2 nadph[c] + malcoa[c] + fa16p1n7coa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + vacccoa[c] '	'5 h2o[c] + 5 nad[c] + 5 coa[c] + 4 o2[c] + vacccoa[c]  -> 5 h[c] + 5 nadh[c] + 5 accoa[c] + 4 h2o2[c] + occoa[c] '	6
% 'lgnccoa[c]'	'RCC0124'	'RCC0125'	'P/C'	0.412245389710732	4x4 table	'9 h[c] + 6 nadph[c] + 3 malcoa[c] + stcoa[c]  -> 3 h2o[c] + 3 co2[c] + 6 nadp[c] + 3 coa[c] + lgnccoa[c] '	'3 h2o[c] + 3 nad[c] + 3 coa[c] + 3 o2[c] + lgnccoa[c]  -> 3 h[c] + 3 nadh[c] + 3 accoa[c] + 3 h2o2[c] + stcoa[c] '	5

% (13) malcoa is implied to be coupled with very long chain FA elongation
% instead of FASN. see how model reacts
% 'malcoa[c]'	'RCC0019'	'RCC0063'	'P/C'	0.674902391693208	1x5 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'3 h[c] + 2 nadph[c] + malcoa[c] + strdnccoa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + eicostetcoa[c] '	13
% 'malcoa[c]'	'RCC0019'	'RCC0064'	'P/C'	0.674902391693208	1x5 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'3 h[c] + 2 nadph[c] + malcoa[c] + lnlncgcoa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + dlnlcgcoa[c] '	13
% 'malcoa[c]'	'RCC0019'	'RCC0070'	'P/C'	0.761844021227047	1x13 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'3 h[c] + 2 nadph[c] + malcoa[c] + fa16p1n7coa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + vacccoa[c] '	13
% 'malcoa[c]'	'RCC0019'	'RCC0072'	'P/C'	0.761844021227047	1x5 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'3 h[c] + 2 nadph[c] + malcoa[c] + fa13p0iso[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + fa15p0iso[c] '	13
% 'malcoa[c]'	'RCC0019'	'RCC0073'	'P/C'	0.761844021227047	1x7 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'3 h[c] + 2 nadph[c] + malcoa[c] + fa15p0iso[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + fa17p0iso[c] '	13
% 'malcoa[c]'	'RCC0019'	'RCC0124'	'P/C'	0.674902391693208	1x4 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'9 h[c] + 6 nadph[c] + 3 malcoa[c] + stcoa[c]  -> 3 h2o[c] + 3 co2[c] + 6 nadp[c] + 3 coa[c] + lgnccoa[c] '	13
% 'malcoa[c]'	'RCC0019'	'RCC0210'	'P/C'	0.674902391693208	1x4 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'6 h[c] + 4 nadph[c] + 2 malcoa[c] + lgnccoa[c]  -> 2 h2o[c] + 2 co2[c] + 4 nadp[c] + 2 coa[c] + fa28p0coa[c] '	13

% (14) nucleitide interconversion is constrained with a few gene
% similarities. and also related to DNA syn. See how it influences the
% model 
% 'ins[c]'	'RC01863'	'RC01560'	'P/C'	0.508164616176076	1x1 table	'pi[c] + ins[c]  <=> r1p[c] + hxan[c] '	'h[c] + h2o[c] + adn[c]  -> nh4[c] + ins[c] '	4
% 'din[c]'	'RC02748'	'RC02556'	'P/C'	0.508164616176076	1x1 table	'pi[c] + din[c]  <=> hxan[c] + 2dr1p[c] '	'h[c] + h2o[c] + dad-2[c]  -> nh4[c] + din[c] '	3
% 'dad-2[c]'	'RC02557'	'RC02556'	'P/C'	0.508164616176076	1x1 table	'pi[c] + dad-2[c]  <=> ade[c] + 2dr1p[c] '	'h[c] + h2o[c] + dad-2[c]  -> nh4[c] + din[c] '	4
% 'datp[c]'	'RC01137'	'BIO0004'	'P/C'	0.309397287122737	1x3 table	'atp[c] + dadp[c]  <=> adp[c] + datp[c] '	'0.368 atp[c] + 0.368 h2o[c] + 0.325 datp[c] + 0.175 dgtp[c] + 0.175 dctp[c] + 0.325 dttp[c]  -> 0.368 adp[c] + 0.368 h[c] + 0.368 pi[c] + ppi[c] + DNA[c] '	3
% 'xmp[c]'	'RC01130'	'RC01230'	'P/C'	0.239644988901667	1x1 table	'h2o[c] + nad[c] + imp[c]  -> h[c] + nadh[c] + xmp[c] '	'atp[c] + nh4[c] + xmp[c]  -> 2 h[c] + amp[c] + ppi[c] + gmp[c] '	6
% 'xmp[c]'	'RC01130'	'RC01231'	'P/C'	0.239644988901667	1x1 table	'h2o[c] + nad[c] + imp[c]  -> h[c] + nadh[c] + xmp[c] '	'atp[c] + h2o[c] + gln-L[c] + xmp[c]  -> 2 h[c] + amp[c] + ppi[c] + glu-L[c] + gmp[c] '	6
% 'dgtp[c]'	'RC01857'	'BIO0004'	'P/C'	0.309397287122737	1x3 table	'atp[c] + dgdp[c]  -> adp[c] + dgtp[c] '	'0.368 atp[c] + 0.368 h2o[c] + 0.325 datp[c] + 0.175 dgtp[c] + 0.175 dctp[c] + 0.325 dttp[c]  -> 0.368 adp[c] + 0.368 h[c] + 0.368 pi[c] + ppi[c] + DNA[c] '	3


% good predictions 
% TCA cycle 
% 'oaa[m]'	'RM00351'	'RM00342'	'P/C'	0.411474101112489	1x1 table	'accoa[m] + h2o[m] + oaa[m]  -> coa[m] + h[m] + cit[m] '	'nad[m] + mal-L[m]  <=> nadh[m] + h[m] + oaa[m] '	4
% 'cit[m]'	'RM00351'	'RMC0001'	'P/C'	0.404391889848461	1x1 table	'accoa[m] + h2o[m] + oaa[m]  -> coa[m] + h[m] + cit[m] '	'cit[m]  <=> icit[m] '	5
% 'fum[m]'	'RM02164'	'RM01082'	'P/C'	0.325678505431294	5x1 table	'succ[m] + q[m]  -> qh2[m] + fum[m] '	'mal-L[m]  <=> h2o[m] + fum[m] '	8
% ETC
% 'q[m]'	'RM02164'	'RM02161'	'P/C'	0.805107311734685	5x8 table	'succ[m] + q[m]  -> qh2[m] + fum[m] '	'2 h[m] + qh2[m] + 2 ficytC[m]  -> 4 h[c] + q[m] + 2 focytC[m] '	12
% 'q[m]'	'RMC0006'	'RM02161'	'P/C'	0.725574680926681	29x8 table	'nadh[m] + 5 h[m] + q[m]  -> 4 h[c] + nad[m] + qh2[m] '	'2 h[m] + qh2[m] + 2 ficytC[m]  -> 4 h[c] + q[m] + 2 focytC[m] '	12
% 'q[m]'	'RM02161'	'RMC0007'	'P/C'	0.725574680926681	8x29 table	'2 h[m] + qh2[m] + 2 ficytC[m]  -> 4 h[c] + q[m] + 2 focytC[m] '	'h[c] + nadh[m] + q[m]  -> nad[m] + qh2[m] '	12
% 'q[m]'	'RM02161'	'RMC0008'	'P/C'	0.725574680926681	8x29 table	'2 h[m] + qh2[m] + 2 ficytC[m]  -> 4 h[c] + q[m] + 2 focytC[m] '	'nadh[m] + h[m] + q[m]  -> nad[m] + qh2[m] '	12
% PPP
% '6pgl[c]'	'RC02736'	'RC02035'	'P/C'	0.310545598955383	1x2 table	'g6p-B[c] + nadp[c]  -> h[c] + nadph[c] + 6pgl[c] '	'h2o[c] + 6pgl[c]  -> h[c] + 6pgc[c] '	2
% '6pgc[c]'	'RC02035'	'RC01528'	'P/C'	0.236161974244792	2x1 table	'h2o[c] + 6pgl[c]  -> h[c] + 6pgc[c] '	'nadp[c] + 6pgc[c]  -> co2[c] + nadph[c] + ru5p-D[c] '	3
% and others
% chitin syn (reveals hxk-1 coupling)
% 'gam6p[c]'	'RC00768'	'RC02058'	'P/C'	0.612605434668005	2x2 table	'f6p[c] + gln-L[c]  -> glu-L[c] + gam6p[c] '	'accoa[c] + gam6p[c]  -> h[c] + coa[c] + acgam6p[c] '	6
% 'gam6p[c]'	'RC01961'	'RC02058'	'P/C'	0.613175977760319	1x2 table	'atp[c] + gam[c]  -> adp[c] + h[c] + gam6p[c] '	'accoa[c] + gam6p[c]  -> h[c] + coa[c] + acgam6p[c] '	6
% folate cycle 
% '5mthf[c]'	'RC00946'	'RC01224'	'P/C'	0.609931610436158	3x2 table	'hcys-L[c] + 5mthf[c]  -> h[c] + met-L[c] + thf[c] '	'2 h[c] + nadph[c] + mlthf[c]  -> nadp[c] + 5mthf[c] '	4
% '5mthf[c]'	'RC00946'	'RC07168'	'P/C'	0.609931610436158	3x2 table	'hcys-L[c] + 5mthf[c]  -> h[c] + met-L[c] + thf[c] '	'2 h[c] + nadh[c] + mlthf[c]  -> nad[c] + 5mthf[c] '	4
% met/sams to pc 
% 'ahcys[c]'	'RC00192'	'RC02037'	'P/C'	0.465621903520335	1x1 table	'h2o[c] + ahcys[c]  -> hcys-L[c] + adn[c] '	'amet[c] + ethamp[c]  -> h[c] + ahcys[c] + methamp[c] '	15
% 'ahcys[c]'	'RC00192'	'RCC0021'	'P/C'	0.507813358622038	1x1 table	'h2o[c] + ahcys[c]  -> hcys-L[c] + adn[c] '	'2 amet[c] + methamp[c]  -> 2 h[c] + 2 ahcys[c] + cholp[c] '	15
% 'amet[c]'	'RC00177'	'RC02037'	'P/C'	0.569789300587654	4x1 table	'atp[c] + h2o[c] + met-L[c]  -> pi[c] + ppi[c] + amet[c] '	'amet[c] + ethamp[c]  -> h[c] + ahcys[c] + methamp[c] '	16
% 'amet[c]'	'RC00177'	'RCC0021'	'P/C'	0.460663269582792	4x1 table	'atp[c] + h2o[c] + met-L[c]  -> pi[c] + ppi[c] + amet[c] '	'2 amet[c] + methamp[c]  -> 2 h[c] + 2 ahcys[c] + cholp[c] '	16
% BCAA deg 
% 'ivcoa[m]'	'RMC0011'	'RM04096'	'P/C'	0.366269350862776	4x2 table	'coa[m] + nad[m] + 4mop[m]  -> co2[m] + nadh[m] + ivcoa[m] '	'ivcoa[m] + etfox[m]  -> 3mb2coa[m] + etfrd[m] '	4
% 'ivcoa[m]'	'RM04096'	'RMC0105'	'P/C'	0.204113332978664	2x4 table	'ivcoa[m] + etfox[m]  -> 3mb2coa[m] + etfrd[m] '	'2 coa[m] + 2 nad[m] + 2 h2o[m] + 2 fad[m] + fa9p0isocoa[m]  -> 2 accoa[m] + 2 nadh[m] + 2 h[m] + 2 fadh2[m] + ivcoa[m] '	4
% '3mb2coa[m]'	'RM04096'	'RM04138'	'P/C'	0.213892985033434	2x2 table	'ivcoa[m] + etfox[m]  -> 3mb2coa[m] + etfrd[m] '	'atp[m] + hco3[m] + 3mb2coa[m]  -> h[m] + adp[m] + pi[m] + 3mgcoa[m] '	2
% tyr and phe deg 
% '34hpp[c]'	'RC00734'	'RC02521'	'P/C'	0.437593126043506	2x1 table	'akg[c] + tyr-L[c]  <=> glu-L[c] + 34hpp[c] '	'o2[c] + 34hpp[c]  -> co2[c] + hgentis[c] '	5
% 'phpyr[c]'	'RC00694'	'RC01372'	'P/C'	0.437593126043506	2x1 table	'akg[c] + phe-L[c]  <=> glu-L[c] + phpyr[c] '	'o2[c] + phpyr[c]  -> co2[c] + 2hyoxplac[c] '	5
% mevolonate, udgp, gdpmannn to n-glycan 
% 'ipdp[c]'	'RC01121'	'RC05556'	'P/C'	0.247223549792247	1x1 table	'atp[c] + 5dpmev[c]  -> adp[c] + pi[c] + co2[c] + ipdp[c] '	'19 ipdp[c] + frdp[c]  -> 19 ppi[c] + dedoldp[c] '	7
% 'gdpmann[c]'	'RC05972'	'RC00885'	'P/C'	0.888561090465890	1x1 table	'chito2pdol[c] + gdpmann[c]  -> h[c] + gdp[c] + mpdol[c] '	'h[c] + gtp[c] + man1p[c]  -> ppi[c] + gdpmann[c] '	7
% 'gdpmann[c]'	'RC05973'	'RC00885'	'P/C'	0.779410925866119	1x1 table	'gdpmann[c] + mpdol[c]  -> h[c] + gdp[c] + m1mpdol[c] '	'h[c] + gtp[c] + man1p[c]  -> ppi[c] + gdpmann[c] '	7
% 'gdpmann[c]'	'RC06238'	'RC00885'	'P/C'	0.779410925866119	1x1 table	'gdpmann[c] + m1mpdol[c]  -> h[c] + gdp[c] + m2mpdol[c] '	'h[c] + gtp[c] + man1p[c]  -> ppi[c] + gdpmann[c] '	7
% 'gdpmann[c]'	'RCC0037'	'RC00885'	'P/C'	0.420349242047196	1x1 table	'2 gdpmann[c] + m2mpdol[c]  -> 2 h[c] + 2 gdp[c] + m4mpdol[c] '	'h[c] + gtp[c] + man1p[c]  -> ppi[c] + gdpmann[c] '	7
% 'gdpmann[c]'	'RC01009'	'RC00885'	'P/C'	0.586272423831442	1x1 table	'dolp[c] + gdpmann[c]  -> gdp[c] + dolmanp[c] '	'h[c] + gtp[c] + man1p[c]  -> ppi[c] + gdpmann[c] '	7
% 'udpg[c]'	'RC01005'	'RC00289'	'P/C'	0.811145061649767	1x1 table	'dolp[c] + udpg[c]  -> udp[c] + dolglcp[c] '	'g1p[c] + h[c] + utp[c]  -> ppi[c] + udpg[c] '	10
% 'man6p[c]'	'RC01326'	'RC01818'	'P/C'	0.234741099054269	1x1 table	'atp[c] + man[c]  -> adp[c] + h[c] + man6p[c] '	'man1p[c]  <=> man6p[c] '	3
% 'man1p[c]'	'RC01818'	'RC00885'	'P/C'	0.634268369450233	1x1 table	'man1p[c]  <=> man6p[c] '	'h[c] + gtp[c] + man1p[c]  -> ppi[c] + gdpmann[c] '	2
% FA oxidation and the oxidation reagent 
% 'ficytb[c]'	'RC00100'	'RC02222'	'P/C'	0.580774015810685	2x2 table	'h[c] + nadh[c] + 2 ficytb[c]  -> nad[c] + 2 focytb[c] '	'o2[c] + 2 focytb[c] + stcoa[c]  -> 2 h2o[c] + 2 ficytb[c] + odecoa[c] '	12
% 'ficytb[c]'	'RC00100'	'RCC0061'	'P/C'	0.228898886594997	2x1 table	'h[c] + nadh[c] + 2 ficytb[c]  -> nad[c] + 2 focytb[c] '	'o2[c] + 2 focytb[c] + odecoa[c]  -> 2 h2o[c] + 2 ficytb[c] + lnlccoa[c] '	12
% 'ficytb[c]'	'RC00100'	'RCC0069'	'P/C'	0.297508166554840	2x1 table	'h[c] + nadh[c] + 2 ficytb[c]  -> nad[c] + 2 focytb[c] '	'o2[c] + pmtcoa[c] + 2 focytb[c]  -> 2 h2o[c] + 2 ficytb[c] + fa16p1n7coa[c] '	12
% 'focytb[c]'	'RC00100'	'RC02222'	'P/C'	0.580774015810685	2x2 table	'h[c] + nadh[c] + 2 ficytb[c]  -> nad[c] + 2 focytb[c] '	'o2[c] + 2 focytb[c] + stcoa[c]  -> 2 h2o[c] + 2 ficytb[c] + odecoa[c] '	12
% 'focytb[c]'	'RC00100'	'RCC0061'	'P/C'	0.228898886594997	2x1 table	'h[c] + nadh[c] + 2 ficytb[c]  -> nad[c] + 2 focytb[c] '	'o2[c] + 2 focytb[c] + odecoa[c]  -> 2 h2o[c] + 2 ficytb[c] + lnlccoa[c] '	12
% 'focytb[c]'	'RC00100'	'RCC0069'	'P/C'	0.297508166554840	2x1 table	'h[c] + nadh[c] + 2 ficytb[c]  -> nad[c] + 2 focytb[c] '	'o2[c] + pmtcoa[c] + 2 focytb[c]  -> 2 h2o[c] + 2 ficytb[c] + fa16p1n7coa[c] '	12
% GSH syn 
% 'glucys[c]'	'RC00894'	'RC00497'	'P/C'	0.261174647006175	2x1 table	'atp[c] + glu-L[c] + cys-L[c]  -> adp[c] + h[c] + pi[c] + glucys[c] '	'atp[c] + gly[c] + glucys[c]  -> adp[c] + h[c] + pi[c] + gthrd[c] '	4

% many cases where degree = 2 is as expected because they are linear
% pathway

