function expressionCol = selectGeneFromGPR_DE_xl(model, gene_names, gene_exp, parsedGPR,responsiveGenes, nonresponsiveGenes)
% Map gene expression to reaction expression using the GPR rules. The
% mapping principles in IMAT++ algorithm is used. This function is
% developed by modifying the original COBRA function "selectGeneFromGPR".
%
% USAGE:
%   expressionCol = selectGeneFromGPR_xl(model, gene_names, gene_exp, parsedGPR)
%
% INPUTS:
%   model:          COBRA model struct
%   gene_names:     gene identifiers corresponding to gene_exp. Names must
%                   be in the same format as model.genes (column vector)
%                   (as returned by "findUsedGeneLevels.m")
%   gene_exp:       quodra-level expression
%   parsedGPR:      GPR matrix as returned by "GPRparser.m"
%
% OUTPUTS:
%   expressionCol:  reaction expression, corresponding to model.rxns.
%                   No gene-expression data and orphan reactions will
%                   be given a value of -1.
%
% Original AUTHOR: Anne Richelle, May 2017
% AUTHOR: Xuhang Li, 08192019 
%       - optimized for speed: remove the GPR parsing part, and requires a pre-parsed field in the input model
%       - added new lines to categorize reactions associated with "high
%         and moderate genes" as high, instead of moderate in the IMAT
%         algorithm
%       - now the parser is only suitable for IMAT++ integration, where
%         genes are categorized and represented as "0"(rarely expressed),
%         "1"(lowly expressed), "2"(moderately) and "3"(highly). IT IS NOT
%         SUITABLE FOR ANY CONTINUOUS EXPRESSION LEVEL INPUT!!

% we first consider all responsive genes to be high gene 
gene_exp(ismember(gene_names,responsiveGenes)) = 3;

expressionCol = -1*ones(length(model.rxns),1); %-1 means unknown/no data
    for i = 1:length(model.rxns)
        mygenes = model.genes(logical(model.rxnGeneMat(i,:)));
        nonresp = nonresponsiveGenes(ismember(nonresponsiveGenes, mygenes));
        if ~isempty(nonresp) % test if it has an effect on classification
            % we put a reaction to low reaction as long as any
            % nonresponsive gene is influencing this reaction by GPR (but
            % if it is originally zero, it stays zero)
            alltests = [];
            for k = 1:length(nonresp)
                gene_exp_test = gene_exp;
                gene_exp_test(strcmp(gene_names, nonresp{k})) = min(1, gene_exp(strcmp(gene_names, nonresp{k}))); % use the lower value between 1 (responsiveness evidence) or orginal level (absolute expression)
                curExprArr=parsedGPR{i};
                curExpr= nan(length(curExprArr),1);
                for j=1:length(curExprArr)
                    if length(curExprArr{j})>=1
                        geneID = ismember(gene_names,curExprArr{j});
                        if any(geneID) %if the gene is measured
                             curExpr(j)= min(gene_exp_test(geneID)); %If there is data for any gene in 'AND' rule, take the minimum value
                        end
                    end
                end
                if any(~isnan(curExpr))
                    alltests(k)=max(curExpr);%If there is data for any gene in the 'OR' rule, take the maximum value
                end
            end
            expressionCol(i)=min(alltests); % if there is one RNAi evidence supporting this reaction not carrying flux, it is not carrying flux
        else % no gene was determined for RNAi evidance, do normal parsing 
            curExprArr=parsedGPR{i};
            curExpr= nan(length(curExprArr),1);
            for j=1:length(curExprArr)
                if length(curExprArr{j})>=1
                    geneID = ismember(gene_names,curExprArr{j});
                    if any(geneID) %if the gene is measured
                         curExpr(j)= min(gene_exp(geneID)); %If there is data for any gene in 'AND' rule, take the minimum value
                    end
                end
            end
            if any(~isnan(curExpr))
                expressionCol(i)=max(curExpr);%If there is data for any gene in the 'OR' rule, take the maximum value
            end
        end
    end
    %now, all the reactions are labeled. 
    %next, we retrieve all the reactions associated with high expression
    %genes(level = 3), label them as high reactions, and only exluded low and zero reactions. Therefore, high
    %reaction will overide moderate level genes.
    HGenes = gene_names(gene_exp == 3);
    HrxnsInd = any(model.rxnGeneMat(:,ismember(model.genes,HGenes)),2);
    NonLowZero = expressionCol ~= 0 & expressionCol ~= 1; %any one that is not a zero or low reaction
    HrxnsCandidate = HrxnsInd & NonLowZero;
    expressionCol(HrxnsCandidate) = 3; %all HGene associated reactions (no conflict with low and zero) are re-labeled as level = 3
    
    % ignore the below..
    % finally, we use the DE evidence for assigning responsive genes 
    % we consider all reactions associated with responsive genes as high,
    % so in the gene-centric fitting, at least one reaction will carry flux
    % ==> this is too aggressive, opens too many reaction. 
%     HrxnsInd = any(model.rxnGeneMat(:,ismember(model.genes,responsiveGenes)),2);
%     expressionCol(HrxnsInd) = 3; 

    % instead, we report the hard conflicts at the GPR level 
    % check if any responsive gene is conflicting with zero genes
    % note: there can be remaining conflicts with the nonresponsive genes,
    % but in these case, we think the responsive gene should overide the
    % non-responsive ones since the latter can be technical artifect or
    % non-tsxp buffering. we dont thraw a warning on this because we will
    % force fitting responsive genes later on (by MILP)
    for i = 1:length(responsiveGenes)
        myrxnInds = any(model.rxnGeneMat(:,ismember(model.genes,responsiveGenes(i))),2);
        if max(expressionCol(myrxnInds)) == 0
            fprintf('responsive gene %s is only associated with off (zero) reaction!\n',responsiveGenes{i});
        end
    end
    
end



