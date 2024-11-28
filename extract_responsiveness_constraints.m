function ExpCateg = extract_responsiveness_constraints(conditionInfo)
%% load required tables
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');
% (specific to worm model)
% worm expression data are often labeled with WormBase ID (WBID). So, we
% provided a simple lookup table for iCEL1314.
lookupTbl = readtable('./input/model/IDtbl.csv');

% load the category for absolute expression
load('input/WPS/categ_expression_only.mat');

% load DE table for refinement use
DEtbl = readtable('./input/WPS/DE_merged_clean_pcutoff_0.005_master_table_FDR2d0.1_FC1.5_FCtype_log2FoldChange_raw_ALL.csv');

%% II: integrate the DE responsiveness
%% first do some data cleaning before use

% seperate the regular samples from phenotypes in the subtable
idx = strcmp(conditionInfo.RNAiID,'PHENOTYPE');
sampled_pheno = conditionInfo.RNAi_WBID(idx); 
conditionInfo = conditionInfo(~idx, :);

% since the WPS data can be noisy (i.e., depends on seeding stage and RNAi
% efficiency), we place a few filters to enhance the data strigency

% (1) genes in responsive category
% First, the responsive conditions should be responsive category
% However, if it is vectorlike, recomb vector, multiple etc, we skip the gene 
responsive_WBID = conditionInfo.RNAi_WBID(strcmp(conditionInfo.isResponsive,'TRUE'));
[A,B] = ismember(responsive_WBID,lookupTbl.WormBase_Gene_ID);
responsive_genes = unique(lookupTbl.ICELgene(B(A)));

% Second, all phenotypic genes should be responsive category
pheno = readtable('./input/WPS/not sequenced RNAis because of strong phenotype.csv','ReadVariableNames',1);
[A,B] = ismember(pheno.WBID,lookupTbl.WormBase_Gene_ID);
pheno_genes = unique(lookupTbl.ICELgene(B(A)));
% keep only what is being in the sample
sampled_pheno = lookupTbl.ICELgene(ismember(lookupTbl.WormBase_Gene_ID, sampled_pheno));
pheno_genes = intersect(sampled_pheno, pheno_genes);


% (2) genes in non-responsive category
% First, the nonresponsive conditions should be non-responsive category
% However, as long as a gene is found responsive in one condition, we 
% should consider it  responsive and exclude it from non-resp. category
nonresponsive_WBID = conditionInfo.RNAi_WBID(strcmp(conditionInfo.isResponsive,'FALSE'));  
[A,B] = ismember(nonresponsive_WBID,lookupTbl.WormBase_Gene_ID);
nonresponsive_genes = lookupTbl.ICELgene(B(A));
nonresponsive_genes = setdiff(nonresponsive_genes, responsive_genes);

% Second, to remove potential problems introduced by RNAi efficiency, we 
% remove those genes who are highly expressed, non-responsive, but not sig. 
% down in their targeted gene expression in RNAi condition. 
DEtbl.RNAiID = strcat(DEtbl.RNAi, {' '}, DEtbl.batchID);
myTbl = DEtbl(ismember(DEtbl.RNAiID, conditionInfo.RNAiID(strcmp(conditionInfo.isResponsive, 'FALSE'))),:);
myTbl.RNAi_geneName = myTbl.RNAi;
myTbl.RNAi_geneName = regexprep(myTbl.RNAi_geneName,'^x.','');
myTbl.RNAi_geneName = regexprep(myTbl.RNAi_geneName,'_','-');
myTbl.RNAi_geneName = regexprep(myTbl.RNAi_geneName,'-L[0-9]$','');
myTbl = myTbl(cellfun(@isequal, myTbl.RNAi_geneName, myTbl.Gene_name),:);
self_downs = unique(myTbl.RNAiID(myTbl.log2FoldChange_raw < 0));
% ID conversion 
[A,B] = ismember(self_downs,conditionInfo.RNAiID);
self_downs(A,2) = conditionInfo.RNAi_WBID(B(A));
% if it is vectorlike, recomb vector, multiple etc, we skip the gene 
self_downs(~A,2) = {'NA'};
[A,B] = ismember(self_downs(:,2),lookupTbl.WormBase_Gene_ID);
nonresponsive_self_down = lookupTbl.ICELgene(B(A));
nonresponsive_self_not_down = setdiff(nonresponsive_genes, nonresponsive_self_down);
% we need to remove the highly expressed but not down genes 
nonresponsive_self_not_down_high_expressed = intersect(nonresponsive_self_not_down, ExpCateg.high);

% merge the sets
ToHighGenes = unique([responsive_genes; pheno_genes]); % this is the responsive gene set for integration use
ToZeroGenes = setdiff(nonresponsive_genes, nonresponsive_self_not_down_high_expressed);% this is the non-responsive gene set for integration use

% get a sense of impact - the following can be ignored
length(intersect(ToHighGenes, ExpCateg.high))
length(intersect(ToHighGenes, ExpCateg.dynamic))
length(intersect(ToHighGenes, ExpCateg.low))
length(intersect(ToHighGenes, ExpCateg.zero))

length(intersect(ToZeroGenes, ExpCateg.high))
length(intersect(ToZeroGenes, ExpCateg.dynamic))
length(intersect(ToZeroGenes, ExpCateg.low))
length(intersect(ToZeroGenes, ExpCateg.zero))

ExpCateg.responsive = ToHighGenes;
ExpCateg.nonresponsive = ToZeroGenes;

