%% Overview
% analyze the responsiveness of genes in the iCEL model and convert it into
% iMAT++ categories. The absolute expression of genes in WT animals were
% also analyzed and used to make gene categories. This script produces the
% input files needed for running iMAT-WPS.

mkdir(['input/WPS/']);
addpath scripts/PlotPub/lib/

%% I: the absolute gene expression integration
%% load the absolute expression data
% we load the TPM calculated by combining the entire WPS dataset
TPM = readtable('input/WPS/WT_TPM.csv','ReadRowNames',1);
TPM.Properties.RowNames = TPM.WBID;
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');
% (specific to worm model)
% worm expression data are often labeled with WormBase ID (WBID). So, we
% provided a simple lookup table for iCEL1314.
lookupTbl = readtable('./input/model/IDtbl.csv');
metGenesInd = ismember(TPM.Properties.RowNames,lookupTbl.WormBase_Gene_ID);
wormbaseTbl = readtable('./input/model/wormbasetable.tsv','FileType','text');
proteinCodingGenes = wormbaseTbl.WormBaseID(strcmp(wormbaseTbl.Type,'protein coding'));
proteinCodingInd = ismember(TPM.Properties.RowNames,proteinCodingGenes);

%% gene categorization based on absolute expression
% to better capture the entire distribution of gene expression, we used the
% whole-dataset-averege TPM in the curve fitting (this is similar to what
% we did in the C. elegans tissue flux paper, Yilmaz, Li, et al., 2020,
% MSB)
fitData = TPM.whole_dataset_TPM(proteinCodingInd);% we perform guassian fitting only on metabolic genes
iCELData = TPM.whole_dataset_TPM(metGenesInd);% we perform guassian fitting only on metabolic genes
%fitData = fitData(:);
histogram(log2(fitData))
% fit bimodel guassian
rng(1126)
x = log2(fitData);% this automatically ignored all the zeros
x(isinf(x)) = [];
fit = fitgmdist(x,2,'Options',statset('Display','final','MaxIter',3000,'TolFun',1e-9),'Replicates',10);
% display the iCEL metabolic genes 
iCel = log2(iCELData);% this automatically ignored all the zeros
iCel(isinf(iCel)) = [];
% note: the sigma in output is sigma^2
%% visualization of the guassian distribution fitted
% user can inspect how well the fitting is and whether it seperates into
% two subpopulation as expected.
figure;
bins = -15:.5:15;
h = bar(bins,histc(x,bins)/(length(x)*.5),'histc');
h.FaceColor = [0 0.4470 0.7410];%[.9 .9 .9];
h.FaceAlpha = 0.25;
xgrid = linspace(1.1*min(x),1.1*max(x),200);
pdfgrid = pdf(fit,xgrid');
hold on
h2 = bar(bins,histc(iCel,bins)/(length(iCel)*.5),'histc');
h2.FaceColor = [0.8500 0.3250 0.0980];
h2.FaceAlpha = 0.25;
legend({'Protein coding genes','iCEL genes'})
plot(xgrid,pdfgrid,'-')
hold off
xlabel('TPM')
ylabel('Probability Density')

% label the desired cutoff - which is sigma1 and which is sigma2 may change
% if the data is changed
xline(fit.mu(1) + 1*sqrt(fit.Sigma(1)),'--k');
fit.mu(1) + 1*sqrt(fit.Sigma(1))
xline(fit.mu(2),'--k');
fit.mu(2)
xline(fit.mu(1),'--k');
fit.mu(1)


% save figure
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [4.5, 3];
plt.YLim = [0 0.13];
plt.XLim = [-17 17];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.LegendLoc = 'NorthWest';
plt.export(['figures/absolute_expression_categorization.pdf']);

% zero2low =
% 
%    -2.5759
% 
% 
% low2dynamic =
% 
%     1.6661
% 
% 
% dynamic2high =
% 
%     5.4007

%% build the gene catagories
zero2low = fit.mu(1) %set thresholds
low2dynamic = fit.mu(1) + 1*sqrt(fit.Sigma(1)) %set thresholds
dynamic2high = fit.mu(2) %set thresholds
metgenes = model.genes;
TPM = TPM(ismember(TPM.Properties.RowNames,lookupTbl.WormBase_Gene_ID),:);%we only keep the genes in the model
 % specifically, we need to convert the gene ID for C. elegans model
[A,B] = ismember(TPM.Properties.RowNames,lookupTbl.WormBase_Gene_ID);
TPM.GeneID(A) = lookupTbl.ICELgene(B(A));
myTPM = log2(TPM.vector_average_TPM);
GeneID = TPM.GeneID;
ExpCateg.zero = GeneID(myTPM < zero2low);
ExpCateg.low = GeneID(myTPM >= zero2low & myTPM < low2dynamic);
ExpCateg.dynamic = GeneID(myTPM >= low2dynamic & myTPM < dynamic2high);
ExpCateg.high = GeneID(myTPM >= dynamic2high);
% the uncalled genes (i.e., NA and ND) are in dynamic (moderately expressed)
ExpCateg.dynamic = [ExpCateg.dynamic; metgenes(~ismember(metgenes,GeneID))];
% save for use in iMAT-WPS
save(['input/WPS/categ_expression_only.mat'],'ExpCateg');

%% VISUALIZATION: save the pie chart of gene category
text = strcat(fieldnames(ExpCateg),'=',num2str(structfun(@length,ExpCateg)));
% create pie chart
pie(structfun(@length,ExpCateg), false(length(fieldnames(ExpCateg)),1), text);

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
plt.export(['figures/absolute_expression_categorization_piechart.pdf']);

%% II: integrate the DE responsiveness
%% first do some data cleaning before use
% the genes carrying flux
DEtbl = readtable('input/WPS/DE_merged_clean_pcutoff_0.005_master_table_FDR2d0.1_FC1.5_FCtype_log2FoldChange_raw_ALL.csv');
% conditionInfo = readtable('./../../MetabolicLibrary/2_DE/output/RNAi_condition_metaInfo.csv');
% 12302023: update to final publish-version data in which four conditions
% was wrong targeted gene annotations (all nonresponsive) were corrected
conditionInfo = readtable('input/WPS/RNAi_condition_information.csv');
% to be consistent with the model annotation, we still use the mrpl-44 as
% mccc-2
conditionInfo.RNAi_WBID(strcmp(conditionInfo.RNAiID,'x.mrpl_44 met3_lib1')) = {'WBGene00008514'};
conditionInfo.RNAi_geneName(strcmp(conditionInfo.RNAiID,'x.mrpl_44 met3_lib1')) = {'mrpl-44'};


% since the WPS data can be noisy (i.e., depends on seeding stage and RNAi
% efficiency), we place a few filters to enhance the data strigency

% (1) genes in responsive category
% First, the responsive conditions should be responsive category
% However, if it is vectorlike, recomb vector, multiple etc, we skip the gene 
responsive_WBID = conditionInfo.RNAi_WBID(strcmp(conditionInfo.isResponsive,'TRUE'));
[A,B] = ismember(responsive_WBID,lookupTbl.WormBase_Gene_ID);
responsive_genes = unique(lookupTbl.ICELgene(B(A)));

% Second, all phenotypic genes should be responsive category
pheno = readtable('input/WPS/not sequenced RNAis because of strong phenotype.csv','ReadVariableNames',1);
[A,B] = ismember(pheno.WBID,lookupTbl.WormBase_Gene_ID);
pheno_genes = unique(lookupTbl.ICELgene(B(A)));

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
ExpCateg
length(intersect(ToHighGenes, ExpCateg.high))
length(intersect(ToHighGenes, ExpCateg.dynamic))
length(intersect(ToHighGenes, ExpCateg.low))
length(intersect(ToHighGenes, ExpCateg.zero))

length(intersect(ToZeroGenes, ExpCateg.high))
length(intersect(ToZeroGenes, ExpCateg.dynamic))
length(intersect(ToZeroGenes, ExpCateg.low))
length(intersect(ToZeroGenes, ExpCateg.zero))

ExpCateg_expression = ExpCateg;
ExpCateg.responsive = ToHighGenes;
ExpCateg.nonresponsive = ToZeroGenes;
% save for use in iMAT-WPS
save(['input/WPS/categ_expression_and_WPS.mat'],'ExpCateg');

%% VISUALIZATION: plot the responsiveness pie chart
% pie chart of responsiveness 
figure;
tmp = rmfield(ExpCateg,{'zero','low','dynamic','high'});
text = strcat(fieldnames(tmp),'=',num2str(structfun(@length,tmp)));
% create pie chart
pie(structfun(@length,tmp), false(length(fieldnames(tmp)),1), text);
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
plt.export(['figures/responsiveness_piechart.pdf']);

% barplot of responsiveness integrated with the absolute expression
% category (this might be better where we use an arrow from expression
% fitting and from N_DE distribution merged to this figure)

data = structfun(@length,ExpCateg);
text = fieldnames(ExpCateg);
tbl = table('RowNames',text);
tbl.total = data;
tbl.responsive = structfun(@(x) length(intersect(x, ExpCateg.responsive)), ExpCateg);
tbl.nonresponsive = structfun(@(x) length(intersect(x, ExpCateg.nonresponsive)), ExpCateg);

% barplot is not informative. Let's make the stacked piechart in R. 
writetable(tbl,'input/WPS/gene_category_summary.csv','WriteRowNames',true);


%% Notes about conditions with <5 DEG but very high FC
% we didnt see a great value doing this - those conditions are often noisy and
% not very interpretable. so we don't specially process these conditions
% % if there is any very significant fc, we dont put to zero in case it is a
% % point-point response
% myTbl = DEtbl(ismember(DEtbl.RNAiID, conditionInfo.RNAiID(strcmp(conditionInfo.isResponsive, 'FALSE'))),:);
% % remove self
% myTbl.RNAi_geneName = myTbl.RNAi;
% myTbl.RNAi_geneName = regexprep(myTbl.RNAi_geneName,'^x.','');
% myTbl.RNAi_geneName = regexprep(myTbl.RNAi_geneName,'_','-');
% myTbl.RNAi_geneName = regexprep(myTbl.RNAi_geneName,'-L[0-9]$','');
% myTbl = myTbl(~cellfun(@isequal, myTbl.RNAi_geneName, myTbl.Gene_name),:);
% point2point = unique(myTbl.RNAiID(abs(myTbl.log2FoldChange_raw) > log2(10) & (strcmp(myTbl.isICEL,'TRUE') | strcmp(myTbl.isMetabolic,'TRUE'))));
% % ID conversion 
% [A,B] = ismember(point2point,conditionInfo.RNAiID);
% point2point(A,2) = conditionInfo.RNAi_WBID(B(A));
% % for now, if it is vectorlike, recomb vector, multiple etc, we skip the gene 
% point2point(~A,2) = {'NA'};
% [A,B] = ismember(point2point(:,2),lookupTbl.WormBase_Gene_ID);
% point2point_genes = lookupTbl.ICELgene(B(A));
% % point2point(~A,:)
