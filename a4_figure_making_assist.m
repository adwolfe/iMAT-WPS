%% NOTE
% the confidence level of PPP back flux and msa degradation was assigned
% manually to high because they are buffered by alternative reactions
% (nadph vs nadh and g6p-A vs g6p-B). we manually confirmed that the FVA is
% bounded when one of the alternative is blocked

%% load
load('myCSM_triple.mat')
load('./../input/makeWormModel/iCEL1314_withUptakes.mat');
flux_triple = myCSM_triple.OFD ./ abs(myCSM_triple.OFD(strcmp(model.rxns,'EXC0050')));
annTbl = readtable('output\fluxTable_annotated.csv','ReadRowNames',true);

%% check by metabolite
met = 'g1p[c]'; % 'accoa[c]' leu-L
tbl3 = listRxn(model,flux_triple,met); % myCSM_merged.OFD
tbl3 = cell2table(tbl3);
tbl3.length10 = log10(1000 * abs(tbl3.tbl35));
tbl3.length2 = 0.5 * log2(200 * abs(tbl3.tbl35)); % anchor 0.01 at 0.5 with an increment of 0.5 per 2-fold increase
% add confidence level 
tbl3.confidence_PFD = annTbl.triple_PFD_bounded(tbl3.tbl31)
% the confidence level of lumped reactions are manually determined, mainly
% based on if the input metabolite has to be converted to the output no
% matter what alternative path to choose in the middle. 



%% lumped fluxes 
%% PPP - F6P
rxns = {'RC01830','RC01827'};
tbl3 = listRxn(model,flux_triple,'f6p[c]');
tbl3 = cell2table(tbl3);
len = 0.5 * log2(200 * sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, rxns)))))
%% CYS-L FROM AA FOOD
tbl3 = listRxn(model,flux_triple,'cys-L[c]');
tbl3 = cell2table(tbl3);
flux = abs(tbl3.tbl35(strcmp('RC00893', tbl3.tbl31))) - abs(tbl3.tbl35(strcmp('RC01001', tbl3.tbl31)));
len = 0.5 * log2(200 * flux)
%% beta oxidation - accoa 
rxns = {'RMC0106','RM04747','RM01177','RMC0107'};
tbl3 = listRxn(model,flux_triple,'accoa[m]');
tbl3 = cell2table(tbl3);
len = 0.5 * log2(200 * sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, rxns)))))
%% akg[m] efflux
tbl3 = listRxn(model,flux_triple,'akg[m]');
tbl3 = cell2table(tbl3);
flux = abs(tbl3.tbl35(strcmp('RM00267', tbl3.tbl31))) + tbl3.tbl35(strcmp('RM00243', tbl3.tbl31)) - abs(tbl3.tbl35(strcmp('RMC0003', tbl3.tbl31)));
len = 0.5 * log2(200 * flux)
%% BCAA mito total flux 
rxns = {'RM01090', 'RM01214','RM02199'};
flux = sum(abs(flux_triple(ismember(model.rxns,rxns))));
len = 0.5 * log2(200 * flux)

%% MET-L FROM AA FOOD
tbl3 = listRxn(model,flux_triple,'met-L[c]');
tbl3 = cell2table(tbl3);
flux = abs(tbl3.tbl35(strcmp('RC00177', tbl3.tbl31))) - abs(tbl3.tbl35(strcmp('RC00946', tbl3.tbl31)));
len = 0.5 * log2(200 * flux)

%% protein for deg
flux = flux_triple(strcmp(model.rxns,'DGR0008')) - sum(abs(flux_triple(ismember(model.rxns,{'BIO0002'  'BIO0001'}))));
len = 0.5 * log2(200 * flux)

%% glycine FROM AA FOOD
tbl3 = listRxn(model,flux_triple,'gly[c]');
tbl3 = cell2table(tbl3);
flux = abs(tbl3.tbl35(strcmp('DGR0001', tbl3.tbl31))) - abs(tbl3.tbl35(strcmp('RC03654', tbl3.tbl31))) - abs(flux_triple(ismember(model.rxns,'RM03654')));
len = 0.5 * log2(200 * flux)
%% lysine for deg
flux = flux_triple(strcmp(model.rxns,'TCM0089')) - sum(abs(flux_triple(ismember(model.rxns,{'RM03658'}))));
len = 0.5 * log2(200 * flux)
%% tyrosine for deg
flux = abs(flux_triple(strcmp(model.rxns,'TCM1037'))) + (abs(flux_triple(ismember(model.rxns,{'RC00734'}))));
len = 0.5 * log2(200 * flux)

%% accoa from peroxi beta 
rxns = {'RCC0088', 'RCC0090', 'RC03858', 'RC04742', 'RC03778'};
flux = sum(abs(flux_triple(ismember(model.rxns,rxns))));
len = 0.5 * log2(200 * flux)
%% malcoa go to ivcoa --> BCFA
rxns = {'RCC0071', 'RCC0072', 'RCC0073'};
flux = sum(abs(flux_triple(ismember(model.rxns,rxns))));
len = 0.5 * log2(200 * flux)

%% phe for deg
flux = abs(flux_triple(strcmp(model.rxns,'RC01795'))) + (abs(flux_triple(strcmp(model.rxns,'RC07211'))));
len = 0.5 * log2(200 * flux)

%% glycan flux
% uacgam
rxns = {'RC05969','RC05970','RC05916','RC05983','RC05985'};
flux = sum(abs(flux_triple(ismember(model.rxns,rxns))));
len = 0.5 * log2(200 * flux)

% gdpmann all goes to N-glycan

% udgp
flux1 = sum(abs(flux_triple(ismember(model.rxns,'RC01005'))));
flux2 = sum(abs(flux_triple(ismember(model.rxns,'RC05989'))));
len = 0.5 * log2(200 * (flux1 + flux2 * 2))
%% net comsumption of ipdp to dolp and to q9
dolp = {'RC05556'};
q9 = {'RCC0185'};
flux_dolp = sum(abs(flux_triple(ismember(model.rxns,dolp)))) * (19 + 1 * 3); % each dolp uses 19 ipdp and 1 frdp that uses 3 ipdp
flux_q9 = sum(abs(flux_triple(ismember(model.rxns,q9)))) * (6 + 1 * 3);

dolp_len = 0.5 * log2(200 * flux_dolp)
dolp_q9 = 0.5 * log2(200 * flux_q9)

%% collagan akg
rxns = {'RC03219','RC03376'};
tbl3 = listRxn(model,flux_triple,'akg[c]');
tbl3 = cell2table(tbl3);
len = 0.5 * log2(200 * sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, rxns)))))

%% total protein flux
rxns = {'BIO0001','BIO0002','RCC0171'};
flux = sum(abs(flux_triple(ismember(model.rxns,rxns))));
len = 0.5 * log2(200 * flux)


%% RNA to ura
rxns = {'RC01055','RC01876'};
tbl3 = listRxn(model,flux_triple,'ura[c]');
tbl3 = cell2table(tbl3);
len = 0.5 * log2(200 * sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, rxns)))))


%% G1P from glycogen
flux = abs(flux_triple(strcmp(model.rxns,'DGR0028'))) * (0.5 + 0.5*1*0.8*1); % one glycogen makes 0.5 g1p directly and 0.4 glp in the trimming f dxtrn1
len = 0.5 * log2(200 * flux)

%% length legend 
fluxes = 0.01 * 4 .^ (0:6);
lens = 0.5 * log2(200 .* fluxes);



%% protein analysis
% since the protein molelcules has similar MW, we can directly compare
% their flux
MolMass(model.metFormulas{strcmp(model.mets,'protein_BAC[c]')})
MolMass(model.metFormulas{strcmp(model.mets,'prot_mito[m]')})
MolMass(model.metFormulas{strcmp(model.mets,'prot_other[c]')})
log10(1000*(5.03712600000000 - 3.06654855000000 - 0.877364150000000))
% 'DGR0001'	'h2o[c] + protein_BAC[c]  -> 0.049203 gln-L[c] + 0.049203 glu-L[c] + 0.04507 asp-L[c] + 0.096044 ala-L[c] + 0.114544 gly[c] + 0.04507 asn-L[c] + 0.040346 ser-L[c] + 0.047432 thr-L[c] + 0.017123 cys-L[c] + 0.028735 met-L[c] + 0.04133 pro-L[c] + 0.017713 his-L[c] + 0.025782 tyr-L[c] + 0.034639 phe-L[c] + 0.010628 trp-L[c] + 0.055304 arg-L[c] + 0.084235 leu-L[c] + 0.05432 ile-L[c] + 0.064161 lys-L[c] + 0.079118 val-L[c] '	-5.03712600000000	3.70218281455045
% 'BIO0002'	'0.324 atp[c] + 2.324 h2o[c] + 2 gtp[c] + 0.027882 tyrtrna[c] + 0.012181 mettrna[c] + 0.051284 sertrna[c] + 0.084822 asptrna[c] + 0.074691 glytrna[c] + 0.045145 protrna[c] + 0.000208 cystrna[c] + 0.085568 glutrna[c] + 0.026184 glntrna[c] + 0.025739 asntrna[c] + 0.042387 argtrna[c] + 0.007259 trptrna[c] + 0.036544 phetrna[c] + 0.017927 histrna[c] + 0.052305 thrtrna[c] + 0.079867 leutrna[c] + 0.052143 iletrna[c] + 0.079083 lystrna[c] + 0.12769 alatrna[c] + 0.071091 valtrna[c]  -> 0.324 adp[c] + 3.324 h[c] + 2.324 pi[c] + 2 gdp[c] + 0.027882 trnatyr[c] + 0.012181 trnamet[c] + 0.051284 trnaser[c] + 0.084822 trnaasp[c] + 0.074691 trnagly[c] + 0.045145 trnapro[c] + 0.000208 trnacys[c] + 0.085568 trnaglu[c] + 0.026184 trnagln[c] + 0.025739 trnaasn[c] + 0.042387 trnaarg[c] + 0.007259 trnatrp[c] + 0.036544 trnaphe[c] + 0.017927 trnahis[c] + 0.052305 trnathr[c] + 0.079867 trnaleu[c] + 0.052143 trnaile[c] + 0.079083 trnalys[c] + 0.12769 trnaala[c] + 0.071091 trnaval[c] + prot_other[c] '	3.06654855000000	3.48664984488530
% 'BIO0001'	0.877364150000000	'2.324 h2o[m] + 0.324 atp[m] + 2 gtp[m] + 0.026832 tyrtrna[m] + 0.020458 mettrna[m] + 0.045815 sertrna[m] + 0.064946 asptrna[m] + 0.094618 glytrna[m] + 0.043238 protrna[m] + 0.008665 cystrna[m] + 0.067386 glutrna[m] + 0.037694 glntrna[m] + 0.035405 asntrna[m] + 0.048846 argtrna[m] + 0.008944 trptrna[m] + 0.035592 phetrna[m] + 0.01782 histrna[m] + 0.049868 thrtrna[m] + 0.082051 leutrna[m] + 0.053231 iletrna[m] + 0.071622 lystrna[m] + 0.111867 alatrna[m] + 0.075105 valtrna[m]  -> 3.324 h[m] + 0.324 adp[m] + 2.324 pi[m] + 2 gdp[m] + 0.026832 trnatyr[m] + 0.020458 trnamet[m] + 0.045815 trnaser[m] + 0.064946 trnaasp[m] + 0.094618 trnagly[m] + 0.043238 trnapro[m] + 0.008665 trnacys[m] + 0.067386 trnaglu[m] + 0.037694 trnagln[m] + 0.035405 trnaasn[m] + 0.048846 trnaarg[m] + 0.008944 trnatrp[m] + 0.035592 trnaphe[m] + 0.01782 trnahis[m] + 0.049868 trnathr[m] + 0.082051 trnaleu[m] + 0.053231 trnaile[m] + 0.071622 trnalys[m] + 0.111867 trnaala[m] + 0.075105 trnaval[m] + prot_mito[m] '	0.877364150000000	2.94317988471303


%% case studies of metabolite balances

%% nadph balance and PPP
met = ['f6p[c]']; % 'accoa[c]' leu-L
tbl3 = listRxn(model,flux_triple,met); % myCSM_merged.OFD
tbl3 = cell2table(tbl3);
tbl3.length = 3 * abs(tbl3.tbl35); % anchor 0.01 at 0.5 with an increment of 0.5 per 2-fold increase
%% nadph balance 1
nonPPPrxns = {'RC08759','RC02236','RC08383','RC00941'};
PPPrxns = {'RC02736','RC01528'};
tbl3 = listRxn(model,flux_triple,'nadph[c]');
tbl3 = cell2table(tbl3);
nonPPP = 3 * sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, nonPPPrxns))))
PPP = 3 * sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, PPPrxns))))

%% nadph balance 2
mevrnxs = {'RC02082'};
FAsynrxns = {'RCC0020','RCC0064','RC07759','RC07761','RCC0071','RCC0073','RCC0070','RCC0124','RCC0072'};
tbl3 = listRxn(model,flux_triple,'nadph[c]');
tbl3 = cell2table(tbl3);
FAsyn = 3 *sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, FAsynrxns))))
mev = 3 * sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, mevrnxs))))
others = 3 * (sum(abs(tbl3.tbl35))/2 - sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, [FAsynrxns,mevrnxs])))))


%% nadh balance in mitochondria
met = 'nadh[m]'; % 'accoa[c]' leu-L
tbl3 = listRxn(model,flux_triple,met); % myCSM_merged.OFD
tbl3 = cell2table(tbl3);
tbl3.length = abs(tbl3.tbl35); % anchor 0.01 at 0.5 with an increment of 0.5 per 2-fold increase
%% etfrd balance
AAdeg = {'RM04432','RM02661','RM04096','RM03172','RM02488'};
fadeg = {'RM04751','RM03777','RM01175','RM04754'};
tbl3 = listRxn(model,flux_triple,'etfrd[m]');
tbl3 = cell2table(tbl3);
AAdeg_len = sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, AAdeg))))
fadeg_len = sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, fadeg))))
%% fadh2 from FA deg
rxns = {'RMC0106','RMC0107'};
tbl3 = listRxn(model,flux_triple,'fadh2[m]');
tbl3 = cell2table(tbl3);
f1 = sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, rxns))))
fadeg_len_sum = fadeg_len + f1



%% atp flux
met = 'atp[c]'; % 'accoa[c]' leu-L
tbl3 = listRxn(model,flux_triple,met); % myCSM_merged.OFD
tbl3 = cell2table(tbl3);
tbl3.length = 0.5 * abs(tbl3.tbl35); % anchor 0.01 at 0.5 with an increment of 0.5 per 2-fold increase

%% total ATP for tRNA 
% each tRNA syn comsumes 2 ATP since 1 atp goes to 1 amp and then 1 amp
% needs another atp to be back to adp
rxns = {'RC03656','RC03663','RC03662','RC03661','RC03665','RC03658','RC03657','RC05577','RC05578','RC03654','RC03038','RC03648','RC02918','RC03652','RC03660','RC03646','RC03655','RC03659','RC03664','RC03650'};
flux = 2 * sum(abs(flux_triple(ismember(model.rxns,rxns))));
len1 = 0.5 * flux
%% gtp and atp to protein syn 
rxns = {'BIO0002','RCC0170'}; % mito protein ignored since we are concerning cytosol now
met = 'gtp[c]'; % 'accoa[c]' leu-L
tbl3 = listRxn(model,flux_triple,met); % myCSM_merged.OFD
tbl3 = cell2table(tbl3);
len2 = 0.5 * sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, rxns))))
%%
met = 'atp[c]'; % 'accoa[c]' leu-L
tbl3 = listRxn(model,flux_triple,met); % myCSM_merged.OFD
tbl3 = cell2table(tbl3);
len3 = 0.5 * sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, rxns))))
%% total ATP for lipid synthesis 
rxns = {'RC00352','RCC0019','RC01280'};
tbl3 = listRxn(model,flux_triple,'atp[c]');
tbl3 = cell2table(tbl3);
len4 = 0.5 * sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, rxns))))
%% others
included = {'BIO0010','RCC0005','DGR0007'};
tbl3 = listRxn(model,flux_triple,'atp[c]');
tbl3 = cell2table(tbl3);
len5 = 0.5 * sum(abs(tbl3.tbl35(ismember(tbl3.tbl31, included))))
others = sum(abs(tbl3.tbl35)) / 2 - (len1 + len2 + len3 + len4 + len5)
others / (sum(abs(tbl3.tbl35)) / 2)



