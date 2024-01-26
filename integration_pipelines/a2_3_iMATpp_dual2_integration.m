function a2_3_iMATpp_dual2_integration(yield, relCap_minLow, relCap_metFit)
%% summary
% this integrates absolute expression + DE similarity information.
% allowing two relative cap parameters to control the interaction between
% the abs_exp/resp. fitting and the metabolite fitting.

%% run

% Load model - this is the same across all five integrations 
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');
load('./input/model/epsilon_generic_withUptakes.mat'); 
load('input/WPS/categ_expression_only.mat');
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
% we use the biomass yield rate to constrain the bacteria waste
% this is to force the nutrient to be efficiently used instead of wasted in
% bulk
% yield * V(EXC0050) + V(BIO0010) >= 0 (V(EXC0050) is a negative number)
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


% set up an empty responsiveness field to avoid error in integration
% function
ExpCateg.responsive = {};
ExpCateg.nonresponsive = {};

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
myCSM.wasteDW,...
myCSM.RHNames,...
myCSM.branchMets,...
myCSM.minimizedLowRxns]...
= IMATplusplus_wiring_triple_inetgration_final(model_coupled,epsilon_f,epsilon_r, ExpCateg, branchTbl, modelType,speedMode,...
relCap_minLow, relCap_metFit, 1, 1, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);
% NOTE
% since the absolute expression/responsiveness is fitted in sequential
% steps with the DE similarity (metabolite fit), we need to set the fitting
% cap properly to handle the interaction. If we use a very tight cap to
% keep perfect fit of the first fitted one (aka the abs_exp/resp), the next
% one can not be optimally fitted. Therefore, we need to adjust the cap
% according to specific needs for a simulation.

% here, since the primary aim is to fit the WPS data, we give DE similarity
% fitting a higher priority. We loosen the cap for the first fitted
% absolute expression (minLow) to 5% relative level to give sufficient 
% space for the second to be fitted. We set the reletive cap for the second
% to 1e-6, which essentially makes the abs minMetFit to be capped around
% 1e-5, thus, perfectly fitted. Applying rigid fitting of similarity data
% here is to focus on the effect of similarity integration in the benchmark
% analysis. 

myCSM_exp_simi = myCSM;
save('output/integration_output/myCSM_exp_simi.mat','myCSM_exp_simi');

%% codes for manual interactions
% %% check the fit of coupling and change of flux 
% load('myCSM_expression_only.mat')
% a = table(myCSM_exp.OFD, myCSM.OFD);
% a.delta = abs(a.Var1 - a.Var2);
% a.rxn = model.rxns;
% a.formula = printRxnFormula(model,model.rxns,0);
% a.bigJump = a.Var1 .* a.Var2 <= 0;
% 
% branchMets = unique(branchTbl.mets);
% MetFlux = [];
% weight = [];
% for i = 1:length(branchMets)
%     % first locate the reactions in the module 
%     myrxns = unique([branchTbl.rxn1(strcmp(branchTbl.mets, branchMets{i})); branchTbl.rxn2(strcmp(branchTbl.mets, branchMets{i}))]);
%     MetFlux(i,1) = sum(model.S(strcmp(model.mets, branchMets{i}), ismember(model.rxns, myrxns)) .* myCSM.OFD(ismember(model.rxns, myrxns))');
%     weight(i,1) = max(branchTbl.maxCosine(strcmp(branchTbl.mets, branchMets{i})));
% end
% metFit = table(branchMets, MetFlux, weight);
% metFit.contribution = abs(MetFlux) .* weight;
% abs(MetFlux)' * weight
% myCSM.minMetBalanceLoss + 1e-5
% 
% %% check cases of metabolite fitting
% % (1) for the apperent problematic couplings, the model tolerates them very
% % well (near perfectly)
% % it may need curation before applied. eg focytC[m] couplings OXPHOS with focytC
% % synthesis, which is obviously a bad coupling for flux wiring
% % the followings are the apparently unreasonable couplings 
% % (low cosine and unintuitive)
% % check if they can be buffered or they had a negative impact on the flux
% % integration
% % (VERY BAD) 'ala-L[c]'	'RC03038'	'RCC0142'	'P/C'	0.306444282724250	1x1 table	'atp[c] + ala-L[c] + trnaala[c]  -> amp[c] + ppi[c] + alatrna[c] '	'2 atp[c] + h2o[c] + 2 cys-L[c] + 2 trdrd[c] + cpmp[c]  -> 4 h[c] + 2 amp[c] + 2 ppi[c] + 2 ala-L[c] + 2 trdox[c] + molybd[c] '	15
% metFit(strcmp(metFit.branchMets,'ala-L[c]'),:)
% a(ismember(a.rxn, {'RC03038'	'RCC0142'}),:)
% % ==> tolerated
% 
% % (likely both producing) 'msa[m]'	'RM00908'	'RM01608'	'P/C'	0.217099278516508	1x1 table	'akg[m] + ala-B[m]  <=> glu-L[m] + msa[m] '	'nad[m] + 3hpp[m]  -> nadh[m] + h[m] + msa[m] '	4
% metFit(strcmp(metFit.branchMets,'msa[m]'),:)
% a(ismember(a.rxn, {'RM00908'	'RM01608'}),:)
% % ==> tolerated
% 
% % 'utp[c]'	'RC00416'	'BIO0003'	'P/C'	0.230503569066053	1x14 table	'h[c] + utp[c] + acgam1p[c]  <=> ppi[c] + uacgam[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	8
% metFit(strcmp(metFit.branchMets,'utp[c]'),:)
% a(ismember(a.rxn, {'RC00416'	'BIO0003'}),:)
% % ==> tolerated
% 
% % 'glu-L[c]'	'RC07396'	'RC00573'	'P/C'	0.243073574138845	1x1 table	'glu-L[c] + 2kmb[c]  -> akg[c] + met-L[c] '	'atp[c] + h2o[c] + gln-L[c] + utp[c]  -> adp[c] + 2 h[c] + pi[c] + glu-L[c] + ctp[c] '	33
% % 'glu-L[c]'	'RC00734'	'RC00573'	'P/C'	0.243073574138845	2x1 table	'akg[c] + tyr-L[c]  <=> glu-L[c] + 34hpp[c] '	'atp[c] + h2o[c] + gln-L[c] + utp[c]  -> adp[c] + 2 h[c] + pi[c] + glu-L[c] + ctp[c] '	33
% % 'glu-L[c]'	'RC00694'	'RC00573'	'P/C'	0.243073574138845	2x1 table	'akg[c] + phe-L[c]  <=> glu-L[c] + phpyr[c] '	'atp[c] + h2o[c] + gln-L[c] + utp[c]  -> adp[c] + 2 h[c] + pi[c] + glu-L[c] + ctp[c] '	33
% % 'glu-L[m]'	'RM01648'	'RM00258'	'P/C'	0.231206999741968	1x1 table	'akg[m] + 4abut[m]  -> sucsal[m] + glu-L[m] '	'akg[m] + ala-L[m]  <=> pyr[m] + glu-L[m] '	27
% % 'glu-L[m]'	'RM04188'	'RM00258'	'P/C'	0.231206999741968	1x1 table	'glu-L[m] + 2mop[m]  -> akg[m] + 3aib[m] '	'akg[m] + ala-L[m]  <=> pyr[m] + glu-L[m] '	27
% % 'glu-L[m]'	'RM00908'	'RM00258'	'P/C'	0.231206999741968	1x1 table	'akg[m] + ala-B[m]  <=> glu-L[m] + msa[m] '	'akg[m] + ala-L[m]  <=> pyr[m] + glu-L[m] '	27
% metFit(strcmp(metFit.branchMets,'glu-L[c]'),:)
% a(ismember(a.rxn, {'RC07396' 'RC00734' 'RC00694'	'RC00573'}),:)
% metFit(strcmp(metFit.branchMets,'glu-L[m]'),:)
% a(ismember(a.rxn, {'RM01648' 'RM00258' 'RM04188'	'RM00908'}),:)
% % ==> partially balanced and mostly tolerated
% 
% % likely flux coupled but not major flux wiring (check the fitting result
% % on these) (big flux --> small flux)
% % 'amp[c]'	'RC00127'	'RC00185'	'P/C'	0.277577285608170	1x2 table	'atp[c] + amp[c]  <=> 2 adp[c] '	'atp[c] + adn[c]  -> adp[c] + h[c] + amp[c] '	69
% % 'amp[c]'	'RC00127'	'RC05578'	'P/C'	0.255730102681089	1x1 table	'atp[c] + amp[c]  <=> 2 adp[c] '	'atp[c] + glu-L[c] + trnaglu[c]  -> amp[c] + ppi[c] + glutrna[c] '	69
% % 'amp[c]'	'RC00127'	'RC03652'	'P/C'	0.227106926045622	1x1 table	'atp[c] + amp[c]  <=> 2 adp[c] '	'atp[c] + gln-L[c] + trnagln[c]  -> amp[c] + ppi[c] + glntrna[c] '	69
% % 'amp[c]'	'RC00127'	'RC03658'	'P/C'	0.305276841143129	1x1 table	'atp[c] + amp[c]  <=> 2 adp[c] '	'atp[c] + lys-L[c] + trnalys[c]  -> amp[c] + ppi[c] + lystrna[c] '	69
% % 'amp[c]'	'RC00127'	'RC00181'	'P/C'	0.252367952983668	1x1 table	'atp[c] + amp[c]  <=> 2 adp[c] '	'h[c] + amp[c] + h2o[c]  -> nh4[c] + imp[c] '	69
% % 'amp[c]'	'RC05578'	'RC00181'	'P/C'	0.221284644491998	1x1 table	'atp[c] + glu-L[c] + trnaglu[c]  -> amp[c] + ppi[c] + glutrna[c] '	'h[c] + amp[c] + h2o[c]  -> nh4[c] + imp[c] '	69
% % 'amp[c]'	'RC03658'	'RC00181'	'P/C'	0.287638620806091	1x1 table	'atp[c] + lys-L[c] + trnalys[c]  -> amp[c] + ppi[c] + lystrna[c] '	'h[c] + amp[c] + h2o[c]  -> nh4[c] + imp[c] '	69
% metFit(strcmp(metFit.branchMets,'amp[c]'),:)
% a(ismember(a.rxn, {'RC00127' 'RC00185' 'RC05578'	'RC03652' 'RC03658', 'RC00181'}),:)
% % ==> tolaterated with huge misfit; but this individual constraint can be
% % fitted if fitting alone
% 
% % 'udpg[c]'	'RC01005'	'RC00289'	'P/C'	0.811145061649767	1x1 table	'dolp[c] + udpg[c]  -> udp[c] + dolglcp[c] '	'g1p[c] + h[c] + utp[c]  <=> ppi[c] + udpg[c] '	10
% metFit(strcmp(metFit.branchMets,'udpg[c]'),:)
% a(ismember(a.rxn, {'RC01005' 'RC00289'}),:)
% % ==> tolerated
% 
% % 'amp[m]'	'RM03658'	'RM00127'	'P/C'	0.305276841143129	1x1 table	'atp[m] + lys-L[m] + trnalys[m]  -> ppi[m] + amp[m] + lystrna[m] '	'atp[m] + amp[m]  -> 2 adp[m] '	23
% metFit(strcmp(metFit.branchMets,'amp[m]'),:)
% a(ismember(a.rxn, {'RM03658' 'RM00127'}),:)
% % ==> tolaterated with huge misfit; strong conflict, even individual alone
% % cannot be fitted; this is well expected as this has to be the case with
% % biomass production
% 
% 
% 
% 
% % interesting to check:
% % (1) whether the g3p constraint eventually lead to the prediction of large
% % flux in gluconeogeneiss 
% % 'g3p[c]'	'RC01070'	'RC01066'	'P/C'	0.282639515581195	1x2 table	'fdp[c]  <=> dhap[c] + g3p[c] '	'2dr5p[c]  <=> g3p[c] + acald[c] '
% % 'g3p[c]'	'RC01641'	'RC01827'	'P/C'	0.597500826723647	1x1 table	'g3p[c] + s7p[c]  <=> xu5p-D[c] + r5p[c] '	'g3p[c] + s7p[c]  <=> f6p[c] + e4p[c] '
% % 'g3p[c]'	'RC01641'	'RC01066'	'P/C'	0.223896084077802	1x2 table	'g3p[c] + s7p[c]  <=> xu5p-D[c] + r5p[c] '	'2dr5p[c]  <=> g3p[c] + acald[c] '
% % 'g3p[c]'	'RC01830'	'RC01827'	'P/C'	0.597500826723647	1x1 table	'f6p[c] + g3p[c]  <=> xu5p-D[c] + e4p[c] '	'g3p[c] + s7p[c]  <=> f6p[c] + e4p[c] '
% % 'g3p[c]'	'RC01830'	'RC01066'	'P/C'	0.223896084077802	1x2 table	'f6p[c] + g3p[c]  <=> xu5p-D[c] + e4p[c] '	'2dr5p[c]  <=> g3p[c] + acald[c] '
% metFit(strcmp(metFit.branchMets,'g3p[c]'),:)
% a(ismember(a.rxn, {'RC01070' 'RC01066','RC01641','RC01827','RC01830'}),:)
% % ==> gluconeogenesis flux is predicted; not large but bounded. looks
% % interesting
% 
% 
% % (2) if gluconeogenesis flux is correctly maintained and if RNA syn is
% % disrupted 
% % 'gtp[c]'	'RC00431'	'BIO0003'	'P/C'	0.429026643869237	2x14 table	'gtp[c] + oaa[c]  <=> pep[c] + gdp[c] + co2[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	15
% % 'gtp[c]'	'RC00330'	'BIO0003'	'P/C'	0.322909515508412	1x14 table	'atp[c] + gdp[c]  <=> adp[c] + gtp[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	15
% % 'oaa[c]'	'RC00431'	'RC00352'	'P/C'	0.211766335917045	2x2 table	'gtp[c] + oaa[c]  <=> pep[c] + gdp[c] + co2[c] '	'atp[c] + coa[c] + cit[c]  -> adp[c] + pi[c] + oaa[c] + accoa[c] '
% % 'oaa[c]'	'RC00431'	'RC00355'	'P/C'	0.229297644288476	2x1 table	'gtp[c] + oaa[c]  <=> pep[c] + gdp[c] + co2[c] '	'akg[c] + asp-L[c]  <=> oaa[c] + glu-L[c] '
% % 'gdp[c]'	'RC00431'	'RC00328'	'P/C'	0.348251665154161	2x1 table	'gtp[c] + oaa[c]  <=> pep[c] + gdp[c] + co2[c] '	'h2o[c] + gdp[c]  -> h[c] + pi[c] + gmp[c] '	18
% metFit(strcmp(metFit.branchMets,'g3p[c]'),:)
% a(ismember(a.rxn, {'RC01070' 'RC01066','RC01641','RC01827','RC01830'}),:)
% % ==> gluconeogenesis flux is predicted; not large but bounded. looks
% % interesting
% % however, still created a large misfit; seems like there are some conclict
% 
% 
% % (3) pyruvate flux is coupled to alanine (although very weak) how would
% % this influence the flux preditcion
% % 'pyr[m]'	'RM00209'	'RM00258'	'P/C'	0.207368207147373	4x1 table	'pyr[m] + coa[m] + nad[m]  -> accoa[m] + co2[m] + nadh[m] '	'akg[m] + ala-L[m]  <=> pyr[m] + glu-L[m] '	12
% metFit(strcmp(metFit.branchMets,'pyr[m]'),:)
% a(ismember(a.rxn, {'RM00209' 'RM00258'}),:)
% % ala was largely converted to pyr without the constriant; no sig influence
% % and leaves a significant misfit
% 
% 
% % (4) accoa largely coupled  to mito beta and PDH; check how this influence
% % the aa/shunt sourced accoa
% % 'accoa[m]'	'RM00209'	'RM00351'	'P/C'	0.340519083595258	4x1 table	'pyr[m] + coa[m] + nad[m]  -> accoa[m] + co2[m] + nadh[m] '	'accoa[m] + h2o[m] + oaa[m]  -> coa[m] + h[m] + cit[m] '	35
% % 'accoa[m]'	'RM00351'	'RMC0088'	'P/C'	0.257813234766866	1x4 table	'accoa[m] + h2o[m] + oaa[m]  -> coa[m] + h[m] + cit[m] '	'4 coa[m] + 4 nad[m] + 4 h2o[m] + 3 fad[m] + fa16p1n7coa[m]  -> 4 accoa[m] + 4 nadh[m] + 4 h[m] + 3 fadh2[m] + occoa[m] '	35
% % and many others
% % the branching points will be integrated in the iMAT as one more step of
% % optimization;
% metFit(strcmp(metFit.branchMets,'accoa[m]'),:)
% a(ismember(a.rxn, {'RM00209' 'RM00351','RMC0088'}),:)
% % increased the contribution of pyruvate; no suprise - a large misfit still
% % remains; didnt influence the high shunt contribution (but now equal to
% % that of pyr)
% 
% 
% % (5) cytosolic accoa largely coupled to histone acetylation; this is very
% % surprisingly predicted by iMAT++ itself; double check if it changes 
% % 'accoa[c]'	'RC01978'	'RC03552'	'P/C'	0.363271878610942	1x1 table	'h2o[c] + accoa[c] + aacoa[c]  -> h[c] + coa[c] + hmgcoa[c] '	'accoa[c] + Nlys[c]  <=> h[c] + coa[c] + Naclys[c] '	40
% % 'accoa[c]'	'RCC0019'	'RC03552'	'P/C'	0.363704525616254	1x1 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'accoa[c] + Nlys[c]  <=> h[c] + coa[c] + Naclys[c] '	40
% % 'accoa[c]'	'RC02058'	'RC03552'	'P/C'	0.270926143382104	2x1 table	'accoa[c] + gam6p[c]  -> h[c] + coa[c] + acgam6p[c] '	'accoa[c] + Nlys[c]  <=> h[c] + coa[c] + Naclys[c] '	40
% metFit(strcmp(metFit.branchMets,'accoa[c]'),:)
% a(ismember(a.rxn, {'RC01978' 'RC03552','RCC0019','RC02058'}),:)
% % --> no change 
% 
% % (6) cytC syn is coupled with ETC reduction of cytc; see if it has a negative
% % impact on flux (should be correctly buffered)
% % 'focytC[m]'	'RM00197'	'RM00081'	'P/C'	0.221484378466578	1x4 table	'lac-D[m] + 2 ficytC[m]  -> pyr[m] + 2 h[m] + 2 focytC[m] '	'8 h[m] + 4 focytC[m] + o2[m]  -> 4 h[c] + 2 h2o[m] + 4 ficytC[m] '	7
% % 'focytC[m]'	'RM02161'	'RM00081'	'P/C'	0.701882771479615	8x4 table	'2 h[m] + qh2[m] + 2 ficytC[m]  -> 4 h[c] + q[m] + 2 focytC[m] '	'8 h[m] + 4 focytC[m] + o2[m]  -> 4 h[c] + 2 h2o[m] + 4 ficytC[m] '	7
% % 'focytC[m]'	'RM00081'	'RM02480'	'P/C'	0.315962093857402	4x2 table	'8 h[m] + 4 focytC[m] + o2[m]  -> 4 h[c] + 2 h2o[m] + 4 ficytC[m] '	'apocytc[m] + pheme[m]  -> focytC[m] '	7
% metFit(strcmp(metFit.branchMets,'focytC[m]'),:)
% a(ismember(a.rxn, {'RM00197' 'RM00081','RM02161','RM00081','RM02480'}),:)
% % --> easily and fully buffered bc the syn flux is very small 
% 
% 
% % (7) methionine influx: the weak corr between mtrr-1 and sams-1 adds
% % coupling between folate-met influx and met/sam cycle. This may strongly
% % change the flux to cyclic pattern that will be problemmatic; see how it
% % will go! (to fix it, we may have to identify all no-measure pairs for all
% % metabolites considered here, and always include these no-measure pairs in
% % the total metabolite flux optimalziation; this will add the dietary met
% % flux back into the influxes considered for met-L) But first see if it can
% % be buffered.
% % 'met-L[c]'	'RC00946'	'RC00177'	'P/C'	0.293263547288988	3x4 table	'hcys-L[c] + 5mthf[c]  -> h[c] + met-L[c] + thf[c] '	'atp[c] + h2o[c] + met-L[c]  -> pi[c] + ppi[c] + amet[c] '	8
% metFit(strcmp(metFit.branchMets,'met-L[c]'),:)
% a(ismember(a.rxn, {'RC00946' 'RC00177'}),:)
% % fully buffered; the other constaints look very strong such that the met-L
% % couling is bearly fitted
% 
% 
% % (8) ketone supply of hmgcoa (see if this will influence flux prediction)
% % 'hmgcoa[c]'	'RC01978'	'RC02082'	'P/C'	0.203254558460800	1x1 table	'h2o[c] + accoa[c] + aacoa[c]  -> h[c] + coa[c] + hmgcoa[c] '	'2 h[c] + 2 nadph[c] + hmgcoa[c]  -> 2 nadp[c] + coa[c] + mev-R[c] '	3
% metFit(strcmp(metFit.branchMets,'hmgcoa[c]'),:)
% a(ismember(a.rxn, {'RC01978' 'RC02082'}),:)
% % ketone was already the major supply; the fitting only increase it a
% % little and it was still not fully fitted
% 
% 
% % (9) are the many weak couplings around nucleic acid and DNA/RNA syn meaningful?
% % 'ctp[c]'	'RC01799'	'RC00570'	'P/C'	0.259224683868232	1x1 table	'h[c] + pa_pl[c] + ctp[c]  -> ppi[c] + cdpdag[c] '	'atp[c] + cdp[c]  <=> adp[c] + ctp[c] '	10
% % 'ctp[c]'	'RC00571'	'BIO0003'	'P/C'	0.240322150587729	1x14 table	'atp[c] + nh4[c] + utp[c]  -> adp[c] + 2 h[c] + pi[c] + ctp[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	10
% % 'ctp[c]'	'RC00573'	'BIO0003'	'P/C'	0.240322150587729	1x14 table	'atp[c] + h2o[c] + gln-L[c] + utp[c]  -> adp[c] + 2 h[c] + pi[c] + glu-L[c] + ctp[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	10
% % 'ctp[c]'	'RC00570'	'BIO0003'	'P/C'	0.322909515508412	1x14 table	'atp[c] + cdp[c]  <=> adp[c] + ctp[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	10
% % 'utp[c]'	'RC00156'	'BIO0003'	'P/C'	0.322909515508412	1x14 table	'atp[c] + udp[c]  <=> adp[c] + utp[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	8
% % 'utp[c]'	'RC00416'	'BIO0003'	'P/C'	0.230503569066053	1x14 table	'h[c] + utp[c] + acgam1p[c]  <=> ppi[c] + uacgam[c] '	'0.639534 atp[c] + 0.406032 h2o[c] + 0.233503 gtp[c] + 0.228426 ctp[c] + 0.304569 utp[c]  -> 0.406032 adp[c] + 0.406032 h[c] + 0.406032 pi[c] + ppi[c] + RNA[c] '	8
% % 'dctp[c]'	'RC02326'	'BIO0004'	'P/C'	0.309397287122737	1x3 table	'atp[c] + dcdp[c]  -> adp[c] + dctp[c] '	'0.368 atp[c] + 0.368 h2o[c] + 0.325 datp[c] + 0.175 dgtp[c] + 0.175 dctp[c] + 0.325 dttp[c]  -> 0.368 adp[c] + 0.368 h[c] + 0.368 pi[c] + ppi[c] + DNA[c] '	2
% % 'dttp[c]'	'RC02093'	'BIO0004'	'P/C'	0.309397287122737	1x3 table	'atp[c] + dtdp[c]  <=> adp[c] + dttp[c] '	'0.368 atp[c] + 0.368 h2o[c] + 0.325 datp[c] + 0.175 dgtp[c] + 0.175 dctp[c] + 0.325 dttp[c]  -> 0.368 adp[c] + 0.368 h[c] + 0.368 pi[c] + ppi[c] + DNA[c] '	3
% metFit(strcmp(metFit.branchMets,'ctp[c]'),:)
% a(ismember(a.rxn, {'RC01799' 'RC00570','RC00571','BIO0003','RC00573'}),:)
% % --> only minor change of flux bc it is already quite consistent with
% % original flux distribution
% metFit(strcmp(metFit.branchMets,'utp[c]'),:)
% a(ismember(a.rxn, {'RC00156' 'BIO0003','RC00416'}),:)
% % flux didnt change with, leaving a sig (but not large) misfit
% metFit(strcmp(metFit.branchMets,'dctp[c]'),:)
% a(ismember(a.rxn, {'RC02326' 'BIO0004'}),:)
% metFit(strcmp(metFit.branchMets,'dttp[c]'),:)
% a(ismember(a.rxn, {'RC02093' 'BIO0004'}),:)
% % easily fitted since originally was already well consistent
% % --> generally these couplings were well consistent with the originally
% % distribution and most can be easily fitted.
% 
% 
% % (10) dag: not likely correct but see if it is buffered or actually
% % meaningful
% % 'dag_pl[c]'	'RC02239'	'RC02057'	'P/C'	0.221548350905133	2x1 table	'h2o[c] + pa_pl[c]  -> pi[c] + dag_pl[c] '	'dag_pl[c] + cdpea[c]  <=> h[c] + cmp[c] + pe[c] '	7
% % 'dag_pl[c]'	'RC02057'	'RC03435'	'P/C'	0.219641996542157	1x3 table	'dag_pl[c] + cdpea[c]  <=> h[c] + cmp[c] + pe[c] '	'h2o[c] + pail45p[c]  -> h[c] + dag_pl[c] + mi145p[c] '	7
% % 'dag_pl[c]'	'RC02057'	'RC03332'	'P/C'	0.219641996542157	1x3 table	'dag_pl[c] + cdpea[c]  <=> h[c] + cmp[c] + pe[c] '	'h2o[c] + pail[c]  -> h[c] + dag_pl[c] + mi1p-D[c] '	7
% metFit(strcmp(metFit.branchMets,'dag_pl[c]'),:)
% a(ismember(a.rxn, {'RC02239' 'RC02057','RC03435','RC03332'}),:)
% % bearly fit at all and was fully buffered. 
% 
% 
% % (11) the source of ethamp; coupling is likely wrong but interesting to
% % see how the model reacts to it
% % 'ethamp[c]'	'RC02037'	'RC02464'	'P/C'	0.318065858870702	1x1 table	'amet[c] + ethamp[c]  -> h[c] + ahcys[c] + methamp[c] '	'sph1p_ce[c]  -> ethamp[c] + m15hxdcal[c] '	7
% metFit(strcmp(metFit.branchMets,'ethamp[c]'),:)
% a(ismember(a.rxn, {'RC02037' 'RC02464'}),:)
% % fitted a little bit and was buffered 
% 
% 
% % (12) coupling between peroxi oxi and the FA syn is captured. see how
% % model reacts
% % 'pmtcoa[c]'	'RC07758'	'RCC0126'	'P/C'	0.399842366054898	2x4 table	'h[c] + pmtcoa[c] + malcoa[c]  -> co2[c] + coa[c] + 3oodcoa[c] '	'h2o[c] + nad[c] + coa[c] + o2[c] + stcoa[c]  -> h[c] + nadh[c] + accoa[c] + h2o2[c] + pmtcoa[c] '	11
% % 'pmtcoa[c]'	'RCC0165'	'RCC0126'	'P/C'	0.692051466068343	1x4 table	'o2[c] + pmtcoa[c]  -> h2o2[c] + hdd2coa[c] '	'h2o[c] + nad[c] + coa[c] + o2[c] + stcoa[c]  -> h[c] + nadh[c] + accoa[c] + h2o2[c] + pmtcoa[c] '	11
% % 'stcoa[c]'	'RCC0124'	'RCC0125'	'P/C'	0.412245389710732	4x4 table	'9 h[c] + 6 nadph[c] + 3 malcoa[c] + stcoa[c]  -> 3 h2o[c] + 3 co2[c] + 6 nadp[c] + 3 coa[c] + lgnccoa[c] '	'3 h2o[c] + 3 nad[c] + 3 coa[c] + 3 o2[c] + lgnccoa[c]  -> 3 h[c] + 3 nadh[c] + 3 accoa[c] + 3 h2o2[c] + stcoa[c] '	9
% % 'dlnlcgcoa[c]'	'RCC0064'	'RCC0093'	'P/C'	0.420804806043524	5x4 table	'3 h[c] + 2 nadph[c] + malcoa[c] + lnlncgcoa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + dlnlcgcoa[c] '	'2 h2o[c] + 2 nad[c] + 2 coa[c] + 2 o2[c] + dlnlcgcoa[c]  -> 2 h[c] + 2 nadh[c] + 2 accoa[c] + 2 h2o2[c] + fa16p3n6coa[c] '	7
% % 'vacccoa[c]'	'RCC0070'	'RCC0090'	'P/C'	0.420804806043524	13x4 table	'3 h[c] + 2 nadph[c] + malcoa[c] + fa16p1n7coa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + vacccoa[c] '	'5 h2o[c] + 5 nad[c] + 5 coa[c] + 4 o2[c] + vacccoa[c]  -> 5 h[c] + 5 nadh[c] + 5 accoa[c] + 4 h2o2[c] + occoa[c] '	6
% % 'lgnccoa[c]'	'RCC0124'	'RCC0125'	'P/C'	0.412245389710732	4x4 table	'9 h[c] + 6 nadph[c] + 3 malcoa[c] + stcoa[c]  -> 3 h2o[c] + 3 co2[c] + 6 nadp[c] + 3 coa[c] + lgnccoa[c] '	'3 h2o[c] + 3 nad[c] + 3 coa[c] + 3 o2[c] + lgnccoa[c]  -> 3 h[c] + 3 nadh[c] + 3 accoa[c] + 3 h2o2[c] + stcoa[c] '	5
% metFit(strcmp(metFit.branchMets,'dlnlcgcoa[c]'),:)
% a(ismember(a.rxn, {'RCC0064' 'RCC0093'}),:)
% metFit(strcmp(metFit.branchMets,'pmtcoa[c]'),:)
% a(ismember(a.rxn, {'RC07758' 'RCC0126','RCC0165'}),:)
% metFit(strcmp(metFit.branchMets,'stcoa[c]'),:)
% a(ismember(a.rxn, {'RCC0124' 'RCC0125'}),:)
% metFit(strcmp(metFit.branchMets,'vacccoa[c]'),:)
% a(ismember(a.rxn, {'RCC0070' 'RCC0090'}),:)
% metFit(strcmp(metFit.branchMets,'lgnccoa[c]'),:)
% a(ismember(a.rxn, {'RCC0124' 'RCC0125'}),:)
% % -- not fitted, leaving a large misfit
% % 'eicostetcoa[c]'	'RCC0063'	'RCC0096'	'P/C'	0.420804806043524	5x4 table	'3 h[c] + 2 nadph[c] + malcoa[c] + strdnccoa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + eicostetcoa[c] '	'2 h2o[c] + 2 nad[c] + 2 coa[c] + 2 o2[c] + eicostetcoa[c]  -> 2 h[c] + 2 nadh[c] + 2 accoa[c] + 2 h2o2[c] + fa16p4n3coa[c] '	7
% metFit(strcmp(metFit.branchMets,'eicostetcoa[c]'),:)
% a(ismember(a.rxn, {'RCC0063' 'RCC0096'}),:)
% % --> fit by pushing to zero flux
% % these couplings were largely NOT fitted at all; The bulk FA such as
% % pmtcoa was not degraded before using for elongation and seems cannot be
% % even fitted so.
% % --> [IMPORTANT] we found peroxi oxi is active but is not from the
% % orxidation of usable FA for wormBiomass; the FA is already in shortage
% % such that it has to be synthesized therefore, the model cannot oxidaze
% % the FA to fit the coupling; this may be changed by loosing the yeild
% % constraint (but may not if the minLow push it to the limit again). This
% % looks like an interesting observation
% 
% 
% 
% % (13) malcoa is implied to be coupled with very long chain FA elongation
% % instead of FASN. see how model reacts
% % 'malcoa[c]'	'RCC0019'	'RCC0063'	'P/C'	0.674902391693208	1x5 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'3 h[c] + 2 nadph[c] + malcoa[c] + strdnccoa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + eicostetcoa[c] '	13
% % 'malcoa[c]'	'RCC0019'	'RCC0064'	'P/C'	0.674902391693208	1x5 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'3 h[c] + 2 nadph[c] + malcoa[c] + lnlncgcoa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + dlnlcgcoa[c] '	13
% % 'malcoa[c]'	'RCC0019'	'RCC0070'	'P/C'	0.761844021227047	1x13 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'3 h[c] + 2 nadph[c] + malcoa[c] + fa16p1n7coa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + vacccoa[c] '	13
% % 'malcoa[c]'	'RCC0019'	'RCC0072'	'P/C'	0.761844021227047	1x5 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'3 h[c] + 2 nadph[c] + malcoa[c] + fa13p0iso[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + fa15p0iso[c] '	13
% % 'malcoa[c]'	'RCC0019'	'RCC0073'	'P/C'	0.761844021227047	1x7 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'3 h[c] + 2 nadph[c] + malcoa[c] + fa15p0iso[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + fa17p0iso[c] '	13
% % 'malcoa[c]'	'RCC0019'	'RCC0124'	'P/C'	0.674902391693208	1x4 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'9 h[c] + 6 nadph[c] + 3 malcoa[c] + stcoa[c]  -> 3 h2o[c] + 3 co2[c] + 6 nadp[c] + 3 coa[c] + lgnccoa[c] '	13
% % 'malcoa[c]'	'RCC0019'	'RCC0210'	'P/C'	0.674902391693208	1x4 table	'atp[c] + accoa[c] + hco3[c]  -> adp[c] + h[c] + pi[c] + malcoa[c] '	'6 h[c] + 4 nadph[c] + 2 malcoa[c] + lgnccoa[c]  -> 2 h2o[c] + 2 co2[c] + 4 nadp[c] + 2 coa[c] + fa28p0coa[c] '	13
% metFit(strcmp(metFit.branchMets,'malcoa[c]'),:)
% a(ismember(a.rxn, {'RCC0019' 'RCC0063','RCC0064','RCC0070','RCC0072','RCC0073','RCC0124','RCC0210'}),:)
% % [interesting] the model entirely igrored this constaint and left a large
% % unfit value; it cannot be fitted even on itself; it is constrianed by the
% % yeild constraint, and a reasonable lowering is not solving the problem
% % (0.65 --> 0.3 will change the unfit proportion from 90% to 50% but still
% % high). Therefore, this is not likely a real flux wiring (i.e., hypothesis
% % is wrong). 
% 
% 
% % (14) nucleitide interconversion is constrained with a few gene
% % similarities. and also related to DNA syn. See how it influences the
% % model 
% % 'ins[c]'	'RC01863'	'RC01560'	'P/C'	0.508164616176076	1x1 table	'pi[c] + ins[c]  <=> r1p[c] + hxan[c] '	'h[c] + h2o[c] + adn[c]  -> nh4[c] + ins[c] '	4
% metFit(strcmp(metFit.branchMets,'ins[c]'),:)
% a(ismember(a.rxn, {'RC01863' 'RC01560'}),:)
% % --> unfitted
% % 'din[c]'	'RC02748'	'RC02556'	'P/C'	0.508164616176076	1x1 table	'pi[c] + din[c]  <=> hxan[c] + 2dr1p[c] '	'h[c] + h2o[c] + dad-2[c]  -> nh4[c] + din[c] '	3
% metFit(strcmp(metFit.branchMets,'din[c]'),:)
% a(ismember(a.rxn, {'RC02748' 'RC02556'}),:)
% % --> easy fit, small flux
% % 'dad-2[c]'	'RC02557'	'RC02556'	'P/C'	0.508164616176076	1x1 table	'pi[c] + dad-2[c]  <=> ade[c] + 2dr1p[c] '	'h[c] + h2o[c] + dad-2[c]  -> nh4[c] + din[c] '	4
% metFit(strcmp(metFit.branchMets,'dad-2[c]'),:)
% a(ismember(a.rxn, {'RC02557' 'RC02556'}),:)
% % --> unfitted
% % 'datp[c]'	'RC01137'	'BIO0004'	'P/C'	0.309397287122737	1x3 table	'atp[c] + dadp[c]  <=> adp[c] + datp[c] '	'0.368 atp[c] + 0.368 h2o[c] + 0.325 datp[c] + 0.175 dgtp[c] + 0.175 dctp[c] + 0.325 dttp[c]  -> 0.368 adp[c] + 0.368 h[c] + 0.368 pi[c] + ppi[c] + DNA[c] '	3
% metFit(strcmp(metFit.branchMets,'datp[c]'),:)
% a(ismember(a.rxn, {'RC01137' 'BIO0004'}),:)
% % --> easy fit
% % 'xmp[c]'	'RC01130'	'RC01230'	'P/C'	0.239644988901667	1x1 table	'h2o[c] + nad[c] + imp[c]  -> h[c] + nadh[c] + xmp[c] '	'atp[c] + nh4[c] + xmp[c]  -> 2 h[c] + amp[c] + ppi[c] + gmp[c] '	6
% % 'xmp[c]'	'RC01130'	'RC01231'	'P/C'	0.239644988901667	1x1 table	'h2o[c] + nad[c] + imp[c]  -> h[c] + nadh[c] + xmp[c] '	'atp[c] + h2o[c] + gln-L[c] + xmp[c]  -> 2 h[c] + amp[c] + ppi[c] + glu-L[c] + gmp[c] '	6
% metFit(strcmp(metFit.branchMets,'xmp[c]'),:)
% a(ismember(a.rxn, {'RC01130' 'RC01230','RC01231'}),:)
% % --> easy fit, small flux
% % 'dgtp[c]'	'RC01857'	'BIO0004'	'P/C'	0.309397287122737	1x3 table	'atp[c] + dgdp[c]  -> adp[c] + dgtp[c] '	'0.368 atp[c] + 0.368 h2o[c] + 0.325 datp[c] + 0.175 dgtp[c] + 0.175 dctp[c] + 0.325 dttp[c]  -> 0.368 adp[c] + 0.368 h[c] + 0.368 pi[c] + ppi[c] + DNA[c] '	3
% metFit(strcmp(metFit.branchMets,'dgtp[c]'),:)
% a(ismember(a.rxn, {'RC01857' 'BIO0004'}),:)
% % --> easy fit
% % --> these couplings are generally not constrained by other things so they
% % are generally fitted at minimal level (e.g., one of the flux was pushed
% % to epsilon in the total flux minimiazation), or as needed for biomass
% % assembly
% 
% 
% % good predictions 
% % --> big fluxes can be not fully fitted, however, already consistent 
% % --> most of them was originally very consistent
% % TCA cycle 
% % 'oaa[m]'	'RM00351'	'RM00342'	'P/C'	0.411474101112489	1x1 table	'accoa[m] + h2o[m] + oaa[m]  -> coa[m] + h[m] + cit[m] '	'nad[m] + mal-L[m]  <=> nadh[m] + h[m] + oaa[m] '	4
% % 'cit[m]'	'RM00351'	'RMC0001'	'P/C'	0.404391889848461	1x1 table	'accoa[m] + h2o[m] + oaa[m]  -> coa[m] + h[m] + cit[m] '	'cit[m]  <=> icit[m] '	5
% % 'fum[m]'	'RM02164'	'RM01082'	'P/C'	0.325678505431294	5x1 table	'succ[m] + q[m]  -> qh2[m] + fum[m] '	'mal-L[m]  <=> h2o[m] + fum[m] '	8
% metFit(strcmp(metFit.branchMets,'oaa[m]'),:)
% a(ismember(a.rxn, {'RM00351' 'RM00342'}),:)
% metFit(strcmp(metFit.branchMets,'cit[m]'),:)
% a(ismember(a.rxn, {'RM00351' 'RMC0001'}),:)
% metFit(strcmp(metFit.branchMets,'fum[m]'),:)
% a(ismember(a.rxn, {'RM02164' 'RM01082'}),:)
% % --> one heavily constrained; cannot be fitted (oaa) but origianlly relatively
% % consistent; other two was fitted perfectly
% 
% % ETC
% % 'q[m]'	'RM02164'	'RM02161'	'P/C'	0.805107311734685	5x8 table	'succ[m] + q[m]  -> qh2[m] + fum[m] '	'2 h[m] + qh2[m] + 2 ficytC[m]  -> 4 h[c] + q[m] + 2 focytC[m] '	12
% % 'q[m]'	'RMC0006'	'RM02161'	'P/C'	0.725574680926681	29x8 table	'nadh[m] + 5 h[m] + q[m]  -> 4 h[c] + nad[m] + qh2[m] '	'2 h[m] + qh2[m] + 2 ficytC[m]  -> 4 h[c] + q[m] + 2 focytC[m] '	12
% % 'q[m]'	'RM02161'	'RMC0007'	'P/C'	0.725574680926681	8x29 table	'2 h[m] + qh2[m] + 2 ficytC[m]  -> 4 h[c] + q[m] + 2 focytC[m] '	'h[c] + nadh[m] + q[m]  -> nad[m] + qh2[m] '	12
% % 'q[m]'	'RM02161'	'RMC0008'	'P/C'	0.725574680926681	8x29 table	'2 h[m] + qh2[m] + 2 ficytC[m]  -> 4 h[c] + q[m] + 2 focytC[m] '	'nadh[m] + h[m] + q[m]  -> nad[m] + qh2[m] '	12
% metFit(strcmp(metFit.branchMets,'q[m]'),:)
% a(ismember(a.rxn, {'RM02164' 'RM02161','RMC0006','RMC0007','RMC0008'}),:)
% % --> not fully fitted but originallly already consistent; seems RMC0007
% % and RMC0008 is constrained by minlow such that it cannot be fitted; or it
% % is something else; not sure
% 
% % PPP
% % '6pgl[c]'	'RC02736'	'RC02035'	'P/C'	0.310545598955383	1x2 table	'g6p-B[c] + nadp[c]  -> h[c] + nadph[c] + 6pgl[c] '	'h2o[c] + 6pgl[c]  -> h[c] + 6pgc[c] '	2
% % '6pgc[c]'	'RC02035'	'RC01528'	'P/C'	0.236161974244792	2x1 table	'h2o[c] + 6pgl[c]  -> h[c] + 6pgc[c] '	'nadp[c] + 6pgc[c]  -> co2[c] + nadph[c] + ru5p-D[c] '	3
% metFit(strcmp(metFit.branchMets,'6pgl[c]'),:)
% a(ismember(a.rxn, {'RC02736' 'RC02035','RC01528'}),:)
% % fully fitted
% 
% % and others
% % chitin syn (reveals hxk-1 coupling)
% % 'gam6p[c]'	'RC00768'	'RC02058'	'P/C'	0.612605434668005	2x2 table	'f6p[c] + gln-L[c]  -> glu-L[c] + gam6p[c] '	'accoa[c] + gam6p[c]  -> h[c] + coa[c] + acgam6p[c] '	6
% % 'gam6p[c]'	'RC01961'	'RC02058'	'P/C'	0.613175977760319	1x2 table	'atp[c] + gam[c]  -> adp[c] + h[c] + gam6p[c] '	'accoa[c] + gam6p[c]  -> h[c] + coa[c] + acgam6p[c] '	6
% metFit(strcmp(metFit.branchMets,'gam6p[c]'),:)
% a(ismember(a.rxn, {'RC00768' 'RC02058','RC01961'}),:)
% % --> largely fitted
% 
% % folate cycle 
% % '5mthf[c]'	'RC00946'	'RC01224'	'P/C'	0.609931610436158	3x2 table	'hcys-L[c] + 5mthf[c]  -> h[c] + met-L[c] + thf[c] '	'2 h[c] + nadph[c] + mlthf[c]  -> nadp[c] + 5mthf[c] '	4
% % '5mthf[c]'	'RC00946'	'RC07168'	'P/C'	0.609931610436158	3x2 table	'hcys-L[c] + 5mthf[c]  -> h[c] + met-L[c] + thf[c] '	'2 h[c] + nadh[c] + mlthf[c]  -> nad[c] + 5mthf[c] '	4
% metFit(strcmp(metFit.branchMets,'5mthf[c]'),:)
% a(ismember(a.rxn, {'RC00946' 'RC01224','RC07168'}),:)
% % (almost) fully fitted
% 
% 
% % met/sams to pc 
% % 'ahcys[c]'	'RC00192'	'RC02037'	'P/C'	0.465621903520335	1x1 table	'h2o[c] + ahcys[c]  -> hcys-L[c] + adn[c] '	'amet[c] + ethamp[c]  -> h[c] + ahcys[c] + methamp[c] '	15
% % 'ahcys[c]'	'RC00192'	'RCC0021'	'P/C'	0.507813358622038	1x1 table	'h2o[c] + ahcys[c]  -> hcys-L[c] + adn[c] '	'2 amet[c] + methamp[c]  -> 2 h[c] + 2 ahcys[c] + cholp[c] '	15
% metFit(strcmp(metFit.branchMets,'ahcys[c]'),:)
% a(ismember(a.rxn, {'RC00192' 'RC02037','RCC0021'}),:)
% % 'amet[c]'	'RC00177'	'RC02037'	'P/C'	0.569789300587654	4x1 table	'atp[c] + h2o[c] + met-L[c]  -> pi[c] + ppi[c] + amet[c] '	'amet[c] + ethamp[c]  -> h[c] + ahcys[c] + methamp[c] '	16
% % 'amet[c]'	'RC00177'	'RCC0021'	'P/C'	0.460663269582792	4x1 table	'atp[c] + h2o[c] + met-L[c]  -> pi[c] + ppi[c] + amet[c] '	'2 amet[c] + methamp[c]  -> 2 h[c] + 2 ahcys[c] + cholp[c] '	16
% metFit(strcmp(metFit.branchMets,'amet[c]'),:)
% a(ismember(a.rxn, {'RC00177' 'RC02037','RCC0021'}),:)
% % --> seems constrained so not fitted, but generally consistent (still sig
% % unfit)
% 
% 
% % BCAA deg 
% % 'ivcoa[m]'	'RMC0011'	'RM04096'	'P/C'	0.366269350862776	4x2 table	'coa[m] + nad[m] + 4mop[m]  -> co2[m] + nadh[m] + ivcoa[m] '	'ivcoa[m] + etfox[m]  -> 3mb2coa[m] + etfrd[m] '	4
% % 'ivcoa[m]'	'RM04096'	'RMC0105'	'P/C'	0.204113332978664	2x4 table	'ivcoa[m] + etfox[m]  -> 3mb2coa[m] + etfrd[m] '	'2 coa[m] + 2 nad[m] + 2 h2o[m] + 2 fad[m] + fa9p0isocoa[m]  -> 2 accoa[m] + 2 nadh[m] + 2 h[m] + 2 fadh2[m] + ivcoa[m] '	4
% % '3mb2coa[m]'	'RM04096'	'RM04138'	'P/C'	0.213892985033434	2x2 table	'ivcoa[m] + etfox[m]  -> 3mb2coa[m] + etfrd[m] '	'atp[m] + hco3[m] + 3mb2coa[m]  -> h[m] + adp[m] + pi[m] + 3mgcoa[m] '	2
% metFit(strcmp(metFit.branchMets,'ivcoa[m]'),:)
% a(ismember(a.rxn, {'RMC0011' 'RM04096','RMC0105'}),:)
% % --> largely fitted
% metFit(strcmp(metFit.branchMets,'3mb2coa[m]'),:)
% a(ismember(a.rxn, {'RM04138' 'RM04096'}),:)
% % --> fully fitted
% 
% 
% % tyr and phe deg 
% % '34hpp[c]'	'RC00734'	'RC02521'	'P/C'	0.437593126043506	2x1 table	'akg[c] + tyr-L[c]  <=> glu-L[c] + 34hpp[c] '	'o2[c] + 34hpp[c]  -> co2[c] + hgentis[c] '	5
% % 'phpyr[c]'	'RC00694'	'RC01372'	'P/C'	0.437593126043506	2x1 table	'akg[c] + phe-L[c]  <=> glu-L[c] + phpyr[c] '	'o2[c] + phpyr[c]  -> co2[c] + 2hyoxplac[c] '	5
% metFit(strcmp(metFit.branchMets,'34hpp[c]'),:)
% a(ismember(a.rxn, {'RC00734' 'RC02521'}),:)
% % fully fitted 
% metFit(strcmp(metFit.branchMets,'phpyr[c]'),:)
% a(ismember(a.rxn, {'RC00694' 'RC01372'}),:)
% % fitted by zero flux
% 
% 
% % mevolonate, udgp, gdpmannn to n-glycan 
% % 'ipdp[c]'	'RC01121'	'RC05556'	'P/C'	0.247223549792247	1x1 table	'atp[c] + 5dpmev[c]  -> adp[c] + pi[c] + co2[c] + ipdp[c] '	'19 ipdp[c] + frdp[c]  -> 19 ppi[c] + dedoldp[c] '	7
% % 'gdpmann[c]'	'RC05972'	'RC00885'	'P/C'	0.888561090465890	1x1 table	'chito2pdol[c] + gdpmann[c]  -> h[c] + gdp[c] + mpdol[c] '	'h[c] + gtp[c] + man1p[c]  -> ppi[c] + gdpmann[c] '	7
% % 'gdpmann[c]'	'RC05973'	'RC00885'	'P/C'	0.779410925866119	1x1 table	'gdpmann[c] + mpdol[c]  -> h[c] + gdp[c] + m1mpdol[c] '	'h[c] + gtp[c] + man1p[c]  -> ppi[c] + gdpmann[c] '	7
% % 'gdpmann[c]'	'RC06238'	'RC00885'	'P/C'	0.779410925866119	1x1 table	'gdpmann[c] + m1mpdol[c]  -> h[c] + gdp[c] + m2mpdol[c] '	'h[c] + gtp[c] + man1p[c]  -> ppi[c] + gdpmann[c] '	7
% % 'gdpmann[c]'	'RCC0037'	'RC00885'	'P/C'	0.420349242047196	1x1 table	'2 gdpmann[c] + m2mpdol[c]  -> 2 h[c] + 2 gdp[c] + m4mpdol[c] '	'h[c] + gtp[c] + man1p[c]  -> ppi[c] + gdpmann[c] '	7
% % 'gdpmann[c]'	'RC01009'	'RC00885'	'P/C'	0.586272423831442	1x1 table	'dolp[c] + gdpmann[c]  -> gdp[c] + dolmanp[c] '	'h[c] + gtp[c] + man1p[c]  -> ppi[c] + gdpmann[c] '	7
% % 'udpg[c]'	'RC01005'	'RC00289'	'P/C'	0.811145061649767	1x1 table	'dolp[c] + udpg[c]  -> udp[c] + dolglcp[c] '	'g1p[c] + h[c] + utp[c]  -> ppi[c] + udpg[c] '	10
% % 'man6p[c]'	'RC01326'	'RC01818'	'P/C'	0.234741099054269	1x1 table	'atp[c] + man[c]  -> adp[c] + h[c] + man6p[c] '	'man1p[c]  <=> man6p[c] '	3
% % 'man1p[c]'	'RC01818'	'RC00885'	'P/C'	0.634268369450233	1x1 table	'man1p[c]  <=> man6p[c] '	'h[c] + gtp[c] + man1p[c]  -> ppi[c] + gdpmann[c] '	2
% metFit(strcmp(metFit.branchMets,'ipdp[c]'),:)
% a(ismember(a.rxn, {'RC01121' 'RC05556'}),:)
% % --> only largely fitted, as expected bc of stoichemitry for entire
% % n-glycan syn
% metFit(strcmp(metFit.branchMets,'gdpmann[c]'),:)
% a(ismember(a.rxn, {'RC05973' 'RC05972','RC00885','RC06238','RCC0037','RC01009'}),:)
% % largely fitted (as so in original)
% metFit(strcmp(metFit.branchMets,'udpg[c]'),:)
% a(ismember(a.rxn, {'RC01005' 'RC00289'}),:)
% % [INTERESTING] not fitted and will be interesting to find out what is the
% % constainnt that prevents the fit 
% metFit(strcmp(metFit.branchMets,'man6p[c]'),:)
% a(ismember(a.rxn, {'RC01326' 'RC01818'}),:)
% % [INTERESTING] not fitted and will be interesting to find out what is the
% % constainnt that prevents the fit 
% metFit(strcmp(metFit.branchMets,'man1p[c]'),:)
% a(ismember(a.rxn, {'RC00885' 'RC01818'}),:)
% % --> fully fitted 
% 
% 
% % FA oxidation and the oxidation reagent 
% % 'ficytb[c]'	'RC00100'	'RC02222'	'P/C'	0.580774015810685	2x2 table	'h[c] + nadh[c] + 2 ficytb[c]  -> nad[c] + 2 focytb[c] '	'o2[c] + 2 focytb[c] + stcoa[c]  -> 2 h2o[c] + 2 ficytb[c] + odecoa[c] '	12
% % 'ficytb[c]'	'RC00100'	'RCC0061'	'P/C'	0.228898886594997	2x1 table	'h[c] + nadh[c] + 2 ficytb[c]  -> nad[c] + 2 focytb[c] '	'o2[c] + 2 focytb[c] + odecoa[c]  -> 2 h2o[c] + 2 ficytb[c] + lnlccoa[c] '	12
% % 'ficytb[c]'	'RC00100'	'RCC0069'	'P/C'	0.297508166554840	2x1 table	'h[c] + nadh[c] + 2 ficytb[c]  -> nad[c] + 2 focytb[c] '	'o2[c] + pmtcoa[c] + 2 focytb[c]  -> 2 h2o[c] + 2 ficytb[c] + fa16p1n7coa[c] '	12
% metFit(strcmp(metFit.branchMets,'ficytb[c]'),:)
% a(ismember(a.rxn, {'RC00100' 'RC02222','RCC0061','RCC0069'}),:)
% % [IMPORTANT] not largely fitted and will be interesting to find out what is the
% % constainnt that prevents the fit 
% 
% 
% % GSH syn 
% % 'glucys[c]'	'RC00894'	'RC00497'	'P/C'	0.261174647006175	2x1 table	'atp[c] + glu-L[c] + cys-L[c]  -> adp[c] + h[c] + pi[c] + glucys[c] '	'atp[c] + gly[c] + glucys[c]  -> adp[c] + h[c] + pi[c] + gthrd[c] '	4
% metFit(strcmp(metFit.branchMets,'glucys[c]'),:)
% a(ismember(a.rxn, {'RC00894' 'RC00497'}),:)
% % fully fitted although at epsilon
% 
% % unsaturated FA syn
% % 'fa16p1n7coa[c]'	'RCC0069'	'RCC0070'	'P/C'	0.493184449254954	1x13 table	'o2[c] + pmtcoa[c] + 2 focytb[c]  -> 2 h2o[c] + 2 ficytb[c] + fa16p1n7coa[c] '	'3 h[c] + 2 nadph[c] + malcoa[c] + fa16p1n7coa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + vacccoa[c] '	8
% % 'fa16p1n7coa[c]'	'RCC0070'	'RCC0089'	'P/C'	0.420804806043524	13x4 table	'3 h[c] + 2 nadph[c] + malcoa[c] + fa16p1n7coa[c]  -> h2o[c] + co2[c] + 2 nadp[c] + coa[c] + vacccoa[c] '	'h2o[c] + nad[c] + coa[c] + o2[c] + odecoa[c]  -> h[c] + nadh[c] + accoa[c] + h2o2[c] + fa16p1n7coa[c] '	8
% metFit(strcmp(metFit.branchMets,'fa16p1n7coa[c]'),:)
% a(ismember(a.rxn, {'RCC0069' 'RCC0070','RCC0089'}),:)
% % --> fully fitted
% 
% % many cases where degree = 2 is as expected because they are linear
% % pathway
% 
% %% check by metabolite
% met = 'accoa[m]'; % 'accoa[c]' leu-L
% tbl3 = listRxn(model,myCSM_merged_alt.OFD,met); % myCSM_merged.OFD
% tbl4 = listRxn(model,myCSM_exp.OFD,met); % myCSM_merged.OFD
% 
% %% more visualization and inspection
% fluxTbl = table(model.rxns, myCSM_merged_alt.OFD, printRxnFormula(model, model.rxns,0));
% fluxTbl.type = regexprep(fluxTbl.Var1,'[0-9]+','');
% sortedTbl = sortrows(fluxTbl, [4, 2]);
% 
% %% FVA of a single rxn
% targetRxns = {'RC01070'};
% parforFlag = 0;
% relMipGapTol = 1e-3; % this is to be redone with maximum precision, otherwsie the box will have problem
% verbose = false;
% 
% % run FVA by calling:
% [minval, maxval] = FVA_MILP(myCSM.MILP, model, targetRxns,parforFlag,relMipGapTol,verbose)
% 
% 
% %% make the inspection tabel 
% fluxTbl = table(model.rxns, myCSM_triple.OFD, printRxnFormula(model, model.rxns,0));
% fluxTbl.type = regexprep(fluxTbl.Var1,'[0-9]+','');
% fluxTbl.lb = minval(1:length(model.rxns))';
% fluxTbl.ub = maxval(1:length(model.rxns))';
% fluxTbl.tightness = ismember(model.rxns,tightRxns);
% fluxTbl = sortrows(fluxTbl, [7, 4, 2],"descend");
% 
% fluxTbl_triple = fluxTbl;

