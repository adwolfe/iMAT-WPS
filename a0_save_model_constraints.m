%% save the fully constrianed model that is ready for CR simulation, with propoer notes attached 
initCobraToolbox(false);

%% Load model
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat'); % model 
model.rxnNotes = repmat({'Original iCEL1314 reaction and default constraints.'},length(model.rxns),1);

% setup the model
% stay consistent with regular worm FBA 
model = changeRxnBounds(model,'RMC0005',0,'u');
model = changeRxnBounds(model,'RMC0005',0,'l');
model.rxnNotes(strcmp(model.rxns,'RMC0005')) = {'Block one of the two alternative ATPase reactions, to simplify the solution space.'}

% block a reaction that is a misannotation in the model and should be not
% existing
model = changeRxnBounds(model,'TCM1071',0,'b'); 
model.rxnNotes(strcmp(model.rxns,'TCM1071')) = {'Block TCM1071 that is a misannotation in the model and should be not existing.'};


% constraints to optimize the modeling of redox balance
% the core conclusion such as cyclic PPP flux (bounded large flux) is not
% sensitive to the inclusion or exclusive of these redox constraints

% modify the nnt-1 reaction; assume producing nadph and assume not massive production (since H+ gradient coupling was not reconstructed)
model = changeRxnBounds(model,'RM00112',0,'l');
model = changeRxnBounds(model,'RM00112',0.01,'u'); % nnt-1 reaction; assume producing nadph and assume not massive production (since H+ gradient coupling was not reconstructed)
model.rxnNotes(strcmp(model.rxns,'RM00112')) = {'Constrain the flux through nnt-1 reaction to maintain proper redox balance. Since H+ gradient coupling was not reconstructed in the model, this reaction can convert nadh to nadph without any cost. Therefore, we constain it to minimal level of nadph production to avoid massive thermodynamically infeasible nadph production.'}


% allow the following canonical nadph producing pathways unconstrained:PMID: 24805240
% RM00248: glu --> akg
% RC02736, RC01528: PPP
% RC00267, RM00267: idh-1, icit --> akg
% RM01220, RC01220: folate (MTHFD)
% RM00216: malate --> pyruvate

% To avoid loops, we assume directions for canonical pathways if there is
% a nadh version of them

% assume both the directions of MTHFD and MTHFD2 as reducing direction
% (making nadh or nadph) (PMID: 24805240)
model = changeRxnBounds(model,'RC01218',0,'l'); 
model = changeRxnBounds(model,'RC01220',0,'l');
model.rxnNotes(ismember(model.rxns,{'RC01218','RC01220'})) = {'Assuming both the directions of MTHFD and MTHFD2 as reducing direction (producing NADPH/NADH) in vivo, based on PMID: 24805240'};

% assuming reducing direction of GDH in vivo (PMID: 28208702)
model = changeRxnBounds(model,'RM00248',0,'l'); % assume glutamate dehydrogenase only in forward direction (making nh4) in our study
model = changeRxnBounds(model,'RM00243',0,'l'); % assume glutamate dehydrogenase only in forward direction (making nh4) in our study
model.rxnNotes(ismember(model.rxns,{'RM00248','RM00243'})) = {'Assuming the direction of GDH as reducing direction (producing NADPH/NADH) in vivo, based on PMID: 28208702'};


% the following were constrianed to small flux (all other nadph
% producing reactions)
% of note, in most cases, the capacity of converting major reactant is
% not compromised by the constraint because there is another nadh
% version to support the flux
% RM00112 nnt-1 reaction was already constrained in above
model = changeRxnBounds(model,'RC04360',-0.01,'l');
model = changeRxnBounds(model,'RC04571',-0.01,'l');
model = changeRxnBounds(model,'RC07759',-0.01,'l');
model = changeRxnBounds(model,'RC01787',-0.01,'l');
model = changeRxnBounds(model,'RC05692',-0.01,'l');
model = changeRxnBounds(model,'RC01095',-0.01,'l');
model = changeRxnBounds(model,'RC02577',0.01,'u');
model = changeRxnBounds(model,'RC01041',-0.01,'l');
model = changeRxnBounds(model,'RM02566',0.01,'u');
model = changeRxnBounds(model,'RC00711',0.01,'u');
model = changeRxnBounds(model,'RM00711',0.01,'u');
model = changeRxnBounds(model,'RM00716',-0.01,'l');
model = changeRxnBounds(model,'RM03103',0.01,'u');
model = changeRxnBounds(model,'RC05623',-0.01,'l');
model = changeRxnBounds(model,'RC07140',0.01,'u');
model = changeRxnBounds(model,'RC08539',-0.01,'l');
model = changeRxnBounds(model,'RC00939',-0.01,'l');
model = changeRxnBounds(model,'RC02236',-0.01,'l');
model = changeRxnBounds(model,'RC01224',-0.01,'l');
model = changeRxnBounds(model,'RC00941',0.01,'u');
model = changeRxnBounds(model,'RC01904',-0.01,'l'); 
model = changeRxnBounds(model,'RC01431',-0.01,'l'); 
model = changeRxnBounds(model,'RC01759',-0.01,'l');
model = changeRxnBounds(model,'RC01481',-0.01,'l');
model = changeRxnBounds(model,'RM08759',0.01,'u'); 
model = changeRxnBounds(model,'RC08759',0.01,'u');
model = changeRxnBounds(model,'RM00706',0.01,'u');
model = changeRxnBounds(model,'RC00978',0.01,'u');
model = changeRxnBounds(model,'RC01415',0.01,'u');
model = changeRxnBounds(model,'RC08379',0.01,'u');
model = changeRxnBounds(model,'RC08383',0.01,'u');
model = changeRxnBounds(model,'RC03596',0.01,'u');
model = changeRxnBounds(model,'RC04940',-0.01,'l');
model = changeRxnBounds(model,'RC02082',-0.01,'l');
model = changeRxnBounds(model,'RC02697',0.01,'u');
model = changeRxnBounds(model,'RC03302',0.01,'u');
model = changeRxnBounds(model,'RC10059',0.01,'u');
model = changeRxnBounds(model,'RM03293',0.01,'u');
model.rxnNotes(ismember(model.rxns,{'RC04360','RC04571','RC07759','RC01787','RC05692','RC01095','RC02577',...
                                    'RC01041','RM02566','RC00711','RM00711','RM00716','RM03103','RC05623',...
                                    'RC07140','RC08539','RC00939','RC02236','RC01224','RC00941','RC01904',...
                                    'RC01431','RC01759','RC01481','RM08759','RC08759','RM00706','RC00978',...
                                    'RC01415','RC08379','RC08383','RC03596','RC04940','RC02082','RC02697',...
                                    'RC03302','RC10059','RM03293'})) = ...
                                    {'Restricting non-canonical NADPH production flux. Canonical reactions were selected based on PMID: 24805240 (incl. RM00248, RC02736, RC01528, RC00267, RM00267, RM01220, RC01220 and RM00216).'};


% block the thermodynamically infeasible loops that can use the low-energy
% bound in ATP as a high-energy bound
% the massive PPI produced in tRNA synthesis can be converted to GTP in a
% loop that further support protein synthesis; this is an infeasible loop
% to bypass real energy demand; we prevent this loop
model = changeRxnBounds(model,'RCC0139',0,'l'); % although BRENDA supports reversible, a sig. reverse flux is not likely feasible and this is the setting in human model
model.rxnNotes(strcmp(model.rxns,'RCC0139')) = {'Reaction RCC0139 was constrained to only allow forward flux (bounds = [0,1000]). This is because we found the reverse flux of this reaction can convert diphosphate (ppi) back to GTP without energy cost, forming a thermodynamically infeasible loop that recycles ppi produced in AA-tRNA synthesis to fuel GTP for protein synthesis without using energy (ATP). The modification of reversibility is supported by the reversibility annotation of the corresponding EC family (2.7.7.13) in the human model Human1.'};


% rescue the fitting of a few genes 
model = changeRxnBounds(model,'DMN0033',0.01,'u'); % to fit responsive gene F39B2.3, demand q2 has to be active; we assume minimal flux
model.rxnNotes(strcmp(model.rxns,'DMN0033')) = {'Only allow minimal cytosolic q regeneration since this reaction is not masss-balanced and causes redox leak. It represents oxidative stress per model reconstruction note #31.'}
model = changeRxnBounds(model,'EX00535',-0.01,'l'); % to fit responsive gene F42F12.3, demand tststerone has to be active; we assume minimal flux for hormone
model.rxnNotes(strcmp(model.rxns,'EX00535')) = {'Allowing minimal flux in an exchange reaction EX00535 because its flux is required for internal reactions to carry flux but this reaction was not included in the uptake reactions.'}


model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
model.rxnNotes(strcmp(model.rxns,'EXC0050')) = {'free bacteria uptake for iMAT-WPS integration.'}


%% save
model_out = struct();
model_out.rxnID  = model.rxns;
model_out.rxnFormula  = printRxnFormula(model,model.rxns,false);
model_out.LB  = model.lb;
model_out.UB  = model.ub;
model_out.note  = model.rxnNotes;
model_out = struct2table(model_out);

writetable(model_out,'output/iCEL1314_uptake_constriants.csv');
