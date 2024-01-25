%% Load model - this is the same across all five integrations 
function model = configurateModel(model)

% setting up the model 

% stay consistent with regular worm FBA 
model = changeRxnBounds(model,'RMC0005',0,'u');
model = changeRxnBounds(model,'RMC0005',0,'l');

% block a reaction that is a misannotation in the model and should be not
% existing
model = changeRxnBounds(model,'TCM1071',0,'b');

% constraints to optimize the modeling of redox balance
% the core conclusion such as cyclic PPP flux (bounded large flux) is not
% sensitive to the inclusion or exclusive of these redox constraints

% modify the nnt-1 reaction; assume producing nadph and assume not massive production (since H+ gradient coupling was not reconstructed)
model = changeRxnBounds(model,'RM00112',0,'l');
model = changeRxnBounds(model,'RM00112',0.01,'u');

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

% assuming reducing direction of GDH in vivo (PMID: 28208702)
model = changeRxnBounds(model,'RM00248',0,'l'); % assume glutamate dehydrogenase only in forward direction (making nh4) in our study
model = changeRxnBounds(model,'RM00243',0,'l'); % assume glutamate dehydrogenase only in forward direction (making nh4) in our study

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



% block the thermodynamically infeasible loops that can use the low-energy
% bound in ATP as a high-energy bound
% the massive PPI produced in tRNA synthesis can be converted to GTP in a
% loop that further support protein synthesis; this is an infeasible loop
% to bypass real energy demand; we prevent this loop
model = changeRxnBounds(model,'RCC0139',0,'l'); % although BRENDA supports reversible, a sig. reverse flux is not likely feasible and this is the setting in human model


% rescue the fitting of a few genes 
model = changeRxnBounds(model,'DMN0033',0.01,'u'); % to fit responsive gene F39B2.3, demand q2 has to be active; we assume minimal flux
model = changeRxnBounds(model,'EX00535',-0.01,'l'); % to fit responsive gene F42F12.3, demand tststerone has to be active; we assume minimal flux for hormone

