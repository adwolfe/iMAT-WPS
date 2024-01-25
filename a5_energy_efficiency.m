%% load model 
% Load model - this is the same across all five integrations 
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');

% setup the model
model = configurateModel(model);
model = changeRxnBounds(model,'EXC0050',-1,'l');% free bacteria uptake for integration

% to correctly model the energy balance given our model does not have inner
% membrance space (so we have proton leak issues), we did some fix here

% we couple the flux though oxidative PPP (they will anyway be coupled by
% inner membrane H+ in the model with inner space)
model.S(end+1,:) = zeros(1,length(model.rxns));
model.S(end, strcmp('RMC0006',model.rxns)) = 4; 
model.S(end, strcmp('RM02161',model.rxns)) = 4; 
model.S(end, strcmp('RMC0184',model.rxns)) = 4; 
model.S(end, strcmp('RM00081',model.rxns)) = 4; 
model.S(end, strcmp('RMC0004',model.rxns)) = -4; 
model.csense(end+1) = 'E';
model.b(end+1) = 0;
model.BiGG(end+1) = {'NA'};
model.metCharges(end+1) = 0;
model.metFormulas(end+1) = {'NA'};
model.metNames(end+1) = {'Pseudo-metabolite that represents a constraint coupling the ETC fluxes'};
model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};

% second, block the reverse flux through glyoxylate shunt (which will make
% a more energy-efficient shortcut TCA cycle that is not likely feasible in
% physiological conditions)
model = changeRxnBounds(model,'RM00479', 0, 'l');

model = changeRxnBounds(model,'RCC0005',0,'l');% remove the ATPm

%% stoichemitry-focused FBA analysis of theoretical ATP production
% calculate the maximum ATP production of candidate nutrients

ATPyield = [];
checkMets = {'glc-D[e]','rib-D[e]','ppa[e]','hdca[e]','rna_BAC[e]'};

for i = 1:length(checkMets)
    myMet = checkMets{i};% ppcoa[m]
    
    testing_model = changeRxnBounds(model,'EXC0050',0,'l');% block bacteria uptake
    testing_model = changeObjective(testing_model, 'RCC0005');
    
    % add a controlled influx
    testing_model = addReaction(testing_model,'mySink','reactionFormula',...
        [myMet, ' <=> '],'geneRule', 'NA','printLevel',1);
    testing_model = changeRxnBounds(testing_model, 'mySink',-1,'l');
    testing_model = changeRxnBounds(testing_model, 'mySink',0,'u');
    maxflux = optimizeCbModel(testing_model);
    % testing_model = changeRxnBounds(testing_model, 'RCC0005',maxflux.f,'l');
    % myflux = minimizeModelFlux_XL(testing_model);
    ATPyield(i,1) = maxflux.f;
    i
end

% ATP yield per gram dw or per mole for every metabolite and see if is
% among the top 
% even if it is not, BAC per gram ATP yield can be a reference point 

MW = [];
for i = 1:length(checkMets)
    x = model.metFormulas{strcmp(model.mets,checkMets{i})};
    if ~isempty(regexp(x, 'X|R|n|A','once'))
        MW(i,1) = nan;
    else 
        MW(i,1) = MolMass(x); 
    end
end

ATPyield_per_gram = ATPyield ./ MW;


checkMets
ATPyield % mole ATP yeild per mole
ATPyield_per_gram % mole ATP per gram 

% glc-D[e] 32 mole ATP/mole; 0.1776 mole ATP/gram
% rib-D[e] 26.6667 mole ATP/mole; 0.1776 mole ATP/gram
% ppa[e]  14.5 mole ATP/mole; 0.1984 mole ATP/gram
% hdca[e] 106 mole ATP/mole; 0.4150 mole ATP per gram
% rna_BAC[e] 34.9828 mole ATP/mole (not meaningful); 0.1081 mole ATP per gram

% glucose and ribose are equally well in producing ATP

%% maximize biomass and minimize total flux
% pFBA simulation to see the flux usage in optimal ATP production

obj = 'RCC0005';

testing_model = changeObjective(model, obj);
maxflux = optimizeCbModel(testing_model);
fluxMat = [];
maxFluxMat = [];

for i = 1:100
    testing_model = changeRxnBounds(testing_model, obj, maxflux.f .* i ./ 100,'l');
    myflux = minimizeModelFlux_XL(testing_model);
    fluxMat = [fluxMat, myflux];
%     for j = 1:length(targets)
%         tmp = changeObjective(testing_model, targets{j});
%         f_tmp = optimizeCbModel(tmp);
%         maxFluxMat(j,i) = f_tmp.f;
%     end
end

% plot if these reactions of interest carry flux
hold on 
plot(1:100, fluxMat(strcmp(testing_model.rxns,'RM04432'),:)) % shunt
% plot(1:100, fluxMat(strcmp(testing_model.rxns,'RM00705'),:)) % shunt
% plot(1:100, fluxMat(strcmp(testing_model.rxns,'RM00209'),:)) % PDH
plot(1:100, fluxMat(strcmp(testing_model.rxns,'RC01057'),:)) % r1p - r5p
% plot(1:100, fluxMat(strcmp(testing_model.rxns,'RC00754'),:)) % etoh
plot(1:100, sum(fluxMat(ismember(testing_model.rxns,{'RC03321','RC02740'}),:),1)) % g6p - f6p
% plot(1:100, fluxMat(strcmp(testing_model.rxns,'RC01641'),:)) % tkt-1 (r5p)
% plot(1:100, fluxMat(strcmp(testing_model.rxns,'RMC0005'),:)) % ATPase
hold off
legend({'RM04432:  etfox + ppcoa <=> etfrd + prpncoa'; ...
        % 'RM00705: coa + msa + nad -> accoa + co2 + nadh', ...
        'RC01057: r1p <=> r5p';...
        % 'RC00754: : etoh + nad <=> acald + nadh + h';...
        'RC03321 + RC02740: g6p <=> f6p'...
        % 'RC01641: g3p + s7p <=> xu5p + r5p'...
        },"Location","northwest")
yline(0,'--')
% 'RMC0005: ETC V' 'RM00209: PDH',
xlabel('ATP production (% maximum)')
ylabel('Flux in pFBA (mmole/gDW/h)')
hold off
saveas(gca, 'figures/energy_efficiency.pdf');

%% sensitivity analysis 
% first analzye all bulk nutrients
diet_cmp = {'DGR0001'; % protein
            'DGR0002'; % RNA
            'DGR0003'; % DNA
            'DGR0005'; % phospholipid
            'DGR0013'; % lps
            'DGR0028'; % glycogen
            'DGR0011'}; % peptido
% atp demand 
obj = 'RCC0005';

% run sensitivity analysis
objFluxMat = [];
for k = 1:length(diet_cmp)
    qry = diet_cmp{k};
    
    testing_model = changeObjective(model, qry);
    maxflux = optimizeCbModel(testing_model);
    objFlux = [];
    
    for i = 0:100
        testing_model = changeObjective(model, obj);
        testing_model = changeRxnBounds(testing_model,qry,  maxflux.f .* i ./ 100, 'b');
    %     testing_model.S(end+1,:) = zeros(1,length(testing_model.rxns));
    %     testing_model.S(end, ismember(testing_model.rxns, qry)) = 1; 
    %     testing_model.csense(end+1) = 'E';
    %     testing_model.b(end+1) =  maxflux.f .* i ./ 100;
    %     testing_model.BiGG(end+1) = {'NA'};
    %     testing_model.metCharges(end+1) = 0;
    %     testing_model.metFormulas(end+1) = {'NA'};
    %     testing_model.metNames(end+1) = {'Pseudo-metabolite that represents a constraint'};
    %     testing_model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};
        myflux = optimizeCbModel(testing_model);
        objFlux = [objFlux, myflux.f];
    end

    objFluxMat(k,:) = objFlux;
end

plot(0:100, objFluxMat)
legend(diet_cmp)

% pie chart 
testing_model = changeObjective(model, obj);
myflux = optimizeCbModel(testing_model);
maxATP = myflux.f;
deltaATP = maxATP - objFluxMat(:,1);
maxATPcontribution = [deltaATP; maxATP - sum(deltaATP)];
nutrients = {'protein'; % protein
            'RNA'; % RNA
            'DNA'; % DNA
            'phospholipid'; % phospholipid
            'lps'; % lps
            'glycogen'; % glycogen
            'peptido'; % peptido
            'others'}; 
figure;
pie(maxATPcontribution, nutrients);
title('Contribution to theoretical maximal ATP production');
saveas(gca, 'figures/max_energy_contribution_all_bulk_nutrients.pdf');

% run special analysis of amino acids
all_aa = model.mets(model.S(:,strcmp(model.rxns,'DGR0001')) > 0);

obj = 'RCC0005';

objFlux = [];
for k = 1:length(all_aa)
    qry = all_aa{k};
    % change obj to energy
    testing_model = changeObjective(model, obj);
    % block one aa intake
    testing_model.S(strcmp(model.mets,qry),strcmp(model.rxns,'DGR0001')) = 0;
    myflux = optimizeCbModel(testing_model);
    objFlux = [objFlux, myflux.f];
end

% pie chart 
testing_model = changeObjective(model, obj);
myflux = optimizeCbModel(testing_model);
maxATP = myflux.f;
deltaATP = maxATP - objFlux;
maxATPcontribution_AA = deltaATP; % [deltaATP, maxATPcontribution(1) - sum(deltaATP)];
nutrients_AA = [all_aa];
%            {'others'}];

figure;
pie(maxATPcontribution_AA, nutrients_AA);
title('Contribution to theoretical maximal ATP production');
saveas(gca, 'figures/max_energy_contribution_AA.pdf');

%% make plot that includes r1p
obj = 'RCC0005';

% first block dietary ribose from RNA
% change obj to energy
testing_model = changeObjective(model, obj);
% change amp to ade directly
testing_model.S(strcmp(testing_model.mets,'ade[c]'), ...
                strcmp(testing_model.rxns,'DGR0002') ...
    ) = testing_model.S(strcmp(testing_model.mets,'amp[c]'), ...
                strcmp(testing_model.rxns,'DGR0002') ...
    );
testing_model.S(strcmp(testing_model.mets,'amp[c]'), ...
                strcmp(testing_model.rxns,'DGR0002') ...
    ) = 0;
% change gmp to gua 
testing_model.S(strcmp(testing_model.mets,'gua[c]'), ...
                strcmp(testing_model.rxns,'DGR0002') ...
    ) = testing_model.S(strcmp(testing_model.mets,'gmp[c]'), ...
                strcmp(testing_model.rxns,'DGR0002') ...
    );
testing_model.S(strcmp(testing_model.mets,'gmp[c]'), ...
                strcmp(testing_model.rxns,'DGR0002') ...
    ) = 0;
% change cmp and ump to ura
testing_model.S(strcmp(testing_model.mets,'ura[c]'), ...
                strcmp(testing_model.rxns,'DGR0002') ...
    ) = sum(testing_model.S(ismember(testing_model.mets,{'cmp[c]','ump[c]'}), ...
                strcmp(testing_model.rxns,'DGR0002')) ...
    );
testing_model.S(ismember(testing_model.mets,{'cmp[c]','ump[c]'}), ...
                strcmp(testing_model.rxns,'DGR0002') ...
    ) = 0;
printRxnFormula(testing_model,'DGR0002');

myflux = optimizeCbModel(testing_model);
deltaATP_r1p = maxATP - myflux.f;

mixContribution = [deltaATP_r1p;
                   maxATPcontribution_AA(strcmp(nutrients_AA, 'phe-L[c]'));
                   maxATPcontribution_AA(strcmp(nutrients_AA, 'tyr-L[c]'));
                   maxATPcontribution_AA(strcmp(nutrients_AA, 'lys-L[c]'));
                   maxATPcontribution_AA(strcmp(nutrients_AA, 'leu-L[c]'));
                   maxATPcontribution_AA(strcmp(nutrients_AA, 'ile-L[c]'));
                   maxATPcontribution_AA(strcmp(nutrients_AA, 'val-L[c]'));
                   maxATPcontribution_AA(strcmp(nutrients_AA, 'met-L[c]'));
                   maxATPcontribution_AA(strcmp(nutrients_AA, 'thr-L[c]'));
                    ];
mixContribution = [mixContribution; maxATP - sum(mixContribution)];
mixNutrients = {'r1p';
                'phe-L';
                'tyr-L';
                'lys-L';
                'leu-L';
                'ile-L';
                'val-L';
                'met-L';
                'thr-L';
                'other'};

figure;
pie(mixContribution, mixNutrients);
title('Contribution to theoretical maximal ATP production');
saveas(gca, 'figures/max_energy_contribution_measured_nutrients.pdf');
%% special analysis of pyruvate and propinate - pyruvate is too hard to block
% diet_cmp = {{'RM04432','RM01859'};
%             {'RM00209'}};
% 
% obj = 'RCC0005';
% 
% objFluxMat = [];
% for k = 1:length(diet_cmp)
%     qry = diet_cmp{k};
%     
%     testing_model = changeObjective(model, qry);
%     maxflux = optimizeCbModel(testing_model);
%     objFlux = [];
%     
%     for i = 0:100
%         testing_model = changeObjective(model, obj);
%         testing_model.S(end+1,:) = zeros(1,length(testing_model.rxns));
%         testing_model.S(end, ismember(testing_model.rxns, qry)) = 1; 
%         testing_model.csense(end+1) = 'E';
%         testing_model.b(end+1) =  maxflux.f .* i ./ 100;
%         testing_model.BiGG(end+1) = {'NA'};
%         testing_model.metCharges(end+1) = 0;
%         testing_model.metFormulas(end+1) = {'NA'};
%         testing_model.metNames(end+1) = {'Pseudo-metabolite that represents a constraint'};
%         testing_model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};
%         myflux = optimizeCbModel(testing_model);
%         objFlux = [objFlux, myflux.f];
%     end
% 
%     objFluxMat(k,:) = objFlux;
% end
% 
% plot(0:100, objFluxMat)
% legend({'propionate degradation'; % peptido
%             'pyruvate dehydrogenase'})
% 
% % pie chart 
% testing_model = changeObjective(model, obj);
% myflux = optimizeCbModel(testing_model);
% maxATP = myflux.f;
% deltaATP = maxATP - objFluxMat(:,1);
% maxATPcontribution = [deltaATP; maxATP - sum(deltaATP)];
% nutrients = {'propionate degradation'; % peptido
%             'pyruvate dehydrogenase'}; 
% 
% bar(maxATPcontribution, nutrients);
% title('Theoretical ATP Contribution by Nutrients');
% 
