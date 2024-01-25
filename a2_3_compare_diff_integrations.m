%% set parameters
sigFlux = 1e-5; 

%% load the results 
load('output/integration_output/myCSM_no_data.mat')
load('output/integration_output/myCSM_expression_only.mat')
load('output/integration_output/myCSM_exp_resp.mat')
load('output/integration_output/myCSM_exp_simi.mat')
load('output/integration_output/myCSM_exp_resp_simi.mat')

FVA = struct();
load('output/integration_output/FVA_PFD_exp_resp_simi.mat')
FVA.PFD.exp_resp_simi = [minval', maxval'];
load('output/integration_output/FVA_PFD_exp_resp.mat')
FVA.PFD.exp_resp = [minval', maxval'];
load('output/integration_output/FVA_PFD_exp_simi.mat')
FVA.PFD.exp_simi = [minval', maxval'];
load('output/integration_output/FVA_PFD_expression_only.mat')
FVA.PFD.exp_only = [minval', maxval'];
load('output/integration_output/FVA_PFD_no_data.mat')
FVA.PFD.no_data = [minval', maxval'];
load('output/integration_output/FVA_OFD_exp_resp_simi.mat')
FVA.OFD.exp_resp_simi = [minval', maxval'];
load('output/integration_output/FVA_OFD_exp_resp.mat')
FVA.OFD.exp_resp = [minval', maxval'];
load('output/integration_output/FVA_OFD_exp_simi.mat')
FVA.OFD.exp_simi = [minval', maxval'];
load('output/integration_output/FVA_OFD_expression_only.mat')
FVA.OFD.exp_only = [minval', maxval'];
load('output/integration_output/FVA_OFD_no_data.mat')
FVA.OFD.no_data = [minval', maxval'];
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');

%% combine the results into a table 
% make the heatmap of flux, FVA, and narrowing/selection of solution space 
% if a reaction runs in both direction in FVA (sig flux), it is unbounded 
% if a reaction is no flux or is unidirectional, it is bounded, so we
% calcualate the reduction of solution space by the corresponding boundary
% finally, check the number of unbounded reactions and reduction of
% solution space as a measure.

% 8/7/23 reduction of solution space is not systematically analyzed, we
% skip this metric

flux_cmp = table();
flux_cmp.rxns = model.rxns;
flux_cmp.formula = printRxnFormula(model, model.rxns,0);
flux_cmp.Properties.RowNames = flux_cmp.rxns;

% normalized flux (flux of per unit bacterial uptake)
flux_cmp.normalized_OFD_exp_resp_simi = myCSM_exp_resp_simi.OFD ./ abs(myCSM_exp_resp_simi.OFD(strcmp(model.rxns,'EXC0050'))); 
flux_cmp.normalized_OFD_exp_resp = myCSM_exp_resp.OFD ./ abs(myCSM_exp_resp.OFD(strcmp(model.rxns,'EXC0050')));
flux_cmp.normalized_OFD_exp_simi = myCSM_exp_simi.OFD ./ abs(myCSM_exp_simi.OFD(strcmp(model.rxns,'EXC0050')));
flux_cmp.normalized_OFD_exp_only = myCSM_exp_only.OFD ./ abs(myCSM_exp_only.OFD(strcmp(model.rxns,'EXC0050')));
flux_cmp.normalized_OFD_no_data = myCSM_no_data.OFD ./ abs(myCSM_no_data.OFD(strcmp(model.rxns,'EXC0050')));

% raw OFD flux
flux_cmp.OFD_exp_resp_simi = myCSM_exp_resp_simi.OFD; 
flux_cmp.OFD_exp_resp = myCSM_exp_resp.OFD;
flux_cmp.OFD_exp_simi = myCSM_exp_simi.OFD;
flux_cmp.OFD_exp_only = myCSM_exp_only.OFD;
flux_cmp.OFD_no_data = myCSM_no_data.OFD;

% add in the FVA interval 
% we cannot normalize FVA interval to per unit bacteria uptake because that
% needs a special FVA objective (i.e., max and min of Vi/Vbacteria_uptake)
flux_cmp.FVA_PFD_exp_resp_simi = FVA.PFD.exp_resp_simi; 
flux_cmp.FVA_PFD_exp_resp = FVA.PFD.exp_resp; 
flux_cmp.FVA_PFD_exp_simi = FVA.PFD.exp_simi; 
flux_cmp.FVA_PFD_exp_only = FVA.PFD.exp_only;
flux_cmp.FVA_PFD_no_data = FVA.PFD.no_data;
% annotate the FVA boundedness 
flux_cmp.PFD_bounded_exp_resp_simi  =  (flux_cmp.FVA_PFD_exp_resp_simi(:,2) < -sigFlux |... % negtive flux
                                        flux_cmp.FVA_PFD_exp_resp_simi(:,1) > sigFlux |... % positive flux 
                                        (flux_cmp.FVA_PFD_exp_resp_simi(:,1) > -sigFlux & flux_cmp.FVA_PFD_exp_resp_simi(:,2) < sigFlux)); % zero flux 
flux_cmp.PFD_bounded_exp_resp  =  (flux_cmp.FVA_PFD_exp_resp(:,2) < -sigFlux |... % negtive flux
                                    flux_cmp.FVA_PFD_exp_resp(:,1) > sigFlux |... % positive flux 
                                    (flux_cmp.FVA_PFD_exp_resp(:,1) > -sigFlux & flux_cmp.FVA_PFD_exp_resp(:,2) < sigFlux)); % zero flux 
flux_cmp.PFD_bounded_exp_simi  =  (flux_cmp.FVA_PFD_exp_simi(:,2) < -sigFlux |... % negtive flux
                                    flux_cmp.FVA_PFD_exp_simi(:,1) > sigFlux |... % positive flux 
                                    (flux_cmp.FVA_PFD_exp_simi(:,1) > -sigFlux & flux_cmp.FVA_PFD_exp_simi(:,2) < sigFlux)); % zero flux 
flux_cmp.PFD_bounded_exp_only  =  (flux_cmp.FVA_PFD_exp_only(:,2) < -sigFlux |... % negtive flux
                                        flux_cmp.FVA_PFD_exp_only(:,1) > sigFlux |... % positive flux 
                                        (flux_cmp.FVA_PFD_exp_only(:,1) > -sigFlux & flux_cmp.FVA_PFD_exp_only(:,2) < sigFlux)); % zero flux 
flux_cmp.PFD_bounded_no_data  =  (flux_cmp.FVA_PFD_no_data(:,2) < -sigFlux |... % negtive flux
                                        flux_cmp.FVA_PFD_no_data(:,1) > sigFlux |... % positive flux 
                                        (flux_cmp.FVA_PFD_no_data(:,1) > -sigFlux & flux_cmp.FVA_PFD_no_data(:,2) < sigFlux)); % zero flux 

% % calculate the tightness of solution space - skipped
% flux_cmp.PFD_tightness_triple = ones(size(flux_cmp,1),1); % 1 = free FVA
% flux_cmp.PFD_tightness_merged = ones(size(flux_cmp,1),1); % 1 = free FVA
% flux_cmp.PFD_tightness_merged_alt = ones(size(flux_cmp,1),1); % 1 = free FVA
% flux_cmp.PFD_tightness_expression = ones(size(flux_cmp,1),1); % 1 = free FVA
% flux_cmp.PFD_tightness_base = ones(size(flux_cmp,1),1); % 1 = free FVA

% same for OFD FVA
flux_cmp.FVA_OFD_exp_resp_simi = FVA.OFD.exp_resp_simi; 
flux_cmp.FVA_OFD_exp_resp = FVA.OFD.exp_resp; 
flux_cmp.FVA_OFD_exp_simi = FVA.OFD.exp_simi; 
flux_cmp.FVA_OFD_exp_only = FVA.OFD.exp_only;
flux_cmp.FVA_OFD_no_data = FVA.OFD.no_data;

flux_cmp.OFD_bounded_exp_resp_simi  =  (flux_cmp.FVA_OFD_exp_resp_simi(:,2) < -sigFlux |... % negtive flux
                                        flux_cmp.FVA_OFD_exp_resp_simi(:,1) > sigFlux |... % positive flux 
                                        (flux_cmp.FVA_OFD_exp_resp_simi(:,1) > -sigFlux & flux_cmp.FVA_OFD_exp_resp_simi(:,2) < sigFlux)); % zero flux 
flux_cmp.OFD_bounded_exp_resp  =  (flux_cmp.FVA_OFD_exp_resp(:,2) < -sigFlux |... % negtive flux
                                    flux_cmp.FVA_OFD_exp_resp(:,1) > sigFlux |... % positive flux 
                                    (flux_cmp.FVA_OFD_exp_resp(:,1) > -sigFlux & flux_cmp.FVA_OFD_exp_resp(:,2) < sigFlux)); % zero flux 
flux_cmp.OFD_bounded_exp_simi  =  (flux_cmp.FVA_OFD_exp_simi(:,2) < -sigFlux |... % negtive flux
                                    flux_cmp.FVA_OFD_exp_simi(:,1) > sigFlux |... % positive flux 
                                    (flux_cmp.FVA_OFD_exp_simi(:,1) > -sigFlux & flux_cmp.FVA_OFD_exp_simi(:,2) < sigFlux)); % zero flux 
flux_cmp.OFD_bounded_exp_only  =  (flux_cmp.FVA_OFD_exp_only(:,2) < -sigFlux |... % negtive flux
                                        flux_cmp.FVA_OFD_exp_only(:,1) > sigFlux |... % positive flux 
                                        (flux_cmp.FVA_OFD_exp_only(:,1) > -sigFlux & flux_cmp.FVA_OFD_exp_only(:,2) < sigFlux)); % zero flux 
flux_cmp.OFD_bounded_no_data  =  (flux_cmp.FVA_OFD_no_data(:,2) < -sigFlux |... % negtive flux
                                        flux_cmp.FVA_OFD_no_data(:,1) > sigFlux |... % positive flux 
                                        (flux_cmp.FVA_OFD_no_data(:,1) > -sigFlux & flux_cmp.FVA_OFD_no_data(:,2) < sigFlux)); % zero flux 

% % calculate the tightness of solution space - skipped
% flux_cmp.OFD_tightness_triple = ones(size(flux_cmp,1),1); % 1 = free FVA
% flux_cmp.OFD_tightness_merged = ones(size(flux_cmp,1),1); % 1 = free FVA
% flux_cmp.OFD_tightness_merged_alt = ones(size(flux_cmp,1),1); % 1 = free FVA
% flux_cmp.OFD_tightness_expression = ones(size(flux_cmp,1),1); % 1 = free FVA
% flux_cmp.OFD_tightness_base = ones(size(flux_cmp,1),1); % 1 = free FVA
% for i = 1:size(flux_cmp,1)
%     for j = [13:17]
%         if flux_cmp{i,j} % bounded 
%             if flux_cmp{i,j-5}(1) > sigFlux % positive flux 
%                flux_cmp{i,j+5} = (flux_cmp{i,j-5}(2) - flux_cmp{i,j-5}(1)) / ... % solution space
%                                     (flux_cmp.base_FVA_PFD(i,2) - max(flux_cmp.base_FVA_PFD(i,1),0)); % data-independent solution space (always use OFD since PFD is too loose)
%             elseif flux_cmp{i,j-5}(2) < -sigFlux % negative flux 
%                flux_cmp{i,j+5} = (flux_cmp{i,j-5}(2) - flux_cmp{i,j-5}(1)) / ... % solution space
%                                     (min(flux_cmp.base_FVA_PFD(i,2),0) - flux_cmp.base_FVA_PFD(i,1)); % data-independent solution space
%             elseif flux_cmp{i,j-5}(2) < sigFlux && flux_cmp{i,j-5}(1) > -sigFlux % zero flux
%                 flux_cmp{i,j+5} = 0;
%             else
%                 error('?')
%             end
%         end
%     end
% 
%     for j = [28:32]
%         if flux_cmp{i,j} % bounded 
%             if flux_cmp{i,j-5}(1) > sigFlux % positive flux 
%                flux_cmp{i,j+5} = (flux_cmp{i,j-5}(2) - flux_cmp{i,j-5}(1)) / ... % solution space
%                                     (flux_cmp.base_FVA_OFD(i,2) - max(flux_cmp.base_FVA_OFD(i,1),0)); % data-independent solution space (always use OFD since PFD is too loose)
%             elseif flux_cmp{i,j-5}(2) < -sigFlux % negative flux 
%                flux_cmp{i,j+5} = (flux_cmp{i,j-5}(2) - flux_cmp{i,j-5}(1)) / ... % solution space
%                                     (min(flux_cmp.base_FVA_OFD(i,2),0) - flux_cmp.base_FVA_OFD(i,1)); % data-independent solution space
%             elseif flux_cmp{i,j-5}(2) < sigFlux && flux_cmp{i,j-5}(1) > -sigFlux % zero flux
%                 flux_cmp{i,j+5} = 0;
%             else
%                 error('?')
%             end
%         end
%     end
% end

% flux_cmp('BIO0010',:)
% flux_cmp('EXC0050',:)
writetable(flux_cmp,'output/fluxTable.csv')

%% compare the number of bounded reactions

% to show the relationship of bounded rxns in different integrations, we
% make an euler plot in R 
% bar plot here is only to display numbers
%-- there are too many categories - bar plot is not working - go for ruler


N_bounded = [sum(flux_cmp.OFD_bounded_exp_resp_simi);
            sum(flux_cmp.OFD_bounded_exp_resp);
            sum(flux_cmp.OFD_bounded_exp_simi);
            sum(flux_cmp.OFD_bounded_exp_only);
            sum(flux_cmp.OFD_bounded_no_data);];

figure;
bar(N_bounded)
hold on 
% plot in R
bar(repmat(sum(flux_cmp.OFD_bounded_exp_resp_simi & flux_cmp.OFD_bounded_exp_simi & flux_cmp.OFD_bounded_exp_resp & flux_cmp.OFD_bounded_exp_only),4,1));
bar(repmat(sum(flux_cmp.OFD_bounded_exp_resp_simi & flux_cmp.OFD_bounded_exp_simi & flux_cmp.OFD_bounded_exp_resp & flux_cmp.OFD_bounded_exp_only & flux_cmp.OFD_bounded_no_data),5,1));
hold off
legend({'others', ...
        'exp+resp+simi ∩  exp+resp. ∩ exp+simi. ∩ exp' ...
        ,'exp+resp+simi ∩ exp+resp. ∩ exp+simi. ∩ exp ∩ no-data'})
ylabel('# rxns')
set(gca,'xticklabel',{'exp+resp+simi','exp+resp.','exp+simi.', 'exp','no-data'})
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85, 2.85];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.LegendLoc = 'SouthOutside';
plt.export(['figures/Number_of_bounded_OFD_reactions.pdf']);

sum(flux_cmp.OFD_bounded_exp_resp_simi | flux_cmp.OFD_bounded_exp_simi | flux_cmp.OFD_bounded_exp_resp | flux_cmp.OFD_bounded_exp_only | flux_cmp.OFD_bounded_no_data)

% PFD
N_bounded = [sum(flux_cmp.PFD_bounded_exp_resp_simi);
            sum(flux_cmp.PFD_bounded_exp_resp);
            sum(flux_cmp.PFD_bounded_exp_simi);
            sum(flux_cmp.PFD_bounded_exp_only);
            sum(flux_cmp.PFD_bounded_no_data);];


figure;
bar(N_bounded)
hold on 
bar(repmat(sum(flux_cmp.PFD_bounded_exp_resp_simi & flux_cmp.PFD_bounded_exp_resp & flux_cmp.PFD_bounded_exp_simi & flux_cmp.PFD_bounded_exp_only),4,1));
bar(repmat(sum(flux_cmp.PFD_bounded_exp_resp_simi & flux_cmp.PFD_bounded_exp_resp & flux_cmp.PFD_bounded_exp_simi & flux_cmp.PFD_bounded_exp_only & flux_cmp.PFD_bounded_no_data),5,1));
hold off
legend({'others', ...
        'exp+resp+simi ∩  exp+resp. ∩ exp+simi. ∩ exp' ...
        ,'exp+resp+simi ∩ exp+resp. ∩ exp+simi. ∩ exp ∩ no-data'})
ylabel('# rxns')
set(gca,'xticklabel',{'exp+resp+simi','exp+resp.','exp+simi.', 'exp','no-data'})
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85, 2.85];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.LegendLoc = 'SouthOutside';
plt.export(['figures/Number_of_bounded_PFD_reactions.pdf']);


% tightness of solution space analysis has been discountinued
% compare the tightness of solution space
% we need a better metric. normalizing to base doesnt seem very
% interpretable
% 
% figure;
% hold on
% h1 = histogram(flux_cmp.tightness_expression);
% h1.FaceColor = '#808080';
% h2 = histogram(flux_cmp.tightness_merged, 'BinEdges',h1.BinEdges);
% h2.FaceColor = '#FFA500';
% h3 = histogram(flux_cmp.tightness_triple, 'BinEdges',h1.BinEdges);
% h3.FaceColor = '#ff5900';
% xlabel('Tightness (space w/ constraints / space w/o const.')
% ylabel('# rxns')
% legend({'expression', 'responsiveness + expression', 'similarity + responsiveness + expression'})
% hold off
% plt = Plot(); % create a Plot object and grab the current figure
% plt.BoxDim = [2.85, 2.35];
% plt.LineWidth = 1;
% plt.FontSize = 7;
% plt.FontName = 'Arial';
% plt.ShowBox = 'off';
% plt.XMinorTick = 'off';
% plt.YMinorTick = 'off';
% plt.TickDir = 'out';
% % plt.LegendLoc = 'SouthEast';
% plt.export(['figures/tightness_of_solution_space.pdf']);
% 
% 
% figure;
% plot(flux_cmp.tightness_expression, flux_cmp.tightness_merged,'.')
% hold on 
% plot(flux_cmp.tightness_expression, flux_cmp.tightness_triple,'.r')
% xlabel('Absolute expression only')
% ylabel('Data integration')
% legend({'responsiveness','responsiveness + similarity'})
% xlim([-0.05,1.05])
% ylim([-0.05,1.05])
% l = refline(1,0);
% l.Color = 'red';
% hold off
% plt = Plot(); % create a Plot object and grab the current figure
% plt.BoxDim = [2.85, 2.35];
% plt.LineWidth = 1;
% plt.FontSize = 7;
% plt.FontName = 'Arial';
% plt.ShowBox = 'off';
% plt.XMinorTick = 'off';
% plt.YMinorTick = 'off';
% plt.TickDir = 'out';
% plt.export(['figures/solution_space_absolute_expression_vs_data_integration.pdf']);

%% make the FVA figures for key reactions 

% show the flux and boundary for each reaction;
% to display OFD and FVA under the same scale, we show the raw OFD flux instead of normalized flux 

% (1) the TCA cycle: only full integration predicts the full cyclic flux in
% TCA (FVA high conf)
TCArxns = {'RMC0001','RM00267','RMC0003','RM00432','RM02164','RM01082','RM00342','RM00351'};
% not carrying large flux: 'RM00709','RMC0183','RM00405'
% Assuming 'height', 'ub', 'lb', and 'text' are 2D matrices
% with each row representing a data point and each column a group
height = flux_cmp{TCArxns,["OFD_exp_resp_simi","OFD_exp_resp","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}; % Heights of the bars
ub_OFD = flux_cmp{TCArxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]); % Upper bounds
lb_OFD = flux_cmp{TCArxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]); % Lower bounds
ub_PFD = flux_cmp{TCArxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]); % Upper bounds
lb_PFD = flux_cmp{TCArxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]); % Lower bounds
mytext = strcat(TCArxns,':',{' '}, printRxnFormula(model,TCArxns,0)); % Annotations

h = plotReactionSet(height, ub_PFD, lb_PFD, ub_OFD, lb_OFD, mytext);
h.Position(4) = 100*length(TCArxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_TCA.pdf');

% (2) pentose phosphate cycle: full integration predicts high flux (in the
% contrary, expression only uses the folate cycle associated with high 
% serine synthesis [we have serine data]) (FVA high conf)
PPPrxns = {'RC02736','RC02035','RC01528','RC01529','RC01641','RC01830','RC01827','RC03321','RC02740','RC02739'};
h = plotReactionSet(flux_cmp{PPPrxns,["OFD_exp_resp_simi","OFD_exp_resp","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                    flux_cmp{PPPrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{PPPrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{PPPrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{PPPrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    strcat(PPPrxns,':',{' '}, printRxnFormula(model,PPPrxns,0)));% Annotations
% to domonstrate the cycling flux of gpi-1, inhibit f6p --> g6p-A interconversion (RC02740)
% and check FVA interval again (in which case we check if f6p has to go
% through the reverse flux of gpi-1 reaction to recycle g6p or has g6p to
% be recycled
h.Position(4) = 100*length(PPPrxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_PPP.pdf');


% (3) propionate shunt: full integration predicts high flux (FVA high conf)
PPshuntrxns = {'RM04432','RM03045','RM03158','RM01608','RM00705','RM00706'};
h= plotReactionSet(flux_cmp{PPshuntrxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                   flux_cmp{PPshuntrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{PPshuntrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{PPshuntrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{PPshuntrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                  strcat(PPshuntrxns,':',{' '}, printRxnFormula(model,PPshuntrxns,0)));% Annotations
h.Position(4) = 100*length(PPshuntrxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_PPshunt.pdf');


% (4) mthgxl: full integration predicts high secrection - this needs
% investigation and reasoning 
mthgxlrxns = {'RCC0131','RC01016','TCE0660','RM02530','RC02531'};
h = plotReactionSet(flux_cmp{mthgxlrxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                   flux_cmp{mthgxlrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{mthgxlrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{mthgxlrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{mthgxlrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                  strcat(mthgxlrxns,':',{' '}, printRxnFormula(model,mthgxlrxns,0)));% Annotations
h.Position(4) = 100*length(mthgxlrxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_mthgxl.pdf');



% (5) RNA degradation: full integration predicts high degradation of both
% ribose and the base ring, featured the beta-alanine production from ura
RNArxns = {'DGR0002','RC01876','RC01057','RC00978','RC02269','RC00905','RM00908'};
h = plotReactionSet(flux_cmp{RNArxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                     flux_cmp{RNArxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{RNArxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{RNArxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{RNArxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                  strcat(RNArxns,':',{' '}, printRxnFormula(model,RNArxns,0)));% Annotations
h.Position(4) = 100*length(RNArxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_RNA_degradation.pdf');


% (6) high PDH flux: although pep--> pyr is limited, full integration
% predicts high PDH flux
pyrrxns = {'RM00209','RC00200'};
h = plotReactionSet(flux_cmp{pyrrxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                   flux_cmp{pyrrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{pyrrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{pyrrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{pyrrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                  strcat(pyrrxns,':',{' '}, printRxnFormula(model,pyrrxns,0)));% Annotations
h.Position(4) = 100*length(pyrrxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_PDH.pdf');



% (7) high h2o2 production: reflected in high h2o2 secrection, full
% integration predicts high h2o2 production driven by peroxisomal beta
h2o2rxns = {'RCC0090','RCC0165','EX00027'};
h = plotReactionSet(flux_cmp{h2o2rxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                   flux_cmp{h2o2rxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{h2o2rxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{h2o2rxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{h2o2rxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                 strcat(h2o2rxns,':',{' '}, printRxnFormula(model,h2o2rxns,0)));% Annotations
h.Position(4) = 100*length(h2o2rxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_h2o2.pdf');


% (8) DNA deg: DNA deg produced acald and eventually produce ac and used as
% accoa
dnarxns = {'DGR0003','RC01663','RC02484','RC02749','RC01066','RC00711','RM00710','RC00710','RCC0163','RC00235'};
h = plotReactionSet(flux_cmp{dnarxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                   flux_cmp{dnarxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{dnarxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{dnarxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{dnarxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                  strcat(dnarxns,':',{' '}, printRxnFormula(model,dnarxns,0)));% Annotations
h.Position(4) = 100*length(dnarxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_DNA_degradation.pdf');


% (8) mito beta: full integration does not use mito beta
mitoBetarxns = {'RM01282','RM04738','RM03990','RM03777','RMC0107','RM00238','RMC0090'};
h = plotReactionSet(flux_cmp{mitoBetarxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                   flux_cmp{mitoBetarxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{mitoBetarxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{mitoBetarxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{mitoBetarxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                  strcat(mitoBetarxns,':',{' '}, printRxnFormula(model,mitoBetarxns,0)));% Annotations
h.Position(4) = 100*length(mitoBetarxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_mito_beta_oxi.pdf');


% (9) intuitive cases: coA syn ('DMN0044'), urea secraction, mito production 
% of co2 ('TCM1287') and a lot of reactions that were
% trapped at epsilon flux likely due to a responsive gene constraint

% trivial mistake by the expression only: pp canonical is bounded; one of
% the purine synthesis is bounded ('RC04560', but secrect right after); 
otherrxns = {'DMN0044','TCM1287','RC04560','RM01859','RM02765','RM00833'};
h = plotReactionSet(flux_cmp{otherrxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                   flux_cmp{otherrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{otherrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{otherrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{otherrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
               strcat(otherrxns,':',{' '}, printRxnFormula(model,otherrxns,0)));% Annotations
h.Position(4) = 100*length(otherrxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_others.pdf');



% (11) glycine cleavage direction predictions
glyRxns = {'RM03425','RM04125','RM03815'};
h = plotReactionSet(flux_cmp{glyRxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                    flux_cmp{glyRxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{glyRxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{glyRxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{glyRxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                  strcat(glyRxns,':',{' '}, printRxnFormula(model,glyRxns,0)));% Annotations
h.Position(4) = 100*length(glyRxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_glycine_cleavage.pdf');


% (11) novel predictions
novelRxns = {'RM00355'};
h = plotReactionSet(flux_cmp{novelRxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                   flux_cmp{novelRxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{novelRxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{novelRxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{novelRxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                  strcat(novelRxns,':',{' '}, printRxnFormula(model,novelRxns,0)));% Annotations
h.Position(4) = 100*length(novelRxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_others2.pdf');


% non-specific predictions
% (13) FA synthesis:  pmtcoa is synthesized instead of burnt bc bacterial pmtcoa 
% is not enough for worm lipid biomass (but both integration predicted)
pmtRxns = {'RCC0020', 'RC07758', 'BIO0006', 'RCC0165'};
h = plotReactionSet(flux_cmp{pmtRxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                  flux_cmp{pmtRxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{pmtRxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{pmtRxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{pmtRxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                  strcat(pmtRxns,':',{' '}, printRxnFormula(model,pmtRxns,0)));% Annotations
h.Position(4) = 100*length(pmtRxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_FA_syn.pdf');


% (14) high flux in mevelonate syn 'RC02082'
mevRxns = {'RC02082'};
h = plotReactionSet(flux_cmp{mevRxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                   flux_cmp{mevRxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{mevRxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{mevRxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{mevRxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                  strcat(mevRxns,':',{' '}, printRxnFormula(model,mevRxns,0)));% Annotations
h.Position(4) = 100*length(mevRxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_mevalonate_syn.pdf');

% (15) mito accoa source
accoaRxns = {'RMC0106','RMC0088','RM00927','RMC0090','RM04747','RM01177','RM00238','RMC0107','RM00209','RM00705','RMC0090','RM03858','RM04742','RM03778'};
h = plotReactionSet(flux_cmp{accoaRxns,["OFD_exp_resp_simi","OFD_exp_resp_simi","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                   flux_cmp{accoaRxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{accoaRxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{accoaRxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{accoaRxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                  strcat(accoaRxns,':',{' '}, printRxnFormula(model,accoaRxns,0)));% Annotations
h.PaperSize = [8.5 20];
h.PaperPosition = [0.05 0.05 8.3 18];  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_accoa_syn.pdf');

% (16) met/sam cycle
MetSAMrxns = {'RC00177','RC00192','RC01290','RC00946'};
h = plotReactionSet(flux_cmp{MetSAMrxns,["OFD_exp_resp_simi","OFD_exp_resp","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                    flux_cmp{MetSAMrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{MetSAMrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{MetSAMrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{MetSAMrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    strcat(MetSAMrxns,':',{' '}, printRxnFormula(model,MetSAMrxns,0)));% Annotations
h.Position(4) = 100*length(MetSAMrxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_Met_SAM_cycle.pdf');

% (17) glycolysis cycle
glycolysisrxns = {'RC01786','RC01600','RC03321','RC02740','RC04779','RC01070','RC01015',...
                    'RC01061','RC01512','RC01518','RC00658','RC00200','RC00431'};
h = plotReactionSet(flux_cmp{glycolysisrxns,["OFD_exp_resp_simi","OFD_exp_resp","OFD_exp_simi","OFD_exp_only","OFD_no_data"]}, ...% Heights of the bars
                    flux_cmp{glycolysisrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{glycolysisrxns,["FVA_PFD_exp_resp_simi","FVA_PFD_exp_resp","FVA_PFD_exp_simi","FVA_PFD_exp_only","FVA_PFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    flux_cmp{glycolysisrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[2,4,6,8,10]),...% Upper bounds
                    flux_cmp{glycolysisrxns,["FVA_OFD_exp_resp_simi","FVA_OFD_exp_resp","FVA_OFD_exp_simi","FVA_OFD_exp_only","FVA_OFD_no_data"]}(:,[1,3,5,7,9]),...% Lower bounds
                    strcat(glycolysisrxns,':',{' '}, printRxnFormula(model,glycolysisrxns,0)));% Annotations
h.Position(4) = 100*length(glycolysisrxns)+100;  % for example, set the height to 800
saveas(h, 'figures/flux_comparison_glycolysis.pdf');

% --> the two constraints, responsiveness and similarity synergistically
% narrow down the solution space (if we claim this, we need a simulation
% with only exp + similarity)
%% metabolite splits
met = 'accoa[m]'; % 'accoa[c]' leu-L
[h1 h2] = plotContribution(model,myCSM_no_data.OFD,met); % myCSM_merged.OFD
saveas(h1, 'figures/flux_contribution_accoa_prod_no_data.pdf');
saveas(h2, 'figures/flux_contribution_accoa_cons_no_data.pdf');
[h1 h2] = plotContribution(model,myCSM_exp_only.OFD,met); % myCSM_merged.OFD
saveas(h1, 'figures/flux_contribution_accoa_prod_exp_only.pdf');
saveas(h2, 'figures/flux_contribution_accoa_cons_exp_only.pdf');
[h1 h2] = plotContribution(model,myCSM_exp_resp.OFD,met); % myCSM_merged.OFD
saveas(h1, 'figures/flux_contribution_accoa_prod_exp_resp.pdf');
saveas(h2, 'figures/flux_contribution_accoa_cons_exp_resp.pdf');
[h1 h2] = plotContribution(model,myCSM_exp_simi.OFD,met); % myCSM_merged.OFD
saveas(h1, 'figures/flux_contribution_accoa_prod_exp_simi.pdf');
saveas(h2, 'figures/flux_contribution_accoa_cons_exp_simi.pdf');
[h1 h2] = plotContribution(model,myCSM_exp_resp_simi.OFD,met); % myCSM_merged.OFD
saveas(h1, 'figures/flux_contribution_accoa_prod_exp_resp_simi.pdf');
saveas(h2, 'figures/flux_contribution_accoa_cons_exp_resp_simi.pdf');

%% codes for manual interaction 
% %% data inspection and identification of point of interest
% met = 'pyr[m]'; % 'accoa[c]' leu-L
% tbl3 = listRxn(model,myCSM_exp.OFD,met); % myCSM_merged.OFD
% 
% %% data visualization
% % it seems in merged integration, many reactions in non-core metabolism
% % (such as pai3p) is unbounded regardless that they are bounded in
% % expression integration; it seems the merged on is good at core metabolism
% % but bad at the non-core.
% close all
% %% flux heatmap 
% fluxMat = flux_cmp{:,3:7};
% fluxMat = normalize(fluxMat', "scale",max(abs(fluxMat),[],2)')';
% fluxMat(isnan(fluxMat)) = 0;
% clustergram(fluxMat, 'Cluster','column')

%% make the figures for data basis of the key predictions - discontinued - use LOO
% addpath scripts/parfor_wait/
% targetRxn = 'EX00027'; % EX00027, RMC0090 RC02749 RC00200
% parforFlag = 1;
% relMipGapTol = 1e-3; % this is to be redone with maximum precision, otherwsie the box will have problem
% verbose = true;
% 
% % run FVA by calling:
% [minval, maxval] = perturb_constraints(myCSM_triple.MILP, model, targetRxn,parforFlag,relMipGapTol,verbose);
% %%
% % if no hit is returned, it may be either multiple constraints or variable
% % boundaries.
% ind = find(minval == 0);
% % first look at relevant metabolites
% model.mets(ind(ind <= length(model.mets)))
% % ==> often difficult to understand 
% 
% % next look at those constriants
% ind = ind(ind > length(model.mets));
% 
% RHnames = model.rxns(find(ismember(model.rxns, myCSM_triple.RHNames)));
% RLnames = model.rxns(find(ismember(model.rxns, myCSM_triple.RLNames)));
% 
% c_v = model.rxns;
% c_yh1 = strcat('RHfit_1_', RHnames);
% c_yl = RLnames;
% c_yh2 = strcat('RHfit_2_', RHnames);
% c_minFluxLowRxns = strcat('abs_flux_', model.rxns);
% c_gene = myCSM_triple.HGenes;
% c_fluxCoupling = myCSM_triple.branchMets;
% latentVar = [strcat('Latentfit_1_', myCSM_triple.latentRxn);strcat('Latentfit_2_', myCSM_triple.latentRxn)];
% variables = [c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns;c_fluxCoupling;latentVar];
% 
% txt = printConst(myCSM_triple.MILP,ind, variables);
% 
% % ==> too many; it seems loosing many constraints alone can create a large
% % feasible space for the target rxns to go zero. It is likely becuase when
% % one constraint is removed, the caps (like minLow, metFit, or MILP obj)
% % were not updated. we may have to change on input data (like remove a
% % responsive gene) at a time and redo the entire thing to check...
% 
% % first check for more rxns and more FVA interval threshold 
% 
% % ==> this methold is misleading as the constraints like minLow/MetFit
% % need to be updated to enable meaningful analysis; otherwise any big flux
% % related constraints can be the hits in the screen (they make room for the
% % target flux to be flexible).
% 
% %%
% targetRxns = {'RC00200'};
% parforFlag = 0;
% relMipGapTol = 1e-3; % this is to be redone with maximum precision, otherwsie the box will have problem
% verbose = false;
% 
% % run FVA by calling:
% [minval, maxval] = FVA_MILP(myCSM_triple.MILP, model, targetRxns,parforFlag,relMipGapTol,verbose)
% 
% 
