%% add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
addpath ./../../MetabolicLibrary/7_FBA_modeling/PlotPub/lib/

load('output/sampling_results.mat')
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat'); % model 
% setup the model
model = configurateModel(model);
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

%% note
% since this sampling analysis re-executed the iMAT-WPS in different
% environment (cluster vs server). there is numerical difference in the
% output result; the objective values of each optimization remains
% identical given the numerical precision of the solver but the flux
% distribution differs minorly (alternative distributions under the same
% optimization precision). We will present the data as is and ignore the
% numenric difference in the flux values and number of constrained reaction
% when doing 100% subsampling compared with the iMAT-WPS reference result.
ori_pred = readtable('output\fluxTable.csv');
sigFlux = 1e-5; % the significant flux threshold to determine bounded vs. unbounded reactions.

%% visualize the OFD and FVA bounds distribution with sampling 

% start with examples
targetRxns = {'RC02736','RC02035','RC01528','RC01529','RC01641','RC01830','RC01827','RC03321',... % PPP reactions
            'RMC0001','RM00267','RMC0003','RM00432','RM02164','RM01082','RM00342','RM00351',... TCA 
            'RM04432','RM03045','RM03158','RM01608','RM00705','RM00706', ... PP shunt
            'RC04560', 'RC01290' % purine and met/sam
            };

for zz = 1:length(targetRxns)
    myrxn = targetRxns{zz};
    ori_flux = ori_pred.normalized_OFD_exp_resp_simi(strcmp(ori_pred.rxns,myrxn));

    % for numerical stability, mask insignificant flux 
    ori_flux(abs(ori_flux) < sigFlux) = 0;

    % grab the OFD flux values 
    flux = [];
    PFM_ub = [];
    PFM_lb = [];
    OFM_lb = [];
    OFM_ub = [];
    depth = [];
    for i = 1:length(samplingDepth)
        for j = 1:samplingReps
            flux = [flux; samples_CSM{i,j}.OFD(strcmp(model.rxns,myrxn)) ./ abs(samples_CSM{i,j}.OFD(strcmp(model.rxns,'EXC0050')))];
            PFM_lb = [PFM_lb; samples_PFM_sol{i,j}(strcmp(model.rxns,myrxn),1)];
            PFM_ub = [PFM_ub; samples_PFM_sol{i,j}(strcmp(model.rxns,myrxn),2)];
            OFM_lb = [OFM_lb; samples_OFM_sol{i,j}(strcmp(model.rxns,myrxn),1)];
            OFM_ub = [OFM_ub; samples_OFM_sol{i,j}(strcmp(model.rxns,myrxn),2)];
            depth = [depth; samplingDepth(i)];
        end
    end
    % for numerical stability, mask insignificant flux 
    flux(abs(flux) < sigFlux) = 0;
    
    figure;
    subplot('Position', [0.1, 0.4, 0.8, 0.55]); % [left, bottom, width, height]
    % Create the boxplot
    boxplot(flux, depth)
    
    % Hold on to overlay the data points
    hold on
    
    % Plot the data points with jitter for better visibility
    jitterAmount = 0.1; % Adjust the jitter amount as needed
    group = unique(depth);
    recall = [];
    for i = 1:length(group)
        currentGroup = flux(depth == group(i));
        % calculate recall 
        delta_flux = (currentGroup - ori_flux)/ori_flux;
        delta_flux(isnan(delta_flux)) = 0; % 0 to 0 is recalled
        recall(i) = sum(delta_flux < 0.3 & delta_flux > -0.3) / length(delta_flux);
        x = repmat(i, size(currentGroup)) + (rand(size(currentGroup)) - 0.5) * jitterAmount;
        scatter(x, currentGroup, 'filled', 'k')
    end
    
    % Add labels and title if needed
    xlabel('Subsample depth (% WPS tested gene)')
    ylabel('Flux')
    title(myrxn)
    
    % Release the hold
    hold off
        
    subplot('Position', [0.1, 0.1, 0.8, 0.25]); % [left, bottom, width, height]
    line(categorical(unique(depth)), recall);
    ylim([0,1])
    xlabel('Subsample depth (% WPS tested gene)')
    ylabel('% Times original prediction is recalled')

    saveas(gca,['figures/sensitivity_analysis/OFD_',myrxn,'.pdf']);
end

% figure;
% % plot boundaries too
% boxplot(PFM_lb, depth)
% hold on
% jitterAmount = 0.1; % Adjust the jitter amount as needed
% group = unique(depth);
% for i = 1:length(group)
%     currentGroup = PFM_lb(depth == group(i));
%     x = repmat(i, size(currentGroup)) + (rand(size(currentGroup)) - 0.5) * jitterAmount;
%     scatter(x, currentGroup, 'filled')
% end
% % Add labels and title if needed
% xlabel('Subsample depth (% WPS tested gene)')
% ylabel('FVA lower bound')
% title(myrxn)
% % Release the hold
% hold off
% 
% figure;
% boxplot(PFM_ub, depth)
% hold on
% jitterAmount = 0.1; % Adjust the jitter amount as needed
% group = unique(depth);
% for i = 1:length(group)
%     currentGroup = PFM_ub(depth == group(i));
%     x = repmat(i, size(currentGroup)) + (rand(size(currentGroup)) - 0.5) * jitterAmount;
%     scatter(x, currentGroup, 'filled')
% end
% % Add labels and title if needed
% xlabel('Subsample depth (% WPS tested gene)')
% ylabel('FVA upper bound')
% title(myrxn)
% % Release the hold
% hold off


% ==> too busy to visualize the bounds so we skip bounds and just focus on
% flux prediction for now

%% visualize the change of directionally constrained reactions 
PFD_bounded_ori = ori_pred.rxns(ori_pred.PFD_bounded_exp_resp_simi==1);
OFD_bounded_ori = ori_pred.rxns(ori_pred.OFD_bounded_exp_resp_simi==1);


N_OFM_bounded = [];
N_PFM_bounded = [];
Jaccard_OFM_bounded = [];
Jaccard_PFM_bounded = [];
depth = [];
for i = 1:length(samplingDepth)
    for j = 1:samplingReps

        % annotate the FVA boundedness 
        PFD_bounded  =  (samples_PFM_sol{i,j}(:,2) < -sigFlux |... % negtive flux
                         samples_PFM_sol{i,j}(:,1) > sigFlux |... % positive flux 
                        (samples_PFM_sol{i,j}(:,1) > -sigFlux & samples_PFM_sol{i,j}(:,2) < sigFlux)); % zero flux 
        OFD_bounded  =  (samples_OFM_sol{i,j}(:,2) < -sigFlux |... % negtive flux
                         samples_OFM_sol{i,j}(:,1) > sigFlux |... % positive flux 
                        (samples_OFM_sol{i,j}(:,1) > -sigFlux & samples_OFM_sol{i,j}(:,2) < sigFlux)); % zero flux 
        
        N_PFM_bounded = [N_PFM_bounded; sum(PFD_bounded)];
        N_OFM_bounded = [N_OFM_bounded; sum(OFD_bounded)];

        Jaccard_PFM_bounded = [Jaccard_PFM_bounded; length(intersect(model.rxns(PFD_bounded), PFD_bounded_ori))/length(union(model.rxns(PFD_bounded), PFD_bounded_ori))];
        Jaccard_OFM_bounded = [Jaccard_OFM_bounded; length(intersect(model.rxns(OFD_bounded), OFD_bounded_ori))/length(union(model.rxns(OFD_bounded), OFD_bounded_ori))];

        depth = [depth; samplingDepth(i)];
    end
end

% the absolute number of constrained reactions is sensitive to modeling
% noise and numerical stability (like extreme maxmization or minimization);
% for sensitivity evaulation, we focus on the overlap (recall) with the
% original full dataset prediction

% % Create the boxplot
% boxplot(N_OFM_bounded, depth)
% 
% % Hold on to overlay the data points
% hold on
% 
% % Plot the data points with jitter for better visibility
% jitterAmount = 0.1; % Adjust the jitter amount as needed
% group = unique(depth);
% for i = 1:length(group)
%     currentGroup = N_OFM_bounded(depth == group(i));
%     x = repmat(i, size(currentGroup)) + (rand(size(currentGroup)) - 0.5) * jitterAmount;
%     scatter(x, currentGroup, 'filled','r')
% end
% 
% 
% % Create the boxplot
% boxplot(N_PFM_bounded, depth)
% % Plot the data points with jitter for better visibility
% jitterAmount = 0.1; % Adjust the jitter amount as needed
% group = unique(depth);
% for i = 1:length(group)
%     currentGroup = N_PFM_bounded(depth == group(i));
%     x = repmat(i, size(currentGroup)) + (rand(size(currentGroup)) - 0.5) * jitterAmount;
%     scatter(x, currentGroup, 'filled','b')
% end
% 
% % Add labels and title if needed
% xlabel('Subsample depth (% WPS tested gene)')
% ylabel('Number of directionally constrained reactions')
% 
% % Release the hold
% hold off

% plot jaccard index 
figure;
% Create the boxplot
boxplot(Jaccard_OFM_bounded, depth, 'Colors','rrrr')
% Hold on to overlay the data points
hold on
% Plot the data points with jitter for better visibility
jitterAmount = 0.1; % Adjust the jitter amount as needed
group = unique(depth);
for i = 1:length(group)
    currentGroup = Jaccard_OFM_bounded(depth == group(i));
    x = repmat(i, size(currentGroup)) + (rand(size(currentGroup)) - 0.5) * jitterAmount;
    scatter(x, currentGroup, 'filled','r')
end
% Add labels and title if needed
xlabel('Subsample depth (% WPS tested gene)')
ylabel('Jaccard index of OFM directionally constrained reactions')
% refline((max(Jaccard_OFM_bounded)-min(Jaccard_OFM_bounded))/(length(samplingDepth) - 1), min(Jaccard_OFM_bounded)-(max(Jaccard_OFM_bounded)-min(Jaccard_OFM_bounded))/(length(samplingDepth) - 1))

% Create the boxplot
boxplot(Jaccard_PFM_bounded, depth, 'Colors','bbbb')
% Plot the data points with jitter for better visibility
jitterAmount = 0.1; % Adjust the jitter amount as needed
group = unique(depth);
for i = 1:length(group)
    currentGroup = Jaccard_PFM_bounded(depth == group(i));
    x = repmat(i, size(currentGroup)) + (rand(size(currentGroup)) - 0.5) * jitterAmount;
    scatter(x, currentGroup, 'filled','b')
end
% Add labels and title if needed
xlabel('Subsample depth (% WPS tested gene)')
ylabel('Jaccard index of PFM directionally constrained reactions')
hold off
% refline((max(Jaccard_PFM_bounded)-min(Jaccard_PFM_bounded))/(length(samplingDepth) - 1), min(Jaccard_PFM_bounded)-(max(Jaccard_PFM_bounded)-min(Jaccard_PFM_bounded))/(length(samplingDepth) - 1))
ylim([min([Jaccard_PFM_bounded; Jaccard_OFM_bounded])*0.95,1])
% Release the hold

saveas(gca,['figures/sensitivity_analysis/sensitivity_of_directionally_constrained_reactions.pdf']);


%% score the robust predictions
% we define a sensitivity score as: at 80% subsampling depth of the WPS
% data, how many time a prediction is recalled as the full dataset
% prediction (normalzied OFD = OFD_ori +/- 30%).
at_depth = 0.8;
    
% recall matrix
recall = [];
for zz = 1:length(model.rxns)
    myrxn = model.rxns{zz};
    ori_flux = ori_pred.normalized_OFD_exp_resp_simi(strcmp(ori_pred.rxns,myrxn));

    % for numerical stability, mask insignificant flux 
    ori_flux(abs(ori_flux) < sigFlux) = 0;

    % grab the OFD flux values 
    flux = [];
    depth = [];
    for i = 1:length(samplingDepth)
        for j = 1:samplingReps
            flux = [flux; samples_CSM{i,j}.OFD(strcmp(model.rxns,myrxn)) ./ abs(samples_CSM{i,j}.OFD(strcmp(model.rxns,'EXC0050')))];
            depth = [depth; samplingDepth(i)];
        end
    end
    % for numerical stability, mask insignificant flux 
    flux(abs(flux) < sigFlux) = 0;
    
    % calculate recall
    group = unique(depth);
    for i = 1:length(group)
        currentGroup = flux(depth == group(i));
        % calculate recall 
        delta_flux = (currentGroup - ori_flux)/ori_flux;
        delta_flux(isnan(delta_flux)) = 0; % 0 to 0 is recalled
        recall(i,zz) = sum(delta_flux < 0.3 & delta_flux > -0.3) / length(delta_flux);
    end
end

% visualize the recall for all bounded reactions originally 
% (other reactions are more of a result of flux minimization so we dont
% wnat to over interpret them)
PFD_bounded_ori = ori_pred.rxns(ori_pred.PFD_bounded_exp_resp_simi==1);
cgo=clustergram(recall(:,ismember(model.rxns,PFD_bounded_ori)),'Cluster', 'row','RowLabels',group,'Symmetric',false);
c=get(cgo,'ColorMap');
n = 100;
tmp = [ones(n,1), linspace(1,0,n)',linspace(1,0,n)'];
tmp = tmp(2:end,:);
cpr=[linspace(0,1,n)',linspace(0,1,n)',ones(n,1);...
    tmp];
set(cgo,'ColorMap',cpr);
set(cgo,'Symmetric',false);
set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 12)
% save
set(h, 'FontSize', 7)
exportgraphics(gcf,'figures/sensitivity_analysis/clustergram_OFD_recalls.tiff','ContentType','vector')
exportgraphics(gcf,'figures/sensitivity_analysis/clustergram_OFD_recalls.pdf','ContentType','vector')

% observation:
% for predictions specific to WPS integration (not recalled at depth of 0)
% there is a spectrum of sensitivity; some predictions breaks very early
% (i.e., 50% recall breaks at even 90% subsampling) but some breaks very
% late. there is a full continuum of sensitivity.

% however, at global scale, the OFD flux of most reactions is insensitive
% to small perturbations of the WPS dataset, aligned with the multi-layer
% integration framework. we will either use the heatmap or use the
% histogram at fixed depth for the paper

% for a representitive showcase, we arbiturily use the recall at 80%
% subsampling as a metric

% plot the distrobution of recall score 
recall80 = recall(group==0.8,ismember(model.rxns,PFD_bounded_ori));
figure;
histogram(recall80,'Normalization','probability')
xlabel('% Times original prediction is recalled')
ylabel('Probability')
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [3.5, 3];
plt.LineWidth = 0.5;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/sensitivity_analysis/OFD_recalls_80_perc_subsample.pdf']);

% of note, many of these highly recalled reaction fluxes (i.e,> 0.9) are in
% fact the fluxes that can be just predicted by the expression integration
% and is also supported by the WPS data (adding WPS did not change the
% prediction). They are still a type of robustness but of course are not
% the uniquely enabled predictions by WPS dataset (but many are uniquely
% directionally constrained in WPS integration although the two
% integrations indicate the same OFD flux).

% also add the sensitivity score to the table

tbl_pred = readtable('output/fluxTable_putative_mechanism_annotated.csv');
[A B] = ismember(tbl_pred.rxns, model.rxns);
tbl_pred.sensitivity_score = recall(group==0.8,B(A))';
writetable(tbl_pred,'output/fluxTable_putative_mechanism_annotated.csv');


% check the extreme case - exclude those reactions fully recalled as OFD in
% exp-only integration and only look at the sensitivity of reactions knwon
% to be sensitive (changing a lot along titration)
% PFD_bounded_exp = model.rxns(recall(group==0,:) > 0.9);
% 
% % plot the distrobution of recall score 
% recall80 = recall(group==0.8,ismember(model.rxns,setdiff(PFD_bounded_ori, PFD_bounded_exp)));
% figure;
% histogram(recall80,'Normalization','probability')
% xlabel('% Times original prediction is recalled')
% ylabel('Probability')
% plt = Plot(); % create a Plot object and grab the current figure
% plt.BoxDim = [3.5, 3];
% plt.LineWidth = 0.5;
% plt.FontSize = 7;
% plt.FontName = 'Arial';
% plt.ShowBox = 'off';
% plt.XMinorTick = 'off';
% plt.YMinorTick = 'off';
% plt.TickDir = 'out';

% observation: it is a spectrum of sensitivity even for these reactions:
% there are reactions greater than 0.8, around 0.8 and lower than 0.8.
% Those greater than 0.8 are fluxes redundantly predicted by multiple WPS
% data thus is insensitive to small perturbations of the dataset; those
% around 0.8 are likely driven by a single gene data so when the gene data
% is gone, it cannot be predicted. Those lower than 0.8 are predictions
% driven by the existance of multiple coordianted constraints where the
% loss of any of them can significantly influence the flux. They all makes
% sense and the ones that is greater than 0.8 (eg >0.9) is slightly more
% than the others, indicating a good robustness. All together, we are
% convinced that the predictions are relatively robust but we will
% certainly avoid overemphasizing the robustness. 


