function cluster_a2_prediction_mechanism_analysis_simi(tmpDir,batchID)
%% about 
% perform the Leave One Out and Leave One In analysis of the responsiveness
% integration. 

% parameters
sigFlux = 1e-5; 
% we use the disassimilation (yield) rate to constrain the bacteria waste
disAssim = 0.65;

%% Step 1: make the OFDs
% add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/10.0.0/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
addpath integration_pipelines/

load([tmpDir,'/environment.mat']);
restoreEnvironment(environment);

%% Load model
load('./input/model/makeWormModel/iCEL1314_withUptakes.mat');
load('./input/model/epsilon_generic_withUptakes.mat'); % see walkthrough_generic.m for guidance on generating the epsilon values
load('input/WPS/categ_expression_and_WPS.mat');
branchTbl = readtable('input/WPS/final_branchPoint_table.csv'); % must have 'mets', 'rxn1','rxn2', and 'maxCosine'

% setting up the model 
model = configurateModel(model);
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

% Set parameters
modelType = 2; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
speedMode = 1;
%% gather stuffs to analyze mechanisms for 
predTbl = readtable('output/fluxTable.csv');
simi_depend = predTbl.rxns(predTbl.PFD_bounded_exp_simi & ~predTbl.PFD_bounded_exp_only);
predTbl.related_responsiveness_constraints = repmat({'ND'}, size(predTbl,1),1);
predTbl.related_similarity_constraints = repmat({'ND'}, size(predTbl,1),1);

%% analyze the mechanisms dependent on similarity
allMets = unique(branchTbl.mets);

model_coupled = model;
% add the disassimilation constraints 
model_coupled.S(end+1,:) = zeros(1,length(model_coupled.rxns));
model_coupled.S(end, strcmp('EXC0050',model_coupled.rxns)) = disAssim; 
model_coupled.S(end, strcmp('BIO0010',model_coupled.rxns)) = 1; 
model_coupled.csense(end+1) = 'G';
model_coupled.b(end+1) = 0;
model_coupled.BiGG(end+1) = {'NA'};
model_coupled.metCharges(end+1) = 0;
model_coupled.metFormulas(end+1) = {'NA'};
model_coupled.metNames(end+1) = {'Pseudo-metabolite that represents a constraint'};
model_coupled.mets(end+1) = {['NonMetConst',num2str(length(model_coupled.mets))]};

% remove the responsiveness
ExpCateg.responsive = {};
ExpCateg.nonresponsive = {};

parforFlag = 0;
relMipGapTol = 1e-12; % this is to be redone with maximum precision, otherwsie the box will have problem
verbose = false;
    

for zz = str2num(batchID)

    targetRxn = simi_depend{zz};  % EX00027, RMC0090 RC02749 RC00200 RM00209 RC00905 RM01082
    
    fprintf('%d: %s\n', zz, targetRxn);

    for yy = 1:2
        % OFDs = [];
        Maxs = [];
        Mins = [];
        for ii = 1:(size(allMets,1)+1)

            branchTbl_tmp = branchTbl;
            if ii > 1
                if yy == 1 % perform LOO analysis
                    branchTbl_tmp(strcmp(branchTbl_tmp.mets,allMets{ii-1}),:) = [];
                elseif yy == 2 % perform Leave One In (LOI) analysis
                    branchTbl_tmp(~strcmp(branchTbl_tmp.mets,allMets{ii-1}),:) = [];
                end
            end
        
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
            myCSM.branchMets]...
            = IMATplusplus_wiring_triple_inetgration_final(model_coupled,epsilon_f,epsilon_r, ExpCateg, branchTbl_tmp, modelType,speedMode,...
            0.05, 0.05, 0, 0, 0.05, [size(model.S,1)-1 0.01],[size(model.S,1) 0.01],10);
            % skip OFD modeling; set the caps to 5%

            [minval, maxval] = FVA_MILP(myCSM.MILP_PFD, model_coupled, {targetRxn},parforFlag,relMipGapTol,verbose);
            
            % OFDs(:,ii) = myCSM.OFD;
            Maxs(ii) = maxval;
            Mins(ii) = minval;

            if ii > 1
                if yy == 1
                    % fprintf('removing %s: v = %f, FVA = [%f, %f]\n', allMets{ii-1},myCSM.OFD(strcmp(model.rxns,targetRxn)), minval, maxval);
                    fprintf('removing %s: FVA = [%f, %f]\n', allMets{ii-1}, minval, maxval);
                elseif yy == 2
                    % fprintf('preserving %s: v = %f, FVA = [%f, %f]\n', allMets{ii-1},myCSM.OFD(strcmp(model.rxns,targetRxn)), minval, maxval);
                    fprintf('preserving %s: FVA = [%f, %f]\n', allMets{ii-1}, minval, maxval);
                end
            end
        
        end
        
    
        % find out the driver constraint
        max_ori = Maxs(1);
        min_ori = Mins(1);
        % first confirm the type of bounded rxn 
        if min_ori > sigFlux 
            boundedType = 'forward_flux';
        elseif max_ori < -sigFlux
            boundedType = 'revserse_flux';
        elseif max_ori < sigFlux && min_ori > -sigFlux 
            boundedType = 'zero_flux';
        else
            error('FVA not reproduced!')
        end

        % remove the first one that is to reproduce the orginal FVA
        Maxs(1) = [];
        Mins(1) = [];

        if strcmp(boundedType, 'forward_flux')
            if yy == 1 % perform LOO analysis
                % check when forward flux can go to zero or lower
                driverMets_LOO = allMets(Mins < sigFlux);
            elseif yy == 2 % perform LOI analysis 
                % check when forward flux can be bounded
                driverMets_LOI = allMets(Mins > sigFlux);
            else
                error('bug occurs!')
            end
        elseif strcmp(boundedType, 'revserse_flux')
            if yy == 1 % perform LOO analysis
                % check when forward flux can go to zero or lower
                driverMets_LOO = allMets(Maxs > -sigFlux);
            elseif yy == 2 % perform LOI analysis 
                % check when forward flux can be bounded
                driverMets_LOI = allMets(Maxs < -sigFlux);
           else
                error('bug occurs!')
            end
        elseif strcmp(boundedType, 'zero_flux')
            if yy == 1 % perform LOO analysis
                % check when forward flux can go to zero or lower
                driverMets_LOO = allMets(Maxs > sigFlux | Mins < -sigFlux);
            elseif yy == 2 % perform LOI analysis 
                % check when forward flux can be bounded
                driverMets_LOI = allMets(Maxs < sigFlux & Mins > -sigFlux);
            else
                error('bug occurs!')
            end
        else
            error('bug occurs!')
        end
    end

    % save results in human readible texts
    if isempty(driverMets_LOO)
        driverMets_LOO_txt = 'no related metabolite constraints';
    else
        driverMets_LOO_txt = ['the following related metabolites: ',strjoin(driverMets_LOO,', ')];
    end

    if isempty(driverMets_LOI)
        driverMets_LOI_txt = 'no related metabolite constraints';
    else
        driverMets_LOI_txt = ['the following related metabolites: ',strjoin(driverMets_LOI,', ')];
    end

    if isempty(driverMets_LOI) && isempty(driverMets_LOO)
        text = ['Model predicts ', regexprep(boundedType,'_',' '),', but no related constraint is found.'];
    else    
        text = ['Model predicts ', regexprep(boundedType,'_',' '),'. ',...
                'Leave-one-out analysis indicates ', driverMets_LOO_txt,'. ',...
                'Leave-one-in analysis indicates ', driverMets_LOI_txt,'.'];
    end



    predTbl.related_similarity_constraints(strcmp(predTbl.rxns,targetRxn)) = {text};
    
end

% checkTbl = branchTbl(ismember(branchTbl.mets, driverMets),:);
writetable(predTbl, [tmpDir,'/fluxtale_ann_',batchID,'.csv']);

end