addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/10.0.0/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
addpath integration_pipelines/

% also add a no-sampling control to make sure codes are correct


% full WPS dataset 
% 12302023: update to final publish-version data in which four conditions
% was wrong targeted gene annotations (all nonresponsive) were corrected
conditionInfo = readtable('./input/WPS/RNAi_condition_information.csv');
% to be consistent with the model annotation, we still use the mrpl-44 as
% mccc-2
conditionInfo.RNAi_WBID(strcmp(conditionInfo.RNAiID,'x.mrpl_44 met3_lib1')) = {'WBGene00008514'};
conditionInfo.RNAi_geneName(strcmp(conditionInfo.RNAiID,'x.mrpl_44 met3_lib1')) = {'mrpl-44'};

% we also include the phenotype genes in the sampling (but will be
% processed seperately)
pheno = readtable('./input/WPS/not sequenced RNAis because of strong phenotype.csv','ReadVariableNames',1);
tmp = repmat({'PHENOTYPE','PHENOTYPE','PHENOTYPE','PHENOTYPE',...
            'PHENOTYPE', 'PHENOTYPE', 'PHENOTYPE', 'PHENOTYPE', 'PHENOTYPE',...
                'PHENOTYPE','PHENOTYPE','PHENOTYPE','PHENOTYPE','PHENOTYPE',nan}, size(pheno,1),1);
tmp(:,4) = pheno.WBID;
conditionInfo(end+1:end+size(pheno,1),:) = tmp;

% restrict the sampling analysis within the iCEL genes because non-iCEL
% will be always filtered and do not form effective dataset 
lookupTbl = readtable('./input/model/IDtbl.csv');
conditionInfo = conditionInfo(ismember(conditionInfo.RNAi_WBID, lookupTbl.WormBase_Gene_ID),:);

% subsample the dataset
samplingDepth = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0];
samplingReps = 30;

rng(52)
sampleData = {};
for i = 1:length(samplingDepth)
    for j = 1:samplingReps
        % we dont do a strict gene-wise sampling; instead, we do a sampling
        % on the inputy data directly
        sampleData{i,j} = conditionInfo(randperm(size(conditionInfo,1), round(samplingDepth(i) * size(conditionInfo,1))),:);
    end
end


% split the jobs and run
initCobraToolbox(false);
environment = getEnvironment();
tmpDir = ['tmp_',num2str(100+fix(rand()*100000))];
mkdir(tmpDir);
uid = tmpDir(5:6);
save([tmpDir,'/environment.mat'],'environment','sampleData','samplingDepth','-v7.3','-nocompression');

% make the target rxn pools
for i = 1:length(samplingDepth)
    for j = 1:samplingReps    
        batchID = ['b_',uid,'_',num2str(i),'_',num2str(j)];
        env_cmd = 'module load gurobi/10.0.0 && module load matlab/r2022b && ';
        cmd = ['matlab -nodisplay -nosplash -nojvm -r \"cluster_a2_sensitivity_analysis_sampling ',tmpDir,' ',num2str(i),' ',num2str(j),'\"'];
        full_cmd = [env_cmd, cmd, ' > ',tmpDir,'/',batchID,'.log && wait'];
        bsub_cmd = ['bsub -q long -W 24:00 -n 1 -R rusage[mem=10000] -e ',tmpDir,'/err_',batchID,'.log ','-J ',batchID,' "'];
        cmd_ready = [bsub_cmd, full_cmd,'"'];
        system(cmd_ready);
        pause(0.1);
    end
end

%% start job monitoring
fprintf('job pooling finished, pausing...\n');
pause(10); % wait for all jobs to be submitted
fprintf('start job monitoring...\n');
runningBatches = {};
for i = 1:length(samplingDepth)
    for j = 1:samplingReps    
        runningBatches = [runningBatches;{['b_',uid,'_',num2str(i),'_',num2str(j)]}];
    end
end

while ~isempty(runningBatches)
    allFiles = dir([tmpDir,'/iMAT_WPS_subsample_*.mat']);
    allFiles = {allFiles.name};
    finished = regexprep(allFiles,'^iMAT_WPS_subsample_|.mat$','');
    runningBatchesID = setdiff(regexprep(runningBatches,'^b_.._',''),finished);
    runningBatches = strcat(['b_',uid,'_'],runningBatchesID);
    % find the failed runs
    [~,out] = system("bjobs -q long -o 'jobid job_name:50 queue'");
    allTerms = strsplit(out,' ');
    failedBatches = setdiff(runningBatches,allTerms);
    % resubmit the failed batches 
    failedBatches = regexprep(failedBatches,'^b_.._','');
 
    for zzz = 1:length(failedBatches)
        str = failedBatches(zzz);
        strs = strsplit(str{:}, '_');
        batchID = ['b_',uid,'_',str{:}];
        env_cmd = 'module load gurobi/10.0.0 && module load matlab/r2022b && ';
        cmd = ['matlab -nodisplay -nosplash -nojvm -r \"cluster_a2_sensitivity_analysis_sampling ',tmpDir,' ',strs{1},' ',strs{2},'\"'];
        full_cmd = [env_cmd, cmd, ' > ',tmpDir,'/',batchID,'.log && wait'];
        bsub_cmd = ['bsub -q long -W 48:00 -n 1 -R rusage[mem=10000] -e ',tmpDir,'/err_',batchID,'.log ','-J ',batchID,' "'];
        cmd_ready = [bsub_cmd, full_cmd,'"'];
        system(cmd_ready);
        fprintf('batch %s_%s was resubmitted...\n',strs{1},strs{2});
        pause(0.25);
    end
    fprintf('job monitoring round finished; %d batches remain unfinished; checkpoint interval is 5 min...\n',length(runningBatches));
    pause(300);
end

%% clean up 
samples_CSM = {};
samples_PFM_sol = {};
samples_OFM_sol = {};
for i = 1:length(samplingDepth)
    for j = 1:samplingReps
        load([tmpDir,'/iMAT_WPS_subsample_',num2str(i),'_', num2str(j),'.mat'])
        samples_CSM{i,j} = myCSM_ori;
        samples_PFM_sol{i,j} = [minval_PFM', maxval_PFM'];
        samples_OFM_sol{i,j} = [minval_OFM', maxval_OFM'];
    end
end

save('output/sampling_results.mat','samples_CSM','samples_PFM_sol','samples_OFM_sol','samplingDepth','samplingReps');


