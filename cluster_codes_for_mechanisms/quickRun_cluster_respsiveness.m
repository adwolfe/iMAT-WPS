addpath integration_pipelines/

rng(1030)

predTbl = readtable('output/fluxTable.csv');
resp_depend = predTbl.rxns(predTbl.PFD_bounded_exp_resp & ~predTbl.PFD_bounded_exp_only);
initCobraToolbox(false);
environment = getEnvironment();
tmpDir = ['tmp_',num2str(100+fix(rand()*100000))];
mkdir(tmpDir);
uid = tmpDir(5:7);
save([tmpDir,'/environment.mat'],'environment','-v7.3','-nocompression');

% make the target rxn pools
splits = 1:length(resp_depend);
for i = 1:length(splits)
    batchID = ['b_',uid,'_',num2str(i)];
    env_cmd = 'module load gurobi/10.0.0 && module load matlab/r2022b && ';
    cmd = ['matlab -nodisplay -nosplash -nojvm -r \"cluster_a2_prediction_mechanism_analysis_resp ',tmpDir,' ',num2str(splits(i)),'\"'];
    full_cmd = [env_cmd, cmd, ' > ',tmpDir,'/',batchID,'.log && wait'];
    bsub_cmd = ['bsub -q long -W 12:00 -n 1 -R rusage[mem=4096] -e ',tmpDir,'/err_',batchID,'.log ','-J ',batchID,' "'];
    cmd_ready = [bsub_cmd, full_cmd,'"'];
    system(cmd_ready);
    pause(0.1);
end

%% start job monitoring
fprintf('job pooling finished, pausing...\n');
pause(10); % wait for all jobs to be submitted
fprintf('start job monitoring...\n');
runningBatches = {};
for i = 1:length(splits)
    runningBatches{i} = ['b_',uid,'_',num2str(i)];
end

while ~isempty(runningBatches)
    allFiles = dir([tmpDir,'/*.csv']);
    allFiles = {allFiles.name};
    finished = regexprep(allFiles,'^fluxtale_ann_|.csv.?$','');
    runningBatchesID = setdiff(regexprep(runningBatches,'^b_..._',''),finished);
    runningBatches = strcat(['b_',uid,'_'],runningBatchesID);
    % find the failed runs
    [~,out] = system('bjobs -q long');
    allTerms = strsplit(out,' ');
    failedBatches = setdiff(runningBatches,allTerms);
    % resubmit the failed batches 
    failedBatches = cellfun(@str2num, regexprep(failedBatches,'^b_..._',''));
 
    for j = 1:length(failedBatches)
        i = failedBatches(j);
        batchID = ['b_',uid,'_',num2str(i)];
        env_cmd = 'module load gurobi/10.0.0 && module load matlab/r2022b && ';
        cmd = ['matlab -nodisplay -nosplash -nojvm -r \"cluster_a2_prediction_mechanism_analysis_resp ',tmpDir,' ',num2str(splits(i)),'\"'];
        full_cmd = [env_cmd, cmd, ' > ',tmpDir,'/',batchID,'.log && wait'];
        bsub_cmd = ['bsub -q long -W 12:00 -n 1 -R rusage[mem=4096] -e ',tmpDir,'/err_',batchID,'.log ','-J ',batchID,' "'];
        cmd_ready = [bsub_cmd, full_cmd,'"'];
        system(cmd_ready);
        fprintf('batch %d was resubmitted...\n',i);
        pause(0.25);
    end
    fprintf('job monitoring round finished; %d batches remain unfinished; checkpoint interval is 5 min...\n',length(runningBatches));
    pause(300);
end

%% clean up 
predTbl = readtable('output/fluxTable.csv');
resp_depend = predTbl.rxns(predTbl.PFD_bounded_exp_resp & ~predTbl.PFD_bounded_exp_only);
predTbl.related_responsiveness_constraints = repmat({'ND'}, size(predTbl,1),1);
predTbl.related_similarity_constraints = repmat({'ND'}, size(predTbl,1),1);

for i = 1:length(splits)
    tmp = readtable([tmpDir,'/fluxtale_ann_',num2str(i),'.csv']);
    predTbl(strcmp(predTbl.rxns, resp_depend{i}),:) = tmp(strcmp(tmp.rxns, resp_depend{i}),:);
end

writetable(predTbl,'output/mechanism_analysis/fluxTable_ann_responsiveness.csv')



