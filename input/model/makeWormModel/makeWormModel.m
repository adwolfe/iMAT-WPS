%% to make the corbra model from an excel table
%% the generic C. elegans model
% initCobraToolbox;
% the <- in rxns table should be change to <==> manually
% ub/lb header to "Lower bound	Upper bound"
% change met header
% add two psudo met to the met table
% the rxns table and met table should be merged into one excel table
%% generic model with uptakes
%% make the dual model with dynamic controls
model = iCEL2model_xl('./iCEL1314_with_uptakes.xlsx');
sideR=0.01; % default is 0.02
sideR_ind=0.001; % default is 0.005
storageR=0.01; % default is 0.01
bacMW=966.28583751;
% dividedBy=100; this is in the stoicheomitry matrix
% add individual met control and scale by the factor
sideRxns =  model.rxns(model.S(strcmp(model.mets,'sideMet[e]'),:)~=0);
sideRxns(strcmp(sideRxns,'EXC9998')) = [];
for i = 1:length(sideRxns)
    model.S(end+1,:) = zeros(1,length(model.rxns));
    model.S(end, strcmp(sideRxns{i},model.rxns)) = model.S(strcmp(model.mets,'sideMet[e]'),strcmp(sideRxns{i},model.rxns)); 
    model.S(end, strcmp('EXC0050',model.rxns)) = sideR_ind*bacMW*0.01;%the 0.01 here is the 1/devidedBy 
    model.csense(end+1) = 'L';
    model.b(end+1) = 0;
    model.BiGG(end+1) = {'NA'};
    model.metCharges(end+1) = 0;
    model.metFormulas(end+1) = {'NA'};
    model.metNames(end+1) = {'Pseudo-metabolite that represents a constraint'};
    model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};
end

% add overall side/storage control
model.S(end+1,:) = zeros(1,length(model.rxns));
model.S(end, strcmp('EXC9999',model.rxns)) = 1; 
model.S(end, strcmp('EXC0050',model.rxns)) = storageR*bacMW*0.01; 
model.csense(end+1) = 'L';
model.b(end+1) = 0;
model.BiGG(end+1) = {'NA'};
model.metCharges(end+1) = 0;
model.metFormulas(end+1) = {'NA'};
model.metNames(end+1) = {'Pseudo-metabolite that represents a constraint'};
model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};

model.S(end+1,:) = zeros(1,length(model.rxns));
model.S(end, strcmp('EXC9998',model.rxns)) = 1; 
model.S(end, strcmp('EXC0050',model.rxns)) = sideR*bacMW*0.01; 
model.csense(end+1) = 'L';
model.b(end+1) = 0;
model.BiGG(end+1) = {'NA'};
model.metCharges(end+1) = 0;
model.metFormulas(end+1) = {'NA'};
model.metNames(end+1) = {'Pseudo-metabolite that represents a constraint'};
model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};

model = changeObjective(model,'BIO0010');
opt = optimizeCbModel(model,'max')
model = changeObjective(model,'RCC0005');
optimizeCbModel(model,'max')
% save the model
model = changeObjective(model,'BIO0010');
model.description = 'iCEL1314 - with uptakes';
save('iCEL1314_withUptakes.mat','model');