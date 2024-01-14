function [h1 h2]= plotContribution(model,flux,myMet)
% This is a simple flux tracker to show the flux around a metabolite
% according to a flux distribution
%
% USAGE:
%
%    mytbl = listRxn(model,flux,myMet)
%
% INPUTS:
%    model:             input RECON2.2 model (COBRA model structure)
%    flux:              the flux distribution to inspect
%    myMet:             the metabolite of interest to list flux around
%
% OUTPUT:
%   mytbl:             	a table showing all the flux that produces or
%                       comsumes myMet
%
% `Yilmaz et al. (2020). Final Tittle and journal.
% .. Author: - Xuhang Li, Mar 2020
myrxns = model.rxns(any(model.S(strcmp(model.mets,myMet),:),1));
myInds = any(model.S(strcmp(model.mets,myMet),:),1);
myfluxP = flux(ismember(model.rxns,myrxns)) .*  model.S(strcmp(model.mets,myMet),myInds)';
for i = 1:length(myrxns)
    mytbl(i,1) = myrxns(i);
    mytbl(i,2) = {flux(strcmp(model.rxns,myrxns{i}))};
    mytbl(i,3) = printRxnFormula_XL(model, myrxns{i},0);
    mytbl(i,4) = printRxnFormula(model, myrxns{i},false);
end
mytbl(:,5) = mat2cell(myfluxP / (sum(abs(myfluxP)) /2),ones(length(myfluxP),1),1);

p_tbl = mytbl(cell2mat(mytbl(:,5)) > 0.05,:); % 5% contribution minimal
if sum(cell2mat(p_tbl(:,5))) < 0.99
    p_tbl(end+1,:) = {'others',nan,'','',1-sum(cell2mat(p_tbl(:,5)))};
end
c_tbl = mytbl(cell2mat(mytbl(:,5)) < -0.05,:); % 5% contribution minimal
if sum(cell2mat(c_tbl(:,5))) < 0.99
    c_tbl(end+1,:) = {'others',nan,'','',1-sum(cell2mat(c_tbl(:,5)))};
end


h1 = figure;
names = strcat(p_tbl(:, 1),':',{' '}, p_tbl(:, 4))  ;
wrapped_names = cell(size(names));
% go through each string in the original cell array
for i = 1:numel(names)
    str = names{i};
    if length(str) > 60
        % wrap the string if its length is more than 20
        num_chunks = ceil(length(str) / 60);
        chunk_sizes = repmat(60, 1, num_chunks);
        chunk_sizes(end) = chunk_sizes(end) - sum(chunk_sizes) + length(str);
        chunks = mat2cell(str, 1, chunk_sizes);
        wrapped_names{i} = strjoin(chunks, '\n');
    else
        % otherwise, just copy the string
        wrapped_names{i} = names{i};
    end
end
values = cell2mat(p_tbl(:, end));
% create pie chart
pie(full(values), false(size(values)), p_tbl(:, 1));
legend(wrapped_names, 'Location', 'southoutside');
title(['Contribution of production flux of ',myMet]);

h2 = figure;
names = strcat(c_tbl(:, 1),':',{' '}, c_tbl(:, 4))  ;
wrapped_names = cell(size(names));
% go through each string in the original cell array
for i = 1:numel(names)
    str = names{i};
    if length(str) > 60
        % wrap the string if its length is more than 20
        num_chunks = ceil(length(str) / 60);
        chunk_sizes = repmat(60, 1, num_chunks);
        chunk_sizes(end) = chunk_sizes(end) - sum(chunk_sizes) + length(str);
        chunks = mat2cell(str, 1, chunk_sizes);
        wrapped_names{i} = strjoin(chunks, '\n');
    else
        % otherwise, just copy the string
        wrapped_names{i} = names{i};
    end
end
values = -cell2mat(c_tbl(:, end));
% create pie chart
pie(full(values), false(size(values)), c_tbl(:, 1));
legend(names, 'Location', 'southoutside');
title(['Contribution of consumption flux of ',myMet]);


end