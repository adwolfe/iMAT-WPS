function const = printConst(MILP,ind, variables)
for ii = 1:length(ind)
    i = ind(ii);
    ids = MILP.A(i,:) ~= 0;
    stoich = MILP.A(i,ids);
    relatedVariables = variables(ids);
    txt = cellfun(@(x, y) [x '*' y], cellstr(num2str(stoich')), relatedVariables, 'UniformOutput', false);
    if strcmp(MILP.csense(i), 'G')
        sign = '>=';
    elseif strcmp(MILP.csense(i), 'L')
        sign = '<=';
    elseif strcmp(MILP.csense(i), 'E')
        sign = '=';
    else
        error('?');
    end

    const(ii,1) = {[strjoin(txt,' + '), ' ', sign, ' ', num2str(MILP.b(i))]};
end
