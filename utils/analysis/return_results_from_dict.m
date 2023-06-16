%return certain fields of structure across all cell indices. Designed to use
%in displaying results of multiple runs

function[arr] = return_results_from_dict(dict_cell,fields)

numfields = length(fields);

arr = [];
for i=1:length(dict_cell)
    if numfields == 1
        results = getfield(dict_cell{i},fields);
    elseif numfields == 2
        results = getfield(dict_cell{i},fields(1),fields(2));
    elseif numfields == 3
        results = getfield(dict_cell{i},fields(1),fields(2),fields(3));
    elseif numfields == 4
        results = getfield(dict_cell{i},fields(1),fields(2),fields(3),fields(4));
    elseif numfields == 5
        results = getfield(dict_cell{i},fields(1),fields(2),fields(3),fields(4),fields(5));
    else
        error('Error: too many fields (not supported) \n')
    end
    
    %check dimension of results
    if iscell(results)
        arr(i,:,:) = flatten_cell(results);
    elseif sum(size(squeeze(results))==1)==0 %is an 2d array
        arr(i,:,:) = results;
    elseif sum(size(squeeze(results))==1)==1 %is a 1d array
        arr(i,:) = results;
    else %scalar
        arr(i) = results;
    end
end

