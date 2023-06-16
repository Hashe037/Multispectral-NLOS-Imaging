%return certain fields of structure across all cell indices. Designed to use
%in displaying results of multiple runs
%
%Specialized for cells as outcome

function[arr] = return_results_from_dict_cellversion(dict_cell,fields)

numfields = length(fields);

arr = [];
for i=1:length(dict_cell)
    if numfields == 1
        arr(i,:,:) = flatten_cell(getfield(dict_cell{i},fields));
    elseif numfields == 2
        arr(i,:,:) = flatten_cell(getfield(dict_cell{i},fields(1),fields(2)));
    elseif numfields == 3
        arr(i,:,:) = flatten_cell(getfield(dict_cell{i},fields(1),fields(2),fields(3)));
    elseif numfields == 4
        arr(i,:,:) = flatten_cell(getfield(dict_cell{i},fields(1),fields(2),fields(3),fields(4)));
    elseif numfields == 5
        arr(i,:,:) = flatten_cell(getfield(dict_cell{i},fields(1),fields(2),fields(3),fields(4),fields(5)));
    else
        error('Error: too many fields (not supported) \n')
    end
end