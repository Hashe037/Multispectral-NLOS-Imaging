%Flatten a 1D cell into 2D matrix where first indice is cell index
%Assumption is that the max dimension of each cell array is also 1D

function[arr] = flatten_cell(cell)

arr = [];
for i=1:length(cell)
    arr(i,:) = cell{i}(:);
end