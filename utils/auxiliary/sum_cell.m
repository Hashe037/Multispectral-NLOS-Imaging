%sum up elements of a cell (mostly for a total L_inc)
%assumption: all cells have the length array in them

function[tarray] = sum_cell(cell,weights)

if nargin < 2 %no weights
    weights =  ones(size(cell));
end

tarray = zeros(size(cell{1}));
for i=1:length(cell)
   tarray = tarray+weights(i)*cell{i} ;
end
