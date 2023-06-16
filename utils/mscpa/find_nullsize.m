%estimate how many clutter elements they are and use that to determine null
%size for MS-DFOV

function[null_size] = find_nullsize(D,ecut)

minD = min(D(:));
maxD = D(1);
totalE = sum(D(:).^2); %total energy
% cutoff = (maxD-minD)*.2+minD;
% cutoff = .01*maxD;
cutoff = ecut*totalE; %cutoff x amount of energy from signal
cumsumD = cumsum(flipud(D).^2);
% null_size = sum(D>sqrt(cutoff));
null_size = length(D)-sum(cumsumD>cutoff);













