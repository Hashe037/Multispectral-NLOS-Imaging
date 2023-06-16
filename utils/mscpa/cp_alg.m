% FINISH COMMENTING
%
%Algorithm that finds the closest points (cp) between line beta_j and
%hyperplane ax+x0=0 solving the langrangian.
%Can find solution in closed form
%
%--------------------------------------------------------------------------
% Inputs
% beta -- beta vector aka revealed slice's structure (Mx1)
% A -- difference matrix which holds span of clutter elements (Mx(N-1))
% a -- nullspace of idfference matrix (nullsize x M)
% x0 -- single spectra measurement to add to span of clutter (Mx1)
% dis_method -- what distance metric to use
%--------------------------------------------------------------------------
% Outputs
% dis -- resulting distance
% L -- recovered radiance coefficient where closest point is L*beta
% a_vals -- values or "strengths" of clutter elements in nullspace
% otherpoint -- closest point in clutter span matrix
%--------------------------------------------------------------------------
%
%Todo
%-how to add distance/regularizers to solution
%
function[dis,L,a_vals,otherpoint] = cp_alg(dalpha_occ,A,a,x0,dis_method)

a = a; %nullspace should be (N-clutter)xN
% L = 0;%L_guess; %start with L=0
P2 = a'*inv(a*a')*a; %projection operator onto nullspace of A


%performing convex opt
L = (dalpha_occ'*P2*x0)/(dalpha_occ'*P2*dalpha_occ);


%defining distance
resid = P2*(L*dalpha_occ-x0); %residual of the subtraction between the two
otherpoint = L*dalpha_occ-resid; %the point with minimum distance between L*beta
a_vals = otherpoint'/a;

ndot = @(x,y) dot(x,y)/norm(x(:))/norm(y(:));

if strcmp(dis_method,"residnorm") %norm of residual
    dis = norm(resid); %z star p2
elseif strcmp(dis_method,"residnorm_totweighted") %norm of residual weighted by length of both
    dis = norm(resid)/(norm(L*dalpha_occ)+norm(x0));
elseif strcmp(dis_method,"residnorm_measweighted") %norm of residual weighted just by measured dalpha
    dis = norm(resid)/(norm(x0)+.00001);
elseif strcmp(dis_method,"residnorm_reconweighted") %norm of residual weighted by reconstruction values
    dis = norm(resid)/(norm(L*dalpha_occ)+norm(otherpoint)+.00001);
elseif strcmp(dis_method,"cp_angle")
    dis = ndot(L*dalpha_occ,otherpoint);
else
    error('Error: dis_method not recognized')
end

% %find projection values
% vec = zeros(size(dalpha_occ'));
% for i=1:size(A,1) %each vector in a
%     projs2(i) = A(:,i)*(L*dalpha_occ-x0)/(A(:,i)*A(:,i)'); %actual projection
%     vec = vec+projs2(i)*A(i,:);
% end
