%remove the component with the largest distance/residual

function[W2,H2] = remove_worseresidual_component(W,H,dist)

sumdist = [];
if iscell(dist)
    for i=1:length(dist)
        sumdist(i) = sum(dist{i});
    end
else
    sumdist = dist;
end

[~,maxi] = max(sumdist);
W2 = horzcat(W(:,1:maxi-1),W(:,maxi+1:end));
H2 = vertcat(H(1:maxi-1,:),H(maxi+1:end,:));
