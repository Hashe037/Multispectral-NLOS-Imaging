%NEED TO RECOMMENT
%
%flip W and H and find best Linc with respect to ground
%allows for matching matrices to multiple grounds and therefore make
%comparison useful in the end.
%
%if not enough components for all the ground, make the rest zero vectors
%
%inputs
%W,H--W,H matrix for BSS method in question 
%lfield--light fields (just in case we need to flip)
%Sinc--cell of comp reconstructions where length(Sinc) = # components
%Sinc_ground--cell of ground reconstructions where length(Sinc_ground) = # sources
%doflip--whether to actually flip or not (don't want to flip BSS on light
%field just on dalpha or zero-mean things)

%outputs
%W2,H2--Flipped W,H
%Sinc_best--cell length # sources of each best component reconstruction.
%In order
%best_i--list of component index which is the best fit for the given ground
%index, so best_i(ground_i) = index in L_inc that fits ground_i 
function[W,H,lfield,Sinc_best,best_i] = find_flip_best_fit(W,H,lfield,Sinc,Sinc_ground,doflip)

total_spots = length(Sinc{1}); %total spots
taken = []; %don't reuse components that were taken 

%find all mses
mses = {}; %holds all the calculated mses
for ground_i = 1:length(Sinc_ground) %each ground
    L_inc_g = Sinc_ground{ground_i};
    for comp=1:length(Sinc) %each component
        Sinc2 = sum(H(comp,:))*Sinc{comp};
        if norm(Sinc2) > 0 %not a null/zero vector
            mses{ground_i}(comp,1) = nansum((-Sinc2/norm(Sinc2)+L_inc_g/norm(L_inc_g)).^2)/total_spots; %not flipped
            mses{ground_i}(comp,2) = nansum((Sinc2/norm(Sinc2)+L_inc_g/norm(L_inc_g)).^2)/total_spots; %flipped
        else %remove zero vector
            mses{ground_i}(comp,1) = 1e10;
            mses{ground_i}(comp,2) = 1e11;
        end
    end
end

%find minimums for each ground
ground_mins = zeros(length(Sinc_ground),1);
for ground_i = 1:length(Sinc_ground)
    [~,ground_mins(ground_i)] = min(mses{ground_i}(:));
end
    
%rearrange to go over ground_i with minimums
[~,correct_indices] = sort(ground_mins,'ascend');
    
%find best fitting for each ground
for ii = 1:min(length(Sinc_ground),length(Sinc)) %min between ground and components
    ground_i = correct_indices(ii);
    notfound = 1; %still search
    noflip = 1;
    while notfound %havent found non-taken minimum yet
        [~,mini] = min(mses{ground_i}(:));
        [minr,minc] = ind2sub(size(mses{ground_i}),mini);
%         minr = mini(1);
%         minc = mini(2);
        if any(minr==taken) == 0 %not taken yet
            notfound = 0;
            best_i(ground_i) = minr;
            taken(end+1) = minr;
            if ~doflip || (minc==1 && noflip)% && sign(sum(L_inc{minr})) == sign(sum(L_inc_gs{ground_i})) %don't flip
                Sinc_best{ground_i} = Sinc{minr};
            else %flip
                Sinc_best{ground_i} = -Sinc{minr};
                lfield{ground_i} = -lfield{ground_i};
%                 W(:,best_i(ground_i)) = -W(:,best_i(ground_i));
%                 H(best_i(ground_i),:) = H(best_i(ground_i),:); %dont flip both? FLIP IF LS GROUND
            end
        else %has been taken so make that value high
            noflip = (minc==1); %did want to flip though so keep track of that
            mses{ground_i}(minr,1) = realmax;
            mses{ground_i}(minr,2) = realmax;
        end
    end
end
            

