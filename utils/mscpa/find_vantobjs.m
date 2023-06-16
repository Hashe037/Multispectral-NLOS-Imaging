%find the indices of the ground objects for each vantage. This is by the
%maximum strength of each object in that vantage, but we have to account
%for spectra used

function[vant_objs,vant_strengths] = find_vantobjs(dalpha_ground,spec_ground)

num_vants = size(dalpha_ground{1},1); %number of spatial vantages
% num_specs = size(spec_ground,2); %spectra in use

vant_objs = zeros(num_vants,1); %index for each vantage
vant_strengths = zeros(num_vants,1); %portion of power the index have over all grounds

for vant = 1:num_vants
    maxPower = -inf; %current max power
    totalPower = 0; %current total power
    for obj_i = 1:length(dalpha_ground)
        obj_vant = dalpha_ground{obj_i}(vant,:);
        obj_power = sum(norm(obj_vant(:))*spec_ground(obj_i,:).^2);
        totalPower = totalPower + obj_power;
        if obj_power > maxPower
            vant_objs(vant) = obj_i;
            maxPower = obj_power;
        end
        vant_strengths(vant) = maxPower/totalPower;
    end
end


