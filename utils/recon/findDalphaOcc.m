%Find the expected scattering structure (after preconditing) given the
%occluder position and the BRDF properties. For an edge-occluder, this is
%found by integrating over the angles non-obstructed for each scattering
%vantage. Upsampled incident angles is ideal to help with the integration.
%
% Output: dalpha is size num_vants x num_vectors

function[dalpha] = findDalphaOcc(V,vant_pos,angs,occ_x,occ_y)

const_factor = 1e5; %arbitrary scaling to affect the end results scale
occ_x = -occ_x; %make occ_x negative

dalpha = zeros(size(V,2),length(vant_pos)-1); %size of singular vectors and vantages
for i=1:(length(vant_pos)-1) %for each vantage except last one
    x_1 = vant_pos(i); %vant 1
    x_2 = vant_pos(i+1); %vant 2
    theta_1 = atand((occ_x-x_1)/occ_y); %cutoff angle for vant 1
    theta_2 = atand((occ_x-x_2)/occ_y); %cutoff angle for vant 2 (should be larger than theta_1)
    
    abo = find(angs<theta_1); %find all angles above angle 1 (larger/more glancing angles)
    belo = find(angs>theta_2); %find all angles below angle 2 (smaller angles)

    delta_ang = intersect(abo,belo); %find INDICES of angles that correspond to the difference in incident light received
    dalpha(:,i) = sum(V(delta_ang,:),1)'/(length(delta_ang))*abs(angs(2)-angs(1))*const_factor; %integrate over singular vectors
    
    %account for not high enough sampling so the beta is nan
    if sum(isnan(dalpha(:,i))) > 0
        if i>1
            dalpha(:,i) = dalpha(:,i-1);
        end
    end    
end
dalpha = dalpha';
    
    
    
    