%construct 3-D forward model for the computational periscopy method
%essentially just BRDF with cut across from edge but have to worry about
%information across space AND angle
%
% Outputs ------------------------------------------------------------------
% forward_mod_true -- not carrying about how each vantage only views a subset of angles (num_vants x upsampled_cam_angs x num_inc_angs)
% forward_mod_trunc -- truncated to correct outgoing angles based on scattering position (num_vants x cam_angs x num_inc_angs)
%
function [forward_mod_true,forward_mod_trunc] = construct_forward_mod(brdf,vant_pos,brdf_inc_angs,vant_angs,cam_angs,occ_x,occ_y,brdf_scat_angs)

forward_mod_true = nan(length(vant_pos),size(brdf,1),length(brdf_inc_angs));
forward_mod_trunc = nan(length(vant_pos),length(cam_angs),length(brdf_inc_angs));
% brdf_scat_angles = linspace(total_min,total_max,reso); %BRDF scattering angles
for i=1:(length(vant_pos)) %for each vantage except last one
    
    %find cutoff angle
    x_1 = vant_pos(i);
    theta_1 = atand((x_1+occ_x)/occ_y);
    abo = find(brdf_inc_angs>(-1*theta_1)); %find all vant_angles indices above angles 

    %find closest scattering angles and popular the forward model
    brdf_ang_inds = [];
    for j = 1:length(vant_angs(i,:)) %each scattering angle
       [~,brdf_ang_inds(j)] = min(abs(vant_angs(i,j)-brdf_scat_angs));
    end
    forward_mod_true(i,brdf_ang_inds,abo) = brdf(brdf_ang_inds,abo);% value is the BRDF all above certain incidence and over outgoing angle range
    
    %interpolate in between values for better numerical results
    for j = 1:length(abo)
        forward_mod_true(i,:,j) = interp1(brdf_ang_inds',squeeze(forward_mod_true(i,brdf_ang_inds,j)),1:size(forward_mod_true,2) ...
            ,'spline',0);%interpolate in between
    end
    forward_mod_trunc(i,:,:) = forward_mod_true(i,round(linspace(1,size(brdf,1),length(cam_angs))),:);
end

%set all nan values to minimum of nonnan
forward_mod_true(isnan(forward_mod_true)) = min(forward_mod_true(:));
forward_mod_trunc(isnan(forward_mod_trunc)) = min(forward_mod_trunc(:));

%debug
debug = 0;
if debug
    figure,imagesc(squeeze(forward_mod_true(50,:,:)))
    xlabel('Incident Angle'),ylabel('Outgoing Angle')
    
    figure,imagesc(squeeze(forward_mod_trunc(200,:,:)))
    xlabel('Incident Angle'),ylabel('Outgoing Angle')
    
%     figure,imagesc(squeeze(forward_mod_temp(:,:,:)))
%     xlabel('Incident Angle'),ylabel('Outgoing Angle')
end
