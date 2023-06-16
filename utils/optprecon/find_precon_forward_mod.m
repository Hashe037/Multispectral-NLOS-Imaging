%find the forward model, specifically for the preconditioner least squares

function[A_3d,A_2d] = find_precon_forward_mod(do_diff,scene_params,brdf_params)

%instantiate parameters
brdf = brdf_params.brdf;
vant_pos = scene_params.vant_pos;
brdf_inc_angs = brdf_params.brdf_inc_angs;
vant_angs = scene_params.vant_angs;
brdf_scat_angs = brdf_params.brdf_scat_angs;
cam_angs = scene_params.cam_angs;
occ_x = scene_params.occ_x;
occ_y = scene_params.occ_y;

%construct forward model basic
[~,A_3d] = construct_forward_mod(brdf,vant_pos,brdf_inc_angs,vant_angs,cam_angs,...
    occ_x,occ_y,brdf_scat_angs);

%do differential (if needed)
if do_diff
    forward_mod_diff = zeros(size(A_3d,1)-1,size(A_3d,2),size(A_3d,3));
    for i=1:size(A_3d,2) %go through each camera angle, smooth image and take differential
        im = squeeze(A_3d(:,i,:));
        forward_mod_diff(:,i,:) = diff(im,1);
    end
    A_3d = forward_mod_diff;
end
A_3d(A_3d<1e-4)=0; %set values to 0 to prevent numerical errors
% A_3d(A_3d<1e-3)=0; %set values to 0 to prevent numerical errors

%convert to 2d
A_2d = zeros(size(A_3d,1)*size(A_3d,2),size(A_3d,3));
for i=1:size(A_3d,3)
    A_2d(:,i) = reshape(A_3d(:,:,i),numel(A_2d(:,i)),1);
end

%get rid of columns of forward_mod with mostly 0 values as those
%correspond to elements not in CFOV
for i=1:size(A_2d,2)
    if sum(A_2d(:,i)~=0) < 30 %.5*size(forward_mod,1)
        A_2d(:,i) = 0;
    end
end
