%Find the BRDF and angles specific to each vantage point on the scattering
%surface since each point has slightly different views (parallax). Only
%changes the scattering angles received (we are assuming the full incident
%scene still affects each vantage).
%
%--------------------------------------------------------------------------
% Inputs
% -vant_pos -- position of vantages on scattering surface
% -brdf -- BRDF of surface (num_scat_angs x num_inc_angs)
% -cam_angs/cam_angs_all -- angles of camera with/without cutting angles
% -ang_range -- total scattered range of BRDF
% -total_min -- minimum scattering angle of BRDF
% -total_max -- maximum scattering angle of BRDF
% -L -- distance of camera to center of scattering surface
%
%--------------------------------------------------------------------------
% Outputs
% -vant_angs -- scattering angles of each vantage wrt camera position
% -vant_angs_all -- same except without cutting angles
% -vant_brdf -- BRDF of each vantage wrt actual scattering angles
% -vant_brdf_cut -- same except without ang_cut
%
% DO I CUT BASED ON VANT POS AS WELL?
function[vant_angs,vant_angs_all,vant_brdf,vant_brdf_all] ...
    = findVantAngs(brdf,vant_pos,cam_angs,cam_angs_all,brdf_scat_angles,cam_dist)

%variables
vant_angs = []; %degree, angles measured by camera of vant
vant_brdf = {}; %brdf of vants
num_cam_angs = length(cam_angs); %number of camera angles
num_vants = length(vant_pos); %number of vantages

%find the scattered angles for each vantage
%NOTE: MAKE SURE THIS IS THE SAME AS IN 
for i=1:num_vants
    scat_angles = calc_scat_angs_geometry(vant_pos(i),cam_dist,cam_angs);
    vant_angs(i,:) = scat_angles;
end

%obtain BRDF for each vantage given the scattering angles
for i=1:num_vants
    brdf_ang_inds = [];
    for j = 1:num_cam_angs %each camera angle
       [~,brdf_ang_inds(j)] = min(abs(vant_angs(i,j)-brdf_scat_angles));
    end
    vant_brdf{i} = brdf(brdf_ang_inds,:);
end




%Find same as above except not affected by ang_min/max cut
%variables
vant_angs_all = []; %degree, angles measured by camera of vant
vant_brdf_all = {}; %brdf of vants
num_cam_angs = length(cam_angs_all); %number of camera angles
num_vants = length(vant_pos); %number of vantages

for i=1:num_vants
    scat_angles = calc_scat_angs_geometry(vant_pos(i),cam_dist,cam_angs_all);
    vant_angs_all(i,:) = scat_angles;
end

%obtain BRDF for that point
for i=1:num_vants  
    brdf_ang_inds = [];
    for j = 1:num_cam_angs %each camera angle
       [~,brdf_ang_inds(j)] = min(abs(vant_angs_all(i,j)-brdf_scat_angles));
    end
    vant_brdf_all{i} = brdf(brdf_ang_inds,:);
end
