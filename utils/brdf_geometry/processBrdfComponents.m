%Process SVD vectors for better numerical performance. Downsampling
%scattering angles to match camera and upsampling incident angles for
%better future numerical processing
%
%--------------------------------------------------------------------------
%Inputs
% -should be explained in main script
%
%--------------------------------------------------------------------------
%Ouputs
% -U_vant_norms -- norms of each U_vant (scattering angles response for
% each vantage). Norms change since each vantage views slightly different
% scattering angles
% -U_vant_t -- truncated U_vant to account for camera angles
% -U_tot_t -- truncated U_tot (no care for vantage) to account for cam
% angles
% -brdf_up -- BRDF with upsampled incident angles to allow for better
% integration in finding beta_inc
% -V_tot_up -- V_tot (incident angles response) upsampled to match BRDF
% -brdf_inc_angs_up -- upsampled BRDF angles to match brdf_up

function[U_vant_norms,U_vant_t,U_tot_t,brdf_up,V_tot_up,brdf_inc_angs_up] ...
    = processBrdfComponents(U_vant,V_tot,brdf,brdf_inc_angs,num_angs)

num_vants = length(U_vant);
vector_num = size(U_vant{1},2);
num_scat_brdf = size(brdf,1);
num_inc_brdf = size(brdf,2);

%truncating left singular vectors but we still want to maintain the unity
%norm during the cutoff
U_vant_norms = {};
for i = 1:num_vants
    for v = 1:vector_num
        U_vant_norms{i}(v) = norm(U_vant{i}(:,v));
    end
end

%% downsample U_vant to match actual number of camera angles rather than BRDF resolution
U_vant_t = {};
U_tot_t = [];

%downsample U_vant(left singular vectors) to match camera angles
for i=1:num_vants
    samp_points = linspace(1,size(U_vant{i},1),num_angs);
    for j=1:vector_num
        U_vant_t{i}(:,j) = interp1(1:size(U_vant{i},1),U_vant{i}(:,j),samp_points);
        U_vant_t{i}(:,j) = U_vant_t{i}(:,j)/norm(U_vant_t{i}(:,j))*U_vant_norms{i}(j); %norm should be same as U_vant{i}
    end
end


%% upsample BRDF to make slits match camera angles and then upsample V to match\
%Might not be upsampling?
brdf_up = [];
V_tot_up = [];
brdf_inc_angs_up = [];

samp_points = linspace(1,num_inc_brdf,num_scat_brdf);
for i=1:num_scat_brdf %each scattering angle
    brdf_up(i,:) = interp1(1:num_inc_brdf,brdf(i,:),samp_points); %upsample incident
end

samp_int = linspace(1,num_inc_brdf,num_scat_brdf);
for i=1:vector_num
   V_tot_up(:,i) = interp1(1:num_inc_brdf,V_tot(:,i),samp_points);
   V_tot_up(:,i) = V_tot_up(:,i)/norm(V_tot_up(:,i));
end

brdf_inc_angs_up = interp1(1:num_inc_brdf,brdf_inc_angs,samp_points);