%Find SVD values for scattering angles in BRDF
%Define # svd vectors ahead of time
%
% U is left singular vectors (scattering response)
% V is right singular vectors (incident angles)
%
function[U_tot,D_tot,V_tot,U_vant,D_vant,V_vant,vector_num]  ...
    = findBrdfComponents(brdf, vant_angs, brdf_scat_angles, brdf_params)

vector_num = brdf_params.vector_num; %how many SVD vectors to keep
num_vants = size(vant_angs,1); %amount of vantage positions

U_tot = [];
D_tot = [];
V_tot = [];
U_vant = {};
D_vant = {};
V_vant = {};

%construct the svd of all brdf
[Us,Ds,Vs] = svd(brdf,'econ');
Ds = diag(Ds);
V_tot = Vs(:,1:vector_num);
U_tot = Us(:,1:vector_num);
D_tot = Ds(1:vector_num);

%truncate svd for each vantage based on that vantage's angular range
for i=1:num_vants 
    brdf_ang_inds = [];
    for j = 1:length(vant_angs(i,:)) %each scattering angle for the vantage
       [~,brdf_ang_inds(j)] = min(abs(vant_angs(i,j)-brdf_scat_angles));
    end
    U_vant{i} = U_tot(brdf_ang_inds,:);
end

%remove vectors (if needed)
if ~isempty(brdf_params.vectors_remove)
    [U_tot,V_tot,D_tot,U_vant,~,vector_num] ...
       = remove_brdf_vectors(brdf_params.vectors_remove,U_tot,V_tot,D_tot,U_vant,V_vant);
end