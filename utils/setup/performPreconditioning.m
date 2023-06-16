function[meas_params,dalpha_occ,diff_occ] = performPreconditioning(scene_params,brdf_params,meas_params,ground_params,recon_params)

brdf = brdf_params.brdf;
V_tot_up = brdf_params.V_tot_up;
U_vant_t = brdf_params.U_vant_t;
D_tot = brdf_params.D_tot;
brdf_inc_angs_up = brdf_params.brdf_inc_angs_up;
brdf_inc_angs = brdf_params.brdf_inc_angs;
brdf_scat_angs = brdf_params.brdf_scat_angs;
vant_pos = scene_params.vant_pos;
vant_angs = scene_params.vant_angs;
cam_angs = scene_params.cam_angs;
occ_x = scene_params.occ_x;
occ_y = scene_params.occ_y;
lfield = meas_params.lfield;
lfield_ground = meas_params.lfield_ground;
lfield_clutter = meas_params.lfield_clutter;

%expected scattering structure from occluder position (s^scat_* in paper)
%size num_positions x vector_num
dalpha_occ = findDalphaOcc(V_tot_up,vant_pos,...
    brdf_inc_angs_up,occ_x,occ_y); 
diff_occ = findDiffOcc(brdf,vant_pos,brdf_inc_angs,...
    vant_angs,cam_angs,brdf_scat_angs,...
    occ_x,occ_y);
diff_occ = diff_occ(:,1:(end-1));

%preconditioned measurements of the rest
[alpha_meas,dalpha_meas,nonsmooth_alpha_meas,nonsmooth_dalpha_meas] = findAlphaDalpha(lfield,U_vant_t,D_tot,recon_params.alpha_smooth);
[alpha_ground,dalpha_ground,~,~] = findAlphaDalpha(lfield_ground,U_vant_t,D_tot,recon_params.alpha_smooth);
[alpha_clutter,dalpha_clutter,~,~] = findAlphaDalpha(lfield_clutter,U_vant_t,D_tot,recon_params.alpha_smooth);
if ground_params.single_ground
    [alpha_sing_ground,dalpha_sing_ground,~,~] = findAlphaDalpha(meas_params.lfield_sing_ground,U_vant_t,D_tot,recon_params.alpha_smooth);
end

%differential measurements
diff_meas = {};
for spec = 1:length(lfield)
    temp = smoothdata(diff(lfield{spec},1),1);
    diff_meas{spec} = temp;
end

%save to meas_params structure
meas_params.alpha_meas = alpha_meas;
meas_params.alpha_ground = alpha_ground;
meas_params.alpha_clutter = alpha_clutter;
meas_params.dalpha_meas = dalpha_meas;
meas_params.dalpha_ground = dalpha_ground;
meas_params.dalpha_clutter = dalpha_clutter;
meas_params.nonsmooth_alpha_meas = nonsmooth_alpha_meas;
meas_params.nonsmooth_dalpha_meas = nonsmooth_dalpha_meas;
meas_params.diff_meas = diff_meas;
% meas_params.spec_ground = spec_ground;
% meas_params.spec_clutter = spec_clutter;
if ground_params.single_ground
    % meas_params.lfield_sing_ground = lfield_sing_ground;
    % meas_params.spec_sing_ground = spec_sing_ground;
    meas_params.alpha_sing_ground = alpha_sing_ground;
    meas_params.dalpha_sing_ground = dalpha_sing_ground;
end
