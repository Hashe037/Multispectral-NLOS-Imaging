
%vant_angs, vant_brdf, U_vant_t, D_tot, brdf_inc_angs_up
function[scene_params,brdf_params] ...
    = defineBrdfVariables(scene_params,brdf_params)

vant_pos = scene_params.vant_pos;
cam_angs = scene_params.cam_angs;
cam_angs_all = scene_params.cam_angs_all;
cam_dist = scene_params.cam_dist;

%loadin brdf
[brdf, brdf_scat_angs, brdf_inc_angs, brdf_scat_range] = loadin_brdf(brdf_params);

%seperate vantages into their own angular extent
[vant_angs,vant_angs_all,vant_brdf,vant_brdf_all] ...
    = findVantAngs(brdf,vant_pos,cam_angs,cam_angs_all,brdf_scat_angs,cam_dist);

% construct SVD for overall and truncate for each vant_brdf
%U is LEFT sing vector matrix which is numScatAngles x numScatAngles
%V is RIGHT sing vector matrix which is numIncAngles x numIncAngles
[U_tot,D_tot,V_tot,U_vant,D_vant,V_vant,brdf_params.vector_num]  ...
    = findBrdfComponents(brdf,vant_angs,brdf_scat_angs,brdf_params);

[U_vant_norms,U_vant_t,U_tot_t,brdf_up,V_tot_up,brdf_inc_angs_up] ...
    = processBrdfComponents(U_vant,V_tot,brdf,brdf_inc_angs,length(cam_angs));

scene_params.vant_angs = vant_angs;
scene_params.vant_brdf = vant_brdf;
brdf_params.brdf = brdf;
brdf_params.brdf_scat_angs = brdf_scat_angs;
brdf_params.brdf_inc_angs = brdf_inc_angs;
brdf_params.brdf_scat_range = brdf_scat_range;
brdf_params.D_tot = D_tot;
brdf_params.U_vant_t = U_vant_t;
brdf_params.V_tot_up = V_tot_up;
brdf_params.brdf_inc_angs_up = brdf_inc_angs_up;
