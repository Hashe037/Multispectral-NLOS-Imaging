%% instantiate some important global variables
%import variables
[num_specs,back_sig,cam_angs_all,cam_fov,cam_dist,ang_cut,vant_cut,vant_pos_all,occ_x,occ_y] ...
    = expandSceneParams(scene_params);
[brdf,brdf_inc_angs,brdf_scat_angs,brdf_scat_range,~,~] ...
    = expandBrdfParams(brdf_params);

%load light fields and ground truth into workspace
%lfields are size num_vants_all x num_angs_all
[lfield, lfield_ground, lfield_clutter, spec_ground, spec_clutter, num_ground, num_clutter] ...
    = loadinLfields(ground_params,lfield_locations,back_sig,num_specs);

fprintf('Light fields and parameters loaded in! \n')

%% account for camera spectral responsivity
%needed if doing real spectral filters, not needed for OLED monitor
%experiments
if scene_params.account_responsivity
    [lfield,spec_ground,spec_clutter] ...
        = accountForResponsivity(lfield,spec_ground,spec_clutter,scene_params);
    fprintf('Accounting for spectral responsivity \n')
else
    fprintf('Note: not accounting for spectral responsivity \n')
end

    

%% cut out unwanted angles and spatial vantages
%changes size of light fields as well
[cam_angs,lfield,lfield_ground,lfield_clutter] ...
    = cutCamAngles(cam_angs_all,ang_cut,lfield,lfield_ground,lfield_clutter);
num_angs = length(cam_angs); %number of total camera angles

[vant_pos,lfield,lfield_ground,lfield_clutter] ...
    = cutVantPos(vant_pos_all,vant_cut,lfield,lfield_ground,lfield_clutter);
num_vants = length(vant_pos); %number of spatial vantages


%% whether to do single non-uniformly-colored ground
%useful for real-life experiments since this allows a single ground source
%but its spectrum can change vantage-by-vantage (non-uniformly colored)
%repeat many of the previous preprocessing
if ground_params.single_ground
    fprintf('Using single ground parameters \n')
    for spec=1:num_specs
        lfield_sing_ground{spec} = csvread(strcat(ground_params.single_ground_location,"\light_fields\lfield_",num2str(spec),".csv"))-back;
    end   

    %account for spectral responsivity
    if scene_params.account_responsivity
        [lfield_sing_ground,~,~] ...
            = accountForResponsivity(lfield_sing_ground,spec_ground,spec_clutter,scene_params);
    end

    %cut the light field as same as before
    [~,lfield_sing_ground,~,~] ...
        = cutCamAngles(cam_angs_all,ang_cut,lfield_sing_ground,lfield_sing_ground,lfield_sing_ground);
    [~,~,lfield_sing_ground,~] ...
        = cutVantPos(vant_pos_all,vant_cut,lfield_sing_ground,lfield_sing_ground,lfield_sing_ground);
else %dummy variable
    lfield_sing_ground = {zeros(5)};
end


%% seperate vantages into their own angular extent
[vant_angs,vant_angs_all,vant_brdf,vant_brdf_all] ...
    = findVantAngs(brdf,vant_pos,cam_angs,cam_angs_all,brdf_scat_angs,cam_dist);

%% construct SVD for overall and truncate for each vant_brdf
%U is LEFT sing vector matrix which is numScatAngles x numScatAngles
%V is RIGHT sing vector matrix which is numIncAngles x numIncAngles
[U_tot,D_tot,V_tot,U_vant,D_vant,V_vant,brdf_params.vector_num]  ...
    = findBrdfComponents(brdf,vant_angs,brdf_scat_angs,brdf_params);

[U_vant_norms,U_vant_t,U_tot_t,brdf_up,V_tot_up,brdf_inc_angs_up] ...
    = processBrdfComponents(U_vant,V_tot,brdf,brdf_inc_angs,num_angs);


%% perform preconditioning (Section 4.2 of paper)
% Projecting light fields onto the BRDF singular vectors and performing differential
% Also detailed in "Passive non-line-of-sight imaging using plenoptic
% information" (Optics Express 2020)

%expected scattering structure from occluder position (s^scat_* in paper)
%size num_positions x vector_num
dalpha_occ = findDalphaOcc(V_tot_up,vant_pos,brdf_inc_angs_up,occ_x,occ_y); 
diff_occ = findDiffOcc(brdf,vant_pos,brdf_inc_angs,vant_angs,cam_angs,brdf_scat_angs,occ_x,occ_y);
diff_occ = diff_occ(:,1:(end-1));

%preconditioned measurements of the rest
[alpha_meas,dalpha_meas,nonsmooth_alpha_meas,nonsmooth_dalpha_meas] = findAlphaDalpha(lfield,U_vant_t,D_tot,recon_params.alpha_smooth);
[alpha_ground,dalpha_ground,~,~] = findAlphaDalpha(lfield_ground,U_vant_t,D_tot,recon_params.alpha_smooth);
[alpha_clutter,dalpha_clutter,~,~] = findAlphaDalpha(lfield_clutter,U_vant_t,D_tot,recon_params.alpha_smooth);

if ground_params.single_ground
    [alpha_sing_ground,dalpha_sing_ground,~,~] = findAlphaDalpha(lfield_sing_ground,U_vant_t,D_tot,recon_params.alpha_smooth);
end

%differential measurements
diff_meas = {};
for spec = 1:num_specs
    temp = smoothdata(diff(lfield{spec},1),1);
    diff_meas{spec} = temp;
end

%% calculate predictive metrics to estimate the difficulty of NLOS imaging problem
[predmet_dict] = ...
    findPredictiveMetrics(lfield,lfield_ground,lfield_clutter,dalpha_ground,dalpha_clutter,spec_ground,spec_clutter);

%% add useful parameters onto existing parameter structures
%MAKE THIS INTO OUTPUT OF PREVIOUS CODE SO FAR
brdf_params.brdf_inc_angs = brdf_inc_angs;
brdf_params.brdf_inc_angs_up = brdf_inc_angs_up;
brdf_params.brdf_inc_angs = brdf_inc_angs;
brdf_params.brdf_scat_angs = brdf_scat_angs;
brdf_params.brdf = brdf;
brdf_params.D_tot = D_tot;
brdf_params.U_vant_t = U_vant_t;
scene_params.vant_angs = vant_angs;
scene_params.vant_brdf = vant_brdf;
scene_params.vant_pos = vant_pos;
scene_params.cam_angs = cam_angs;
recon_params.dalpha_occ = dalpha_occ;
recon_params.diff_occ = diff_occ;
meas_params.lfield = lfield;
meas_params.lfield_ground = lfield_ground;
meas_params.lfield_clutter = lfield_clutter;
meas_params.alpha_meas = alpha_meas;
meas_params.alpha_ground = alpha_ground;
meas_params.alpha_clutter = alpha_clutter;
meas_params.dalpha_meas = dalpha_meas;
meas_params.dalpha_ground = dalpha_ground;
meas_params.dalpha_clutter = dalpha_clutter;
meas_params.nonsmooth_alpha_meas = nonsmooth_alpha_meas;
meas_params.nonsmooth_dalpha_meas = nonsmooth_dalpha_meas;
meas_params.lfield_sing_ground = lfield_sing_ground;
meas_params.diff_meas = diff_meas;
meas_params.spec_ground = spec_ground;
meas_params.spec_clutter = spec_clutter;
if ground_params.single_ground
    meas_params.spec_sing_ground = spec_sing_ground;
    meas_params.alpha_sing_ground = alpha_sing_ground;
    meas_params.dalpha_sing_ground = dalpha_sing_ground;
end

%% mscpa
% 
% [U_tot,D_tot,V_tot,U_vant,D_vant,V_vant,mscpa_params.mscpa_vector_num]  ...
%     = findBrdfComponents(brdf,vant_angs,brdf_scat_angs,brdf_params);
% 
% [U_vant_norms,U_vant_t,U_tot_t,brdf_up,V_tot_up,brdf_inc_angs_up] ...
%     = processBrdfComponents(U_vant,V_tot,brdf,brdf_inc_angs,num_angs);
% 
% %preconditioned measurements of the rest
% [alpha_meas_mscpa,dalpha_meas_mscpa,nonsmooth_alpha_meas_mscpa,nonsmooth_dalpha_meas_mscpa] = findAlphaDalpha(lfield,U_vant_t,D_tot,mscpa_params.alpha_smooth);
% [alpha_ground_mscpa,dalpha_ground_mscpa,~,~] = findAlphaDalpha(lfield_ground,U_vant_t,D_tot,mscpa_params.alpha_smooth);
% [alpha_clutter_mscpa,dalpha_clutter_mscpa,~,~] = findAlphaDalpha(lfield_clutter,U_vant_t,D_tot,mscpa_params.alpha_smooth);
% 
% if ground_params.single_ground
%     [alpha_sing_ground,dalpha_sing_ground,~,~] = findAlphaDalpha(lfield_sing_ground,U_vant_t,D_tot,mscpa_params.alpha_smooth);
% end
% dalpha_occ_mscpa = findDalphaOcc(V_tot_up,vant_pos,brdf_inc_angs_up,occ_x,occ_y); 
