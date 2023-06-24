function[scene_params, brdf_params, recon_params, msbss_params, mscpa_params, lsrecon_params] ...
    = get_params_real_singleobj_satin(basefolder)

scene_params = struct(); %holds information about the overall setup
scene_params.cam_angs_all = 0:1:70; %camera angles used in light field
scene_params.ang_cut = [2,71]; %2 %min angle index to analyze
scene_params.vant_cut = [10,300];%10; %min angle index to analyze
scene_params.num_specs = 5; %number of spectrums being used
scene_params.back_sig = 184; %noise floor to subtract
scene_params.cam_dist = 10.25*2.54; %cm, distance from sample to lens
scene_params.cam_fov = 9.7; %cm, fov of camera
scene_params.num_vants = 300; %number of spatial vantages
scene_params.vant_pos_all = ...
    (linspace(1024-950,1024+950,scene_params.num_vants)-1024)/2048*scene_params.cam_fov; %cm, vant positions
scene_params.occ_x = 10.6;%(6.6+1.4); %(8.5+1.4)(; %cm; horizontal distance from center of scatterrer to wall
scene_params.occ_y = 11.1;%6.92; %7.7  cm; vertical distance from center of scatterrer to wall
scene_params.account_responsivity = 1;
scene_params.responsivity_location = strcat(basefolder, '\calibration_and_data\filt_sensitivity.mat');

brdf_params = struct(); %holds information about the BRDF and light field
brdf_params.brdf_location = strcat(basefolder, '\calibration_and_data\brdf.mat');
brdf_params.vector_num = 3; %number of vectors to use for SVD
brdf_params.vectors_remove = []; %vectors to remove because troublesome

recon_params = struct(); %holds information about the DFOV algorithm
recon_params.dis_method = "residnorm_totweighted"; %metric to use for DFOV
recon_params.att_dist = 1; %distance to start attenuation results linearly (1 means no attenuation)
recon_params.smooth_Sincs = 1; %smooth the resulting lincs or not
recon_params.Sinc_smooth_param = .0025; %smooth the resulting lincs or not
recon_params.abs_Sincs = 0; %take absolute value of lincs or not
recon_params.alpha_smooth = .1;%.1; %parameter to smooth dalpha
recon_params.do_optprecon = 1;

msbss_params = struct(); %holds information about the BSS algorithms
msbss_params.numcomp = 3; %number of components for BSS
msbss_params.num_clutter = 2; %how many clutter objects to remove (loops this many times)
msbss_params.precon = 'dalpha_inverted'; %what preconditioning to use

mscpa_params = struct(); %holds information about the MS-DFOV algorithm
mscpa_params.K = 2; %expected number of clutter elements to decide nullspace size
mscpa_params.dis_method = 'residnorm_reconweighted'; %which distance metric to use
mscpa_params.vant_smooth = .0025; %.001 %.01 %smoothing on Linc
mscpa_params.attenuate_ondist = 1; %whether to attenuate Linc w.r.t. distance metric
mscpa_params.nullsize_method = 0; %how to determine the nullsize of the CP algorithm
mscpa_params.mscpa_dist_min = .45;
mscpa_params.nullcut = .1; %.01 %if method = 0, the percentage of power to consider
mscpa_params.dalpha_smooth = .1; %.01
mscpa_params.spec_choice = 1; %what spec (if method==0)
mscpa_params.spec_choice_method = 1; %how to choose spec for algorithm (0 for predetermined, 1 for min)
mscpa_params.do_diff = 0; %whether to perform with differential or dalpha (0 for dalpha)

lsrecon_params = struct(); %holds information about the LS reconstruction
lsrecon_params.w = 7.25e9; %1e-3%value to multiply precon
% lsrecon_params.precon_smooth = .025; %smoothing for preconditioning (needs less than spectral methods)