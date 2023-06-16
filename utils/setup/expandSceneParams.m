%expand scene_params and save variables to workspace
function[num_specs,back_sig,cam_angs_all,cam_fov,cam_dist,ang_cut,vant_cut,vant_pos_all,occ_x,occ_y] ...
    = expandSceneParams(scene_params)

num_specs = scene_params.num_specs; %number of spectra measurements
back_sig = scene_params.back_sig; %background signal (noise floor) for light fields
cam_angs_all = scene_params.cam_angs; %deg, camera angles for each picture
cam_fov = scene_params.cam_fov; %cm, field-of-view for the camera
cam_dist = scene_params.cam_dist; %cm, distance between camera and samples
ang_cut = scene_params.ang_cut; %indices, range of camera angles to consider
vant_cut = scene_params.vant_cut; % indices, range of samples vantages to consider
vant_pos_all = scene_params.vant_pos_all; %cm, xi positionss (before being cut)
occ_x = scene_params.occ_x; %cm, occluder horizontal position from center of sample
occ_y = scene_params.occ_y; %cm, occluder vertical position from occluder