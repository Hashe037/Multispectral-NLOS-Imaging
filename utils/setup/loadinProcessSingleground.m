
function[lfield_sing_ground] = loadinProcessSingleground(ground_params,meas_params,scene_params)

lfield_sing_ground = {};
for spec=1:scene_params.num_specs
    lfield_sing_ground{spec} = csvread(strcat(ground_params.single_ground_location,"\light_fields\lfield_",num2str(spec),".csv"))-scene_params.back_sig;
end   

%account for spectral responsivity
if scene_params.account_responsivity
    [lfield_sing_ground,~,~] ...
        = accountForResponsivity(lfield_sing_ground,meas_params.spec_ground,meas_params.spec_clutter,scene_params);
end

%cut the light field as same as before
[~,lfield_sing_ground,~,~] ...
    = cutCamAngles(scene_params.cam_angs_all,scene_params.ang_cut,lfield_sing_ground,lfield_sing_ground,lfield_sing_ground);
[~,~,lfield_sing_ground,~] ...
    = cutVantPos(scene_params.vant_pos_all,scene_params.vant_cut,lfield_sing_ground,lfield_sing_ground,lfield_sing_ground);

