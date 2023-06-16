function[scene_params,meas_params] = loadinProcessLfields(ground_params,scene_params,lfield_locations)

[meas_params] ...
    = loadinLfields(ground_params,lfield_locations,scene_params.back_sig,scene_params.num_specs);

% lfield = meas_params.lfield;
% lfield_ground = meas_params.lfield_ground;
% lfield_clutter = meas_params.lfield_clutter;
% spec_ground = meas_params.spec_ground;
% spec_clutter = meas_params.spec_clutter;

%% account for camera spectral responsivity
%needed if doing real spectral filters, not needed for OLED monitor
%experiments
if scene_params.account_responsivity
    [meas_params.lfield,meas_params.spec_ground,meas_params.spec_clutter] ...
        = accountForResponsivity(meas_params.lfield, meas_params.spec_ground, meas_params.spec_clut,scene_params);
    fprintf('Accounting for spectral responsivity \n')
else
    fprintf('Note: not accounting for spectral responsivity \n')
end
    

%% cut out unwanted angles and spatial vantages
%changes size of light fields as well
[scene_params.cam_angs, meas_params.lfield, meas_params.lfield_ground, meas_params.lfield_clutter] ...
    = cutCamAngles(scene_params.cam_angs_all, scene_params.ang_cut, meas_params.lfield, meas_params.lfield_ground, meas_params.lfield_clutter);
[scene_params.vant_pos, meas_params.lfield, meas_params.lfield_ground, meas_params.lfield_clutter] ...
    = cutVantPos(scene_params.vant_pos_all, scene_params.vant_cut, meas_params.lfield, meas_params.lfield_ground, meas_params.lfield_clutter);

