%account for spectral responsivity of the camera by weighting each light
%field and spectra measurement according to a predefined filt_response

function[lfield,obj_spec,clut_spec] ...
    = accountForResponsivity(lfield,obj_spec,clut_spec,scene_params)

num_clutter = size(clut_spec,1);
num_ground = size(obj_spec,1);
num_specs = length(lfield);

load(scene_params.responsivity_location) %loads in filt_response
for spec_c = 1:num_specs
    lfield{spec_c} = lfield{spec_c}/filt_response(spec_c);

    for ground_i = 1:num_ground
        obj_spec(ground_i,spec_c) = obj_spec(ground_i,spec_c)/filt_response(spec_c);
    end
    for clutter_i = 1:num_clutter
        clut_spec(clutter_i,spec_c) = clut_spec(clutter_i,spec_c)/filt_response(spec_c);
    end
end
