function[results_agnostic] = ...
    performAgnosticRecons(scene_params,brdf_params,meas_params,ground_params,recon_params)

dalpha_meas = meas_params.dalpha_meas;
nonsmooth_dalpha_meas = meas_params.nonsmooth_dalpha_meas;
dalpha_ground = meas_params.dalpha_ground;
dalpha_clutter = meas_params.dalpha_clutter;

[Sinc_vant,Sinc,dist_agnostic] = performLfieldRecon(dalpha_meas,scene_params,brdf_params,recon_params);
[nonsmooth_Sinc_vant,nonsmooth_Sinc,~] = performLfieldRecon(nonsmooth_dalpha_meas,scene_params,brdf_params,recon_params);
[Sinc_vant_ground,Sinc_ground,~] = performLfieldRecon(dalpha_ground,scene_params,brdf_params,recon_params);
[Sinc_vant_clut,Sinc_clut,~] = performLfieldRecon(dalpha_clutter,scene_params,brdf_params,recon_params);


%define recon "spectrum"
num_specs = scene_params.num_specs;
spec_agnostic = zeros(1,num_specs);
for i=1:num_specs
    spec_agnostic(i) = norm(Sinc_vant{i});
end
spec_agnostic = spec_agnostic/max(spec_agnostic);

%perform for "single_ground" if needed
if ground_params.single_ground
    dalpha_sing_ground = meas_params.dalpha_sing_ground;
    [Sinc_vant_sing_ground,Sinc_sing_ground,dis_sing_ground] = performLfieldRecon(dalpha_sing_ground,scene_params,brdf_params,recon_params);
    
    %find spectrum (a bit harder since it changes vantage-to-vantage)
    spec_sing_ground = [];
    thresh = .1e-3;
    thresh2 = 1e-3;
    for vant=1:length(Sinc_vant_sing_ground{1})
        for spec=1:num_specs
            if Sinc_vant_sing_ground{1}(vant) > thresh
                spec_sing_ground(vant,spec) = Sinc_vant_sing_ground{spec}(vant)/(max(abs(Sinc_vant_sing_ground{1}(vant)),thresh2));
            else
                spec_sing_ground(vant,spec) = 1;
            end
        end
    end
end


%save results to structure
results_agnostic.Sinc = Sinc;
results_agnostic.Sinc_vant = Sinc_vant;
results_agnostic.nonsmooth_Sinc_vant = nonsmooth_Sinc_vant;
results_agnostic.nonsmooth_Sinc = nonsmooth_Sinc;
results_agnostic.Sinc_vant_ground = Sinc_vant_ground;
results_agnostic.Sinc_ground = Sinc_ground;
results_agnostic.Sinc_vant_clut = Sinc_vant_clut;
results_agnostic.Sinc_clut = Sinc_clut;
results_agnostic.spec_recon = spec_agnostic;
results_agnostic.dist_agnostic = dist_agnostic;
if ground_params.single_ground
    results_agnostic.Sinc_vant_sing_ground = Sinc_vant_sing_ground;
    results_agnostic.Sinc_sing_ground = Sinc_sing_ground;
    results_agnostic.dis_sing_ground = dis_sing_ground;
    results_agnostic.spec_sing_ground = spec_sing_ground;
end