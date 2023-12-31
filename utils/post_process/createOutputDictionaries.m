%MAKE COMPLETE COMMENT
%
%create output dictionaries that are friendlier to analysis
%
function[vars_dict,dist_dict,Sinc_dict,separated_dict] = createOutputDictionaries(ground_params,brdf_params,meas_params,msbss_params,...
    results_agnostic,results_precon,results_jade,results_mscpa)

spec_ground = ground_params.spec_ground;
num_ground = size(spec_ground,1);
num_specs = size(spec_ground,2);

do_precon = isfield(results_precon,'Sinc_precon');
do_jade = isfield(results_jade,'Sinc_jade');
do_mscpa = isfield(results_mscpa,'Sinc_mscpa');

if do_jade
    Hjade = results_jade.Hjade;
    Sinc_jade = results_jade.Sinc_jade;
    jade_best_i = results_jade.jade_best_i;
    Sinc_jade_best = results_jade.Sinc_jade_best;
end

if do_precon
    Sinc_precon_up = results_precon.Sinc_precon_up;
    Sinc_noprecon_up = results_precon.Sinc_noprecon_up;
end

%% holds important variables
vars_dict = struct();
vars_dict.inc_angles = brdf_params.brdf_inc_angs_up; 
vars_dict.spec_ground = meas_params.spec_ground; 
vars_dict.spec_clutter = meas_params.spec_clut; 
if ground_params.single_ground
    vars_dict.spec_sing_ground = results_agnostic.spec_sing_ground;
end

%% holds distance/residual information
dist_dict = struct();
dist_dict.dist_agnostic = results_agnostic.dist_agnostic;
% if do_precon 
%     dist_dict.dist_precon = results_precon.dist_precon;
% end
if do_jade 
    dist_dict.dist_jade = results_jade.dist_jade;
end
if do_mscpa 
    dist_dict.dist_mscpa = results_mscpa.dist_mscpa;
end

%% holds reconstruction radiance distribution information
Sinc_dict = struct();
Sinc_dict.each = struct(); %each individual ground
Sinc_dict.spec = struct(); %saved by spectrum
Sinc_dict.total = struct(); %combined grounds

Sinc_dict.each.ground = results_agnostic.Sinc_ground;
Sinc_dict.total.ground = sum_cell(results_agnostic.Sinc_ground,sum(spec_ground,2)); 
Sinc_dict.each.clutter = results_agnostic.Sinc_clut;
Sinc_dict.total.clutter = sum_cell(results_agnostic.Sinc_clut); 
if ground_params.single_ground
    Sinc_dict.each.groundone = results_agnostic.Sinc_sing_ground;
    Sinc_dict.total.groundone = sum_cell(results_agnostic.Sinc_sing_ground);
end

%each individual Sinc
for ground_i=1:num_ground
    for spec_c=1:num_specs
        coeff_multiple = 1;%/(inten_coeff(ground_i))

        Sinc_dict.each.agnostic{ground_i,spec_c} = results_agnostic.Sinc{spec_c}*coeff_multiple;
        Sinc_dict.each.singfilt{ground_i,spec_c} = results_agnostic.Sinc_singfilt{spec_c}*coeff_multiple;
        if do_jade
            if isfield(msbss_params,'bss_ver') && msbss_params.bss_ver==2
                Sinc_dict.each.jade{ground_i,spec_c} = sum_cell(Sinc_jade,-1*(Hjade(:,spec_c)))*coeff_multiple;
            else
                Sinc_dict.each.jade{ground_i,spec_c} = Sinc_jade_best{ground_i}*Hjade(jade_best_i(ground_i),spec_c)*coeff_multiple;
            end
        else
            Sinc_dict.each.jade{ground_i,spec_c} = nan(size(results_agnostic.Sinc{spec_c}));
        end

        if do_precon
            Sinc_dict.each.precon{ground_i,spec_c} = Sinc_precon_up{spec_c}*coeff_multiple;
            Sinc_dict.each.noprecon{ground_i,spec_c} = Sinc_noprecon_up{spec_c}*coeff_multiple;
        end

        if do_mscpa
            Sinc_dict.each.mscpa{ground_i,spec_c} = results_mscpa.Sinc_mscpa_multispec{spec_c}*coeff_multiple;
        else
            Sinc_dict.each.mscpa{ground_i,spec_c} = nan(size(Sinc{spec_c}));
        end
    end
end

%total Sincs
coeff_multiple = 1;
Sinc_dict.total.agnostic = sum_cell(results_agnostic.Sinc)*coeff_multiple;
Sinc_dict.total.singfilt = sum_cell(results_agnostic.Sinc_singfilt);
if do_jade
    if isfield(msbss_params,'bss_ver') && msbss_params.bss_ver==2
        if recon_params.smooth_sincs
            for i=1:length(Sinc_jade)
                f = fit((1:length(Sinc_jade{i}))',Sinc_jade{i}','smoothingspline','SmoothingParam',Sinc_smooth);
                Sinc_jade{i} = f(1:length(Sinc_jade{i}))';
            end
        end
        Sinc_dict.total.jade = sum_cell(Sinc_jade,sum(Hjade,2))*coeff_multiple;   
    else
        Sinc_jade_total = sum_cell(Sinc_jade_best,sum(Hjade(jade_best_i,:),2));
        Sinc_dict.total.jade = Sinc_jade_total*coeff_multiple;   
    end
end
if do_precon
    Sinc_dict.total.precon = sum_cell(Sinc_precon_up)*coeff_multiple;   
    Sinc_dict.total.noprecon = sum_cell(Sinc_noprecon_up)*coeff_multiple; 
    Sinc_dict.total.noprecon(isnan(Sinc_dict.total.noprecon)) = 0;
end
if do_mscpa
    Sinc_dict.total.mscpa = sum_cell(results_mscpa.Sinc_mscpa_multispec)*coeff_multiple;
end


%% holds important info for debugging BSS and Opt Precon
separated_dict = struct();
separated_dict.lfield.ground = meas_params.lfield_ground;
separated_dict.lfield.clutter = meas_params.lfield_clutter;
separated_dict.lfield.meas = meas_params.lfield;
separated_dict.alpha.ground = meas_params.alpha_ground;
separated_dict.alpha.clutter = meas_params.alpha_clutter;
separated_dict.alpha.meas = meas_params.alpha_meas;
separated_dict.dalpha.ground = meas_params.dalpha_ground;
separated_dict.dalpha.clutter = meas_params.dalpha_clutter;
separated_dict.dalpha.meas = meas_params.dalpha_meas;
separated_dict.Sinc_vant.ground = results_agnostic.Sinc_vant_ground;
separated_dict.Sinc_vant.clutter = results_agnostic.Sinc_vant_clut;
separated_dict.Sinc_vant.meas = results_agnostic.Sinc_vant;

if do_jade
    separated_dict.lfield.jade = results_jade.lfield_jade(jade_best_i);
    separated_dict.alpha.jade = results_jade.alpha_jade(jade_best_i);
    separated_dict.dalpha.jade = results_jade.dalpha_jade(jade_best_i);
end
if do_precon
    separated_dict.optprecon.meas = results_precon.Sinc_precon;
    separated_dict.optprecon.ground = results_precon.Sinc_precon_ground;
    separated_dict.optprecon.clutter = results_precon.Sinc_precon_clutter;
end
if ground_params.single_ground
    separated_dict.lfield.groundone = meas_params.lfield_sing_ground;
    separated_dict.alpha.groundone = meas_params.alpha_sing_ground;
    separated_dict.dalpha.groundone = meas_params.dalpha_sing_ground;
    separated_dict.lincdw.groundone = results_agnostic.Sinc_vant_sing_ground;
end
