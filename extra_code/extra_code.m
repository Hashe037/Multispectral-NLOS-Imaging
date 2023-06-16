%extra_code

%% perform quanitization if wanted to lower dynamic range
%lower resultion
if do_low_res
    fprintf('Lowering Resolution \n')
    [lfield,lfield_ground,lfield_clutter] = ...
    lower_quant_res(lfield,lfield_ground,lfield_clutter,desired_quant_val)
end

%% smooth light fields slightly across vantages (NOT NEEDED)
%be careful: too much smoothing can hurt, but this should diminish
%differences between vantages
perf_smooth = 0;
vant_smooth_param = .1; %.0001
if perf_smooth==1
    fprintf('Smoothing across vantages \n')
    [lfield,lfield_ground] = ...
        smooth_lfield_vantages(lfield,lfield_ground,vant_smooth_param);
end  
% figure,mesh(lfield_ground)

%% account for spatial inhomogeneities
%Note: because we can never see the full reflection across all angles, the
%norm match and amp match will always be off since each vantage is not
%represented equally so comparing between them is also difficult.

obj_noocc = csvread(ground_params.occ_cal);
max_vant = 250;%35; %max vantage to go to
ref_vant = max_vant;

if do_normalize_curves
    [lfield,lfield_ground,curve_multiples] ...
        = normalize_curves_with_reference(obj_noocc,max_vant,ref_vant,do_norm_match,do_amp_match, ...
        lfield,lfield_ground);
end

%% ground/clut all
%summation of lfields across all spectra
l_field_ground_all = sum_cell(lfield_ground); %total ground of all objects with correct intensity
l_field_clutter_all = sum_cell(lfield_clutter); %total clutter

[alpha_ground_all,dalpha_ground_all,~,~] = find_alpha_dalpha({l_field_ground_all},U_vant_t,D_tot,dalpha_smooth);
alpha_ground_all = alpha_ground_all{1}; dalpha_ground_all = dalpha_ground_all{1}; 
[alpha_clutter_all,dalpha_clutter_all,~,~] = find_alpha_dalpha({l_field_clutter_all},U_vant_t,D_tot,dalpha_smooth);
alpha_clutter_all = alpha_clutter_all{1}; dalpha_clutter_all = dalpha_clutter_all{1}; 

[L_inc_g_dw_all,L_inc_g_all,ground_dis] = ...
    perform_dfov(dalpha_ground_all,beta_inc,r_params);
[L_inc_clut_dw_all,L_inc_clut_all,clutter_dis] = ...
    perform_dfov(dalpha_clutter_all,beta_inc,r_params);

%% find "object-filtered" reconstructions
%weight reconstructions based on object spectra
[Sinc_vant_filt,Sinc_filt,lfield_filt,alpha_filt,dalpha_filt] ...
    = find_filtered_lincs(obj_spec,Sinc_vant,Sinc,lfield,alpha_meas,dalpha_meas);

%weight reconstructions based on object spectra but still keep multiple
%spectra reconstructions
[L_inc_filt2_dw,L_inc_filt2,l_field_filt2,alpha_filt2,dalpha_filt2] ...
    = find_filtered_lincs2(obj_spec,Sinc_vant,Sinc,lfield,alpha_meas,dalpha_meas);

Sinc_singfilt = {};
Sinc_vant_singfilt = {};
for spec=1:length(Sinc_vant) %instantiate
    Sinc_singfilt{spec} = zeros(size(Sinc{1}));
    Sinc_vant_singfilt{spec} = zeros(size(Sinc_vant{1}));
end


%% find the best spectra based on true values
if ground_params.single_ground 
    error_list = [];
    for spec=1:length(Sinc_vant)
        error_list(spec) = nansum(abs(Sinc_vant{spec}/norm(Sinc_vant{spec}(:))-Sinc_vant_sing_ground{1}/norm(Sinc_vant_sing_ground{1}(:))));
    end
    bestspec = find(error_list==min(error_list));
    Sinc_singfilt{bestspec} = Sinc{bestspec};
    Sinc_vant_singfilt{bestspec} = Sinc_vant{bestspec};
else
    for obj_i=1:length(Sinc_vant_ground)
        error_list = [];
        for spec=1:length(Sinc_vant)
            error_list(spec) = nansum(abs(Sinc_vant{spec}/norm(Sinc_vant{spec}(:))-Sinc_vant_ground{obj_i}/norm(Sinc_vant_ground{obj_i}(:))));
        end
        bestspec = find(error_list==min(error_list));
        Sinc_singfilt{bestspec} = Sinc{bestspec};
        Sinc_vant_singfilt{bestspec} = Sinc_vant{bestspec};
    end
end


%% smooth precon results
smooth_precon = 0; %whether to smooth results across vantages
dw_smooth = msdfov_params.linc_msdfov_smooth*50;
if smooth_precon 
    for spec_c = 1:length(Sinc_precon)
        f = fit((1:length(Sinc_precon{spec_c}))',Sinc_precon{spec_c},'smoothingspline','SmoothingParam',dw_smooth);
        Sinc_precon{spec_c} = f(1:length(Sinc_precon{spec_c}));

        f = fit((1:length(Sinc_noprecon{spec_c}))',Sinc_noprecon{spec_c},'smoothingspline','SmoothingParam',dw_smooth);
        Sinc_noprecon{spec_c} = f(1:length(Sinc_noprecon{spec_c}));
    end
end

%% jade extra
%if filtering the resulting components
if W_filter
    for i=1:numcomp_jade
       wtemp = imgaussfilt(reshape(Wjade(:,i),size(l_field{1})),[5,5]);
       Wjade(:,i) = wtemp(:);
    end
end

%identify and remove clutter based on power
if strcmp(msbss_params.clutter_removal,"pow") 
    for k=1:(msbss_params.num_clutter+addone) %repeat for num_clutter
        [Wjade,Hjade] = remove_strongest_component(Wjade,Hjade,0);
    end
end

%smooth L_incs (if desired)
smooth_jade_linc = 0;
linc_smooth = .1;
if smooth_jade_linc
    for i=1:length(Sinc_jade_best)
        f = fit((1:length(Sinc_jade_best{i}))',Sinc_jade_best{i}','smoothingspline','SmoothingParam',linc_smooth);
        Sinc_jade_best{i} = f(1:length(Sinc_jade_best{i}));
    end
end

%% mscpa extra
%set which spectrums to use
msdfov_specs = 1:num_specs;
if isfield(msdfov_params,'msdfov_specs') && sum(msdfov_params.msdfov_specs) ~= 0
    msdfov_specs = msdfov_params.msdfov_specs;
end
color_dict = {};
if isfield(msdfov_params,'color_dict')
    color_dict = msdfov_params.color_dict;
end

%% get spatial data only measurements 
L_inc_spat = {};
for spec=1:num_specs
    [Fxr,Fyr] = gradient(lfield{spec});
    spot_ang = 1; %for single angle

    %for single angle
    L_inc_spat_dw{spec} = Fyr(1:(num_vants-1),spot_ang)';
    % spatmod = smooth(Fy(1:299,80));
    
    %find spatial
    L_inc_spat{spec} = dw_to_inc2(L_inc_spat_dw{spec},brdf_inc_angs_up,vant_pos,d,y_l);
end

%spat spec
spat_spec = zeros(1,num_specs);
for i=1:num_specs
    spat_spec(i) = norm(L_inc_spat{i});
end
spat_spec = spat_spec/max(spat_spec);

% figure,imagesc(L_inc_spat{1})

%% ls recon stuff
%NOTE: because of some numerical differences, all of these reconstructions
%have to be rescaled to match DFOV better. This is the cause of the
%comp_peri_coeff variable.

%forward mod: numVants x camera_angles x length(inc_angles)
do_ls = 0;
if do_ls   
    run_ls_script
else
    L_inc_ls{1} = nan;
    L_inc_g_ls = nan;
    ls_spec = nan;
    comp_peri_coeff = nan;
end

%smoothing
if do_ls
    for i=1:length(L_inc_ls)
        f = fit((1:length(L_inc_ls{i}))',L_inc_ls{i}','smoothingspline','SmoothingParam',smooth_param);
        L_inc_ls{i} = f(1:length(L_inc_ls{i}))';
    end
end

%% recon metrics for LS stuff
if do_ls
    L_inc_ls_up = {};
    total_spots2 = length(L_inc_ls{1});
    for spec_c = 1:numSpecs
        L_inc_ls_up{spec_c} = interp1(L_inc_ls{spec_c}, linspace(1,total_spots2,num_vants));
    %     L_inc_g_ls_up = interp1(sum_cell({L_inc_g_ls},inten_coeff'.*sum(obj_spec,2)), linspace(1,total_spots2,total_spots));          
    end
else
    for spec_c = 1:numSpecs
        L_inc_ls_up{spec} = nan;
    end
end

%each
%least-squares/computational periscopy metrics
if do_ls
    Sinc_g_down = Sinc_g{ground_i}(round(linspace(1,length(Sinc_g{ground_i}),length(Sinc_ls{1,1}))));
    total_spots2 = length(Sinc_ls{1,1});

    %need to upsample so that the error metrics align better with each
    %other
    Sinc_ls_tot = sum_cell(Sinc_ls);
    Sinc_ls_tot_up = interp1(Sinc_ls_tot, linspace(1,total_spots2,num_vants));
    Sinc_g_ls_tot_up = interp1(sum_cell({Sinc_g_ls},inten_coeff'.*sum(spec_ground,2)), linspace(1,total_spots2,num_vants));
      %original sample
%         ls_best_mse(spec_c) = nansum((Sinc_ls{spec_c}*comp_peri_coeff*coeff_multiple-Sinc_g_ls*comp_peri_coeff).^2)/total_spots2;
%         ls_best_shape_mse(spec_c) = nansum((Sinc_ls{spec_c}/norm(Sinc_ls{spec_c})-Sinc_g_ls/norm(Sinc_g_ls)).^2)/total_spots2;
%         ls_best_pow_dif(spec_c) = power_diff(Sinc_ls{spec_c}*comp_peri_coeff*coeff_multiple,Sinc_g_ls*comp_peri_coeff);

        %up sample
    recmet_dict.each.mse.ls(ground_i) = nansum((Sinc_ls_tot_up*comp_peri_coeff*coeff_multiple-Sinc_g_ls_tot_up*comp_peri_coeff).^2)/num_vants;
    recmet_dict.each.shape.ls(ground_i) = nansum((Sinc_ls_tot_up/norm(Sinc_ls_tot_up)-Sinc_g_ls_tot_up/norm(Sinc_g_ls_tot_up)).^2)/num_vants;
    recmet_dict.each.powdif.ls(ground_i) = power_diff(Sinc_ls_tot_up*comp_peri_coeff*coeff_multiple,Sinc_g_ls_tot_up*comp_peri_coeff);
else
    recmet_dict.each.mse.ls(ground_i) = nan; 
    recmet_dict.each.shape.ls(ground_i) = nan; 
    recmet_dict.each.
    
    powdif.ls(ground_i) = nan; 
end


%% lfield difs
ldif_dict = struct();
ldif_dict.lfield_normmse = lfield_normmse;
ldif_dict.lfield_nmf_normmse = lfield_nmf_normmse;
ldif_dict.lfield_bss_normmse = lfield_bss_normmse;
ldif_dict.lfield_sobi_normmse = lfield_sobi_normmse;
ldif_dict.lfield_jade_normmse = lfield_jade_normmse;
ldif_dict.lfield_pca_normmse = lfield_pca_normmse;
