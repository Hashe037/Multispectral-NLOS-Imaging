%MAKE COMMENT
%
%find errors for each ground

% function[metric_each_dict,metric_all]  = find_reconstruction_metrics
function[recmet_dict] ...
    = findReconstructionMetrics(ground_params,recon_params,msbss_params,results_agnostic,results_precon,results_jade,results_mscpa)

power_diff = @(x,y) abs(max(abs(x(:)))-max(abs(y(:))))/max(max(abs(x(:))),max(abs(y(:))));
num_vants = length(results_agnostic.Sinc_ground{1}); %total incident divisions

do_precon = isfield(results_precon,'Sinc_precon');
do_jade = isfield(results_jade,'Sinc_jade');
do_mscpa = isfield(results_mscpa,'Sinc_mscpa');

spec_ground = ground_params.spec_ground;
num_ground = size(spec_ground,1);

%% dictionaries to return
recmet_dict = struct(); %reconstruction metrics
recmet_dict.each = struct(); %errors between each object and its best reconstruction
recmet_dict.each.mse = struct(); %non-normalized MSE
recmet_dict.each.shape = struct(); %normalized MSE
recmet_dict.each.powdif = struct(); %difference in "power" where power is defined by us

recmet_dict.total = struct(); %errors between all objects and all reconstructions
recmet_dict.total.mse = struct(); %non-normalized MSE
recmet_dict.total.shape = struct(); %normalized MSE
recmet_dict.total.powdif = struct(); %difference in "power" where power is defined by us

%% upsample Sincls (later is redundant)
if do_precon
    num_vants_precon = length(results_precon.Sinc_precon{1}); %total incident divisions

    Sinc_precon_up = results_precon.Sinc_precon_up;
    Sinc_noprecon_up = results_precon.Sinc_noprecon_up;
    Sinc_precon_ground_up = results_precon.Sinc_precon_ground_up;
end
% if do_precon
%     Sinc_precon_up = {};
%     Sinc_noprecon_up = {};
%     num_vants_precon = length(results_precon.Sinc_precon{1});
%     for spec_c = 1:numSpecs
%         Sinc_precon_up{spec_c} = interp1(results_precon.Sinc_precon{spec_c}, linspace(1,num_vants_precon,num_vants));
%         Sinc_noprecon_up{spec_c} = interp1(results_precon.Sinc_noprecon{spec_c}, linspace(1,num_vants_precon,num_vants));
%     end
% end

%% find metrics between EACH ground and best corresponding reconstruction
%Summed over all spectra

for ground_i=1:num_ground
    Sinc_ground_i = results_agnostic.Sinc_ground{ground_i}.*sum(spec_ground(ground_i,:));
    coeff_multiple = 1;

    %agnostic (sum over all spectra)
    Sinc_tot = sum_cell(results_agnostic.Sinc);
    recmet_dict.each.mse.agnostic(ground_i) = nansum((-Sinc_tot*coeff_multiple+Sinc_ground_i).^2)/num_vants;
    recmet_dict.each.shape.agnostic(ground_i) = nansum((-Sinc_tot/norm(Sinc_tot)+Sinc_ground_i/norm(Sinc_ground_i)).^2)/num_vants;
    recmet_dict.each.powdif.agnostic(ground_i) = power_diff(Sinc_tot*coeff_multiple,Sinc_ground_i);

    %single spectra filtered
    Sinc_singfilt_tot = sum_cell(results_agnostic.Sinc_singfilt);
    recmet_dict.each.mse.singfilt(ground_i) = nansum((-Sinc_singfilt_tot*sum(spec_ground(ground_i,:))*coeff_multiple+Sinc_ground_i).^2)/num_vants;
    recmet_dict.each.shape.singfilt(ground_i) = nansum((-Sinc_singfilt_tot/norm(Sinc_singfilt_tot)+Sinc_ground_i/norm(Sinc_ground_i)).^2)/num_vants;
    recmet_dict.each.powdif.singfilt(ground_i) = power_diff(Sinc_singfilt_tot*sum(spec_ground(ground_i,:))*coeff_multiple,Sinc_ground_i);

    %Jade BSS metrics
    if do_jade
        Hjade = results_jade.Hjade;
        jade_best_i = results_jade.jade_best_i;
        Sinc_jade_best = results_jade.Sinc_jade_best;
        [recmet_dict.each.shape.jade(ground_i),recmet_dict.each.mse.jade(ground_i),recmet_dict.each.powdif.jade(ground_i)] ...
                = calculate_metrics_withbest_single(Sinc_jade_best{ground_i},jade_best_i(ground_i),Hjade,coeff_multiple,Sinc_ground_i,power_diff);
    else
        recmet_dict.each.mse.jade(ground_i) = nan;
        recmet_dict.each.shape.jade(ground_i) = nan;
        recmet_dict.each.powdif.jade(ground_i) = nan;
    end

    %MS-CPA metrics
    %NOTE: MS-CPA does not go object-by-object so this metric is somewhat
    %useless
    if do_mscpa
        Sinc_mscpa = results_mscpa.Sinc_mscpa;
        recmet_dict.each.mse.msdfov(ground_i) = nansum((-Sinc_mscpa*coeff_multiple+Sinc_ground_i).^2)/num_vants;
        recmet_dict.each.shape.msdfov(ground_i) = nansum((-Sinc_mscpa/norm(Sinc_mscpa)+Sinc_ground_i/norm(Sinc_ground_i)).^2)/num_vants;
        recmet_dict.each.powdif.msdfov(ground_i) = power_diff(Sinc_mscpa*coeff_multiple,Sinc_ground_i);
    else
        recmet_dict.each.mse.msdfov(ground_i) = nan;
        recmet_dict.each.shape.msdfov(ground_i) = nan;
        recmet_dict.each.powdif.msdfov(ground_i) = nan;
    end


    %precon metrics
    %remember, solving this in slightly different way so we need to
    %reconstruct it differently as well
    if do_precon        
        recmet_dict.each.mse.precon(ground_i) = nansum((sum_cell(Sinc_precon_up)*coeff_multiple-Sinc_precon_ground_up).^2)/num_vants;
        recmet_dict.each.shape.precon(ground_i) = nansum((sum_cell(Sinc_precon_up)/norm(sum_cell(Sinc_precon_up))-Sinc_precon_ground_up/norm(Sinc_precon_ground_up)).^2)/num_vants;
        recmet_dict.each.powdif.precon(ground_i) = power_diff(sum_cell(Sinc_precon_up)*coeff_multiple,Sinc_precon_ground_up);

        recmet_dict.each.mse.noprecon(ground_i) = nansum((sum_cell(Sinc_noprecon_up)*coeff_multiple-Sinc_precon_ground_up).^2)/num_vants;
        recmet_dict.each.shape.noprecon(ground_i) = nansum((sum_cell(Sinc_noprecon_up)/norm(sum_cell(Sinc_noprecon_up))-Sinc_precon_ground_up/norm(Sinc_precon_ground_up)).^2)/num_vants;
        recmet_dict.each.powdif.noprecon(ground_i) = power_diff(sum_cell(Sinc_noprecon_up)*coeff_multiple,Sinc_precon_ground_up);

    else
        recmet_dict.each.mse.precon(ground_i) = nan; 
        recmet_dict.each.shape.precon(ground_i) = nan; 
        recmet_dict.each.powdif.precon(ground_i) = nan; 
        recmet_dict.each.mse.noprecon(ground_i) = nan; 
        recmet_dict.each.shape.noprecon(ground_i) = nan; 
        recmet_dict.each.powdif.noprecon(ground_i) = nan; 
    end
end


%% Compute MSE between versions and TOTAL ground (all grounds summed together)
%this is better idea of entire CFOV scene
%Summing over all specs

%if we want to combine all non-clutter BSS or not
bss_ver = 1;
if isfield(msbss_params,'bss_ver')
    bss_ver = msbss_params.bss_ver;
end

%choose which ground to use
if ground_params.single_ground
    Sinc_ground_tot = sum_cell(results_agnostic.Sinc_sing_ground);
else
    Sinc_ground_tot = sum_cell(results_agnostic.Sinc_ground,sum(spec_ground,2));
end

%agnostic
Sinc_tot = sum_cell(results_agnostic.Sinc);
recmet_dict.total.mse.agnostic = nansum((-Sinc_tot+Sinc_ground_tot).^2)/num_vants;
recmet_dict.total.shape.agnostic = nansum((-Sinc_tot/norm(Sinc_tot)+Sinc_ground_tot/norm(Sinc_ground_tot)).^2)/num_vants;
recmet_dict.total.powdif.agnostic = power_diff(Sinc_tot,Sinc_ground_tot);

%JADE
if do_jade
    Hjade = results_jade.Hjade;
    jade_best_i = results_jade.jade_best_i;
    if bss_ver==2 %sum over all components not rejected
        Sinc_jade = results_jade.Sinc_jade;
        if recon_params.smooth_Sincs
            for i=1:length(Sinc_jade) %smooth Sinc_jade since we normally do not
                f = fit((1:length(Sinc_jade{i}))',Sinc_jade{i}','smoothingspline','SmoothingParam',recon_params.Sinc_smooth_param);
                Sinc_jade{i} = f(1:length(Sinc_jade{i}))';
            end
        end
        [recmet_dict.total.shape.jade,recmet_dict.total.mse.jade,recmet_dict.total.powdif.jade,Sinc_jade_total] ...
                = calculate_metrics_total_allcomp(Sinc_jade,length(Sinc_jade),1,Hjade,1,Sinc_ground_tot,power_diff,2);
    else %use components that most closely match
        Sinc_jade_best = results_jade.Sinc_jade_best;
        [recmet_dict.total.shape.jade,recmet_dict.total.mse.jade,recmet_dict.total.powdif.jade,Sinc_jade_total] ...
                = calculate_metrics_total(Sinc_jade_best,jade_best_i,1,Hjade,1,Sinc_ground_tot,power_diff,2);
    end
else
    recmet_dict.total.mse.jade = nan;
    recmet_dict.total.shape.jade = nan;
    recmet_dict.total.powdif.jade = nan;   
end

%MS-CPA
if do_mscpa
    Sinc_mscpa_multispec=results_mscpa.Sinc_mscpa_multispec;
    [recmet_dict.total.shape.mscpa,recmet_dict.total.mse.mscpa,recmet_dict.total.powdif.mscpa,Sinc_mscpa_total] ...
        = calculate_msdfov_metrics_total(Sinc_mscpa_multispec,1,1,Sinc_ground_tot,power_diff,2);
else
    recmet_dict.total.mse.msdfov = nan;
    recmet_dict.total.shape.msdfov = nan;
    recmet_dict.total.powdif.msdfov = nan;
end

if do_precon
    if ground_params.single_ground
        Sinc_precon_single_ground = results_precon.Sinc_precon_single_ground;
        Sinc_precon_ground_all_up = interp1(sum_cell(Sinc_precon_single_ground), linspace(1,num_vants_precon,num_vants));
        Sinc_precon_ground_all_up(isnan(Sinc_precon_ground_all_up)) = 0;
    else
        Sinc_precon_ground = results_precon.Sinc_precon_ground;
        Sinc_precon_ground_all_up = interp1(sum_cell(Sinc_precon_ground,sum(spec_ground,2)), linspace(1,num_vants_precon,num_vants));
        Sinc_precon_ground_all_up(isnan(Sinc_precon_ground_all_up)) = 0;
    end
    Sinc_precon_up_summed = sum_cell(Sinc_precon_up);
    Sinc_precon_up_summed(isnan(Sinc_precon_up_summed)) = 0;
    Sinc_noprecon_up_summed = sum_cell(Sinc_noprecon_up);
    Sinc_noprecon_up_summed(isnan(Sinc_noprecon_up_summed)) = 0;

    recmet_dict.total.mse.precon = nansum((Sinc_precon_up_summed-Sinc_precon_ground_all_up).^2)/num_vants;
    recmet_dict.total.shape.precon = nansum((Sinc_precon_up_summed/norm(Sinc_precon_up_summed)-Sinc_precon_ground_all_up/norm(Sinc_precon_ground_all_up)).^2)/num_vants;
    recmet_dict.total.powdif.precon = power_diff(Sinc_precon_up_summed,Sinc_precon_ground_all_up);
    recmet_dict.total.mse.noprecon = nansum((Sinc_noprecon_up_summed-Sinc_precon_ground_all_up).^2)/num_vants;
    recmet_dict.total.shape.noprecon = nansum((Sinc_noprecon_up_summed/norm(Sinc_noprecon_up_summed)-Sinc_precon_ground_all_up/norm(Sinc_precon_ground_all_up)).^2)/num_vants;
    recmet_dict.total.powdif.noprecon = power_diff(Sinc_noprecon_up_summed,Sinc_precon_ground_all_up);
else
    recmet_dict.total.mse.precon = nan;
    recmet_dict.total.shape.precon = nan; 
    recmet_dict.total.powdif.precon = nan; 
    recmet_dict.total.mse.noprecon = nan;
    recmet_dict.total.shape.noprecon = nan; 
    recmet_dict.total.powdif.noprecon = nan; 
end

%special error for single spectra filtered
act_specs = [];
Sinc_singfilt = results_agnostic.Sinc_singfilt;
% for i=1:length(Sinc_vant_singfilt)
%     if abs(Sinc_vant_singfilt{i}) ~= 0
%         act_specs(i) = 1;
%     else
%         act_specs(i) = 0;
%     end
% end

recmet_dict.total.mse.singfilt = nansum((-sum_cell(Sinc_singfilt)+Sinc_ground_tot).^2)/num_vants;
recmet_dict.total.shape.singfilt = nansum((-sum_cell(Sinc_singfilt)/norm(sum_cell(Sinc_singfilt))+Sinc_ground_tot/norm(Sinc_ground_tot)).^2)/num_vants;
recmet_dict.total.powdif.singfilt = power_diff(sum_cell(Sinc_singfilt),Sinc_ground_tot);









