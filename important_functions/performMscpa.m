% MAKE COMMENTS
%
%run MS-DFOV in script form to put variables to workspace

function[results_mscpa] = performMscpa(ground_params,scene_params,brdf_params,meas_params,recon_params,mscpa_params,results_agnostic)

%% set some MS-DFOV parameters
%number of spectrums to use (more is better but slower)

dalpha_ground = meas_params.dalpha_ground_mscpa;
spec_ground = meas_params.spec_ground;

%% find expected object spectra for each change in vantage
%this is gamma_* in the paper (CFOV spectra)

%which ground object is "brightest" at each vantage
[vant_objs,vant_strengths] = find_vantobjs(dalpha_ground,spec_ground);

%given the brightest objects at each vantage, find the "gamma_cfov", which
%is the gamma_* in the paper
gamma_cfov = zeros(size(dalpha_ground{1},1),size(spec_ground,2));
for i=1:size(gamma_cfov,1)
    gamma_cfov(i,:) = spec_ground(vant_objs(i),:);
end

%if we have a non-uniform colored ground source, just take the spectrum of
%the ground 
if ground_params.single_ground
    gamma_cfov = results_agnostic.spec_sing_ground;
end

%% find Sinc using the CPA algorithm
spec_choice = mscpa_params.spec_choice; %what spec to choose as base (if method==0)
spec_choice_method = mscpa_params.spec_choice_method; %how to choose base spec for algorithm
msdfov_do_diff = mscpa_params.do_diff; %whether to perform with differential or dalpha

if spec_choice_method == 1 %find spec that minimizes MS-CPA distance (best fits assumptions)
    temp_dis_method = 'residnorm';
    for spec=1:length(meas_params.lfield) %each spectra
        [Sinc_mscpa_temp(spec,:),Sinc_vant_mscpa_temp(spec,:),dist_mscpa_temp(spec,:),~] = ...
            perform_cpa_reconstruction(gamma_cfov,scene_params,brdf_params,meas_params,...
            recon_params,mscpa_params,spec,temp_dis_method,mscpa_params.do_diff);
        dist_mscpa_summed(spec) = nansum(dist_mscpa_temp(spec,:));
    end
    [~,spec_choice] = min(dist_mscpa_summed);
    % spec_choice = spec_choice(1); %make sure no duplicates

    %run with ideal spec
    [Sinc_mscpa,Sinc_vant_mscpa,dist_mscpa,~] = ...
            perform_cpa_reconstruction(gamma_cfov,scene_params,brdf_params,meas_params,...
            recon_params,mscpa_params,spec_choice,mscpa_params.dis_method,mscpa_params.do_diff);
    
elseif spec_choice_method == 2 %weighted measurements according to distance
    error('Need to still do \n')
else %prechosen spec
    [Sinc_mscpa,Sinc_vant_mscpa,dist_mscpa,~] = ...
            perform_cpa_reconstruction(gamma_cfov,scene_params,brdf_params,meas_params,...
            recon_params,mscpa_params,spec_choice,mscpa_params.dis_method,mscpa_params.do_diff);
end
fprintf('MS-CPA spec for offset: %i \n',spec_choice)

%invert dueto nature of gradient
Sinc_mscpa = -Sinc_mscpa;
Sinc_vant_mscpa = -Sinc_vant_mscpa;


%% attenuate Sinc based on the distance parameter
if isfield(mscpa_params,'mscpa_dist_min')
    mscpa_dist_min = mscpa_params.mscpa_dist_min; %.2 for monitor experiments
else
    mscpa_dist_min = .2;
end
exp_fact = 5; %5 for monitor experiments
if mscpa_params.attenuate_ondist
   msdfov_dist2 = smooth(dist_mscpa)';
   atten_coeff = (min(1-mscpa_dist_min,1-msdfov_dist2)/(1-mscpa_dist_min)).^exp_fact; %assuming max distance is 1
   Sinc_vant_mscpa_noatten = Sinc_vant_mscpa; %save in case of viewing
   Sinc_vant_mscpa = Sinc_vant_mscpa_noatten.*atten_coeff; %attenuate
   Sinc_mscpa = incvant_to_incang(Sinc_vant_mscpa,brdf_params.brdf_inc_angs_up,...
       scene_params.vant_pos,scene_params.occ_x,scene_params.occ_y);
end

%% smooth and construct multispectral version of msdfov
smooth_mscpa = 1; %whether to smooth results across vantages
vant_smooth = mscpa_params.vant_smooth;
Sinc_vant_mscpa_multispec = repmat(Sinc_vant_mscpa,size(spec_ground,2),1).*gamma_cfov(1:length(Sinc_vant_mscpa),:)';    

%smooth L_inc_msdfov data if needed
if smooth_mscpa %only do it for distance relaxation
    for spec_c = 1:size(Sinc_vant_mscpa_multispec,1)
        f = fit((1:size(Sinc_vant_mscpa_multispec,2))',Sinc_vant_mscpa_multispec(spec_c,:)','smoothingspline','SmoothingParam',vant_smooth);
        Sinc_vant_mscpa_multispec(spec_c,:) = f(1:size(Sinc_vant_mscpa_multispec,2));
    end
end

%convert to cells like the rest
Sinc_mscpa_multispec = {};
for spec_c=1:size(Sinc_vant_mscpa_multispec,1) %each spectra
    Sinc_mscpa_multispec{spec_c} = incvant_to_incang(Sinc_vant_mscpa_multispec(spec_c,:),brdf_params.brdf_inc_angs_up,...
       scene_params.vant_pos,scene_params.occ_x,scene_params.occ_y);
end
Sinc_vant_mscpa = sum(Sinc_vant_mscpa_multispec,1); %sum over spectra
Sinc_mscpa = incvant_to_incang(Sinc_vant_mscpa,brdf_params.brdf_inc_angs_up,...
       scene_params.vant_pos,scene_params.occ_x,scene_params.occ_y);

%measure spectrum
spec_mscpa = zeros(1,size(spec_ground,2));
for vant = 1:size(Sinc_vant_mscpa_multispec,2)
    spec_mscpa = spec_mscpa + abs(Sinc_vant_mscpa_multispec(:,vant)'); %??? is this right
end
spec_mscpa = abs(spec_mscpa)/max(abs(spec_mscpa));


%% create results dictionary
results_mscpa = struct();

results_mscpa.Sinc_mscpa = Sinc_mscpa;
results_mscpa.Sinc_vant_mscpa = Sinc_vant_mscpa;
results_mscpa.spec_mscpa = spec_mscpa;
results_mscpa.dist_mscpa = dist_mscpa;
results_mscpa.Sinc_mscpa_multispec = Sinc_mscpa_multispec;
% results_mscpa.Sinc_mscpa = Sinc_mscpa;
% results_mscpa.Sinc_mscpa = Sinc_mscpa;
% results_mscpa.Sinc_mscpa = Sinc_mscpa;
% results_mscpa.Sinc_mscpa = Sinc_mscpa;






