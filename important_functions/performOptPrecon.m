%Run a modified version of "optimized preconditioning" motivated by the
%paper "Fast Computational Periscopy in Challenging Ambient Light 
% Conditions through Optimized Preconditioning" (ICCP 2021).
%
% While this method does not use TV regularization, it does implement the
% optimized preconditioner with knowledge of the light field of the clutter
% a priori. It also has been adapted to light fields rather than a single
% photo, which is not too difficult. This does solve the for the hidden
% scene slightly different than the other methods (previous methods analyze
% vantage-by-vantage while this analyzes all vantages at once), so the
% ground truth has to be reconstructed to compare accurately.
%
%--------------------------------------------------------------------------
% Output (in structure) ---------------------------------------------------
% - Sinc_precon -- reconstructed incident light after preconditioning
% - Sinc_precon_clutter -- reconstructed clutter sources (should be close
% to null if optimized preconditioner is working correctly).
% - S_inc_precon_ground -- reconstructed ground truth
% - Sinc_precon_single_ground -- reconstructed single ground (if using)
% - lfield_precon -- reconstructed light fields after preconditioning
% - lfield_precon_clutter -- reconstructed clutter light fields (should be
% close to null)
%
function[results_precon] = performOptPrecon(meas_params,scene_params,brdf_params,ground_params,lsrecon_params)

lfield = meas_params.lfield;
lfield_ground = meas_params.lfield_ground;
lfield_clutter = meas_params.lfield_clutter;
lfield_sing_ground = meas_params.lfield_sing_ground;

%% beginning params
do_diff = 1; %take differential when performing opt precon (recommended)

%define weights for clutter matrix
w = lsrecon_params.w;

%find background matrix
B_mat = construct_background(lfield_clutter,w,do_diff);
% B_mat(B_mat<1e-4)=0; %set values to 0 to prevent errors

%% find forward model and optimized preconditioner
[~,A_2d] = find_precon_forward_mod(do_diff,scene_params,brdf_params);
A_precon = horzcat(A_2d,B_mat);

%perform inversion of forward model
Ainv_precon = pinv(A_precon);
Ainv_precon(:,1:55) = 0; %prevent numerical errors
Ainv_precon = Ainv_precon(1:end-length(lfield_clutter),:); %hits is the optimized preconditioner!
%% find ground truth to compare to
%Since we are solving for the hidden scene in a way slightly different than
%performLfieldRecon.m, we want to compare to the true ground truth for this
%method to get the best comparisons
Sinc_precon_ground = {};
lfield_precon_ground = {};
for obj_ind=1:length(lfield_ground)
    lfield_mod = lfield_ground{obj_ind};
    if do_diff
        lfield_mod = smoothdata(diff(lfield_mod,1),1);
    end
    lfield_mod = reshape(lfield_mod,numel(lfield_mod),1);
    Sinc_precon_ground{obj_ind} = -1*pinv(A_2d)*lfield_mod;
    lfield_precon_ground{obj_ind} = -1 * reshape((A_2d * Sinc_precon_ground{obj_ind}),size(lfield{1},1)-do_diff,size(lfield{1},2)); 
end

%if just doing single (non-uniformly-colored) ground
if ground_params.single_ground
    Sinc_precon_single_ground = {};
    for spec=1:length(lfield)
        lfield_mod = lfield_sing_ground{spec};
        if do_diff
            lfield_mod = smoothdata(diff(lfield_mod,1),1);
        end
        lfield_mod = reshape(lfield_mod,numel(lfield_mod),1);
        Sinc_precon_single_ground{spec} = -1*pinv(A_2d)*lfield_mod;
    end
else
    Sinc_precon_single_ground = {zeros(5)};
end

%% perform for clutter measurements to make sure it is preconditioning well
%response should be close to null
Sinc_precon_clutter = {};
lfield_precon_clutter = {};
for clut_ind=1:length(lfield_clutter)
    lfield_mod = lfield_clutter{clut_ind};

    %60db means noise should be about 1000th of the signal
    lfield_mod = lfield_mod+3*randn(size(lfield_mod));
    if do_diff
        lfield_mod = smoothdata(diff(lfield_mod,1),1);
    end
    lfield_mod = reshape(lfield_mod,numel(lfield_mod),1);
    Sinc_precon_clutter{clut_ind} = -1*Ainv_precon*lfield_mod;
    lfield_precon_clutter{clut_ind} = -1 * reshape((A_2d * Sinc_precon_clutter{clut_ind}),size(lfield{1},1)-do_diff,size(lfield{1},2)); 
    lfield_precon_clutter{clut_ind} = replace_zeros_fromabove(lfield_precon_clutter{clut_ind});
end

%% perform for each spectral measurement
Sinc_precon = {};
lfield_precon = {};
for spec=1:length(lfield)
    lfield_mod = lfield{spec};
    if do_diff
        lfield_mod = diff(lfield_mod,1);
        % lfield_mod = smoothdata(diff(lfield_mod,1),1);
    end
    lfield_mod = reshape(lfield_mod,numel(lfield_mod),1);
    Sinc_noprecon{spec} = -1*pinv(A_2d)*lfield_mod;
    Sinc_precon{spec} = -1*Ainv_precon*lfield_mod;
    lfield_precon{spec} = reshape((A_2d * Sinc_precon{spec}),size(lfield{1},1)-do_diff,size(lfield{1},2)); 
    lfield_precon{spec} = replace_zeros_fromabove(lfield_precon{spec});
    lfield_noprecon{spec} = reshape((A_2d * Sinc_noprecon{spec}),size(lfield{1},1)-do_diff,size(lfield{1},2)); 
end

% spec_sample = 1;
% figure,plot(sum_cell(Sinc_precon))
% figure,plot(Sinc_precon_g_all),xlim([0,200]),title('Ground')
% figure,plot(Sinc_noprecon{spec_sample}),xlim([0,200]),title('Normal')
% figure,plot(Sinc_precon{spec_sample}),xlim([0,200]),title('Precon')


%% smooth lfields
smooth_lfield_precon = 1; %whether to smooth results across vantages
dw_smooth = .8; %.1
ang_smooth = .5; %.1
if smooth_lfield_precon 
    for spec_c = 1:length(lfield_precon)
        for i=1:size(lfield_precon{spec_c},2)
            f = fit((1:size(lfield_precon{spec_c},1))', lfield_precon{spec_c}(:,i),'smoothingspline','SmoothingParam',ang_smooth);
            lfield_precon{spec_c}(:,i) = f(1:size(lfield_precon{spec_c},1));
        end
        for i=1:size(lfield_precon{spec_c},1)
            f = fit((1:size(lfield_precon{spec_c},2))', lfield_precon{spec_c}(i,:)','smoothingspline','SmoothingParam',dw_smooth);
            lfield_precon{spec_c}(i,:) = f(1:size(lfield_precon{spec_c},2));
        end
    end
end

%% upsample reconstructions to match others
Sinc_precon_up = {};
Sinc_noprecon_up = {};
num_vants_precon = length(Sinc_precon{1});
num_vants = size(brdf_params.brdf_inc_angs_up,2);
for spec_c = 1:length(Sinc_precon)
    Sinc_precon_up{spec_c} = interp1(Sinc_precon{spec_c}, linspace(1,num_vants_precon,num_vants));
    Sinc_noprecon_up{spec_c} = interp1(Sinc_noprecon{spec_c}, linspace(1,num_vants_precon,num_vants));
end

for ground_i = 1:length(Sinc_precon_ground)
    Sinc_precon_ground_up{ground_i} = interp1(Sinc_precon_ground{ground_i}.*sum(meas_params.spec_ground(ground_i,:)), linspace(1,num_vants_precon,num_vants));
end

%% create precon results structure
results_precon = struct();
results_precon.Sinc_precon = Sinc_precon;
results_precon.Sinc_noprecon = Sinc_noprecon;
results_precon.Sinc_precon_clutter = Sinc_precon_clutter;
results_precon.Sinc_precon_ground = Sinc_precon_ground;
results_precon.Sinc_precon_single_ground = Sinc_precon_single_ground;
results_precon.lfield_precon = lfield_precon;
results_precon.lfield_precon_clutter = lfield_precon_clutter;
results_precon.Sinc_precon_up = Sinc_precon_up;
results_precon.Sinc_noprecon_up = Sinc_noprecon_up;
results_precon.Sinc_precon_ground_up = Sinc_precon_ground_up;



