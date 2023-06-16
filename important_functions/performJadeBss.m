% COMMENT FURTHER
% 
%This code implements a pipeline that peforms blind source separation (BSS)
%of multispectral light fields using an indepdendent component analysis
%(ICA) method called Joint Approximate Diagonalization of Eigenmatrices 
%(JADE) which estimates independence by the fourth order cumulant. JADE was
%originally theorized in "Blind beamforming for non-Gaussian Signals" (IEE
%1993) by the authors J. Cardoso and A. Souloumiac. All credit goes to them
%for the unmixing code, and the link to the MATLAB code can be found here:
%http://www2.iap.fr/users/cardoso/compsep_classic.html. 
%
%The pipeline we implement goes as follows (explained in Section 4 of paper):
% 1. Create mixture with predefined preconditioner
% 2. Perform unmixing with JADE algorithm
% 3. Find which components have highest reconstruction residual and remove
% them
%
%--------------------------------------------------------------------------
% Outputs (inputted in structure) -----------------------------------------
% Xjade -- observation matrix for input to JADE
% Wjade,Hjade -- BSS components/coefficients for JADE
% l_field_jade -- light fields of components reconstructed with JADE
% L_inc_dw_jade,L_inc_jade -- reconstructed components in same order as Wjade
% components
% L_inc_jade_best,jade_best_i -- Linc and indice that corresponds to the
% best matches to grounds in L_inc_g
%

function[results_jade] = performJadeBss(meas_params,scene_params,brdf_params,recon_params,msbss_params,results_agnostic)

lfield = meas_params.lfield;
alpha_meas = meas_params.alpha_meas;
dalpha_meas = meas_params.dalpha_meas;
nonsmooth_alpha_meas = meas_params.nonsmooth_alpha_meas;
nonsmooth_dalpha_meas = meas_params.nonsmooth_dalpha_meas;
nonsmooth_Sinc_vant = results_agnostic.nonsmooth_Sinc_vant;
Sinc_ground = results_agnostic.Sinc_ground;

%% preconditioning
Xjade = find_bss_mixture(lfield,nonsmooth_alpha_meas,nonsmooth_dalpha_meas,nonsmooth_Sinc_vant, ...
    msbss_params.precon);

%% run actual unmixing process which separates mixture Xjade = Wjade*Hjade
%Wjade is scattering structures, Hjade is spectral coefficients
[Hjade,~] = jade_modified(Xjade,msbss_params.numcomp); %Here Hjade is the unmixing matrix
Wjade = Hjade*Xjade;
Wjade = Wjade';
Hjade = pinv(Hjade)'; %turn into mixing matrix for scaling

%perform reconstruction of hidden scene
[lfield_jade,alpha_jade,dalpha_jade,Sinc_vant_jade,Sinc_jade,dist_jade] =...
    perform_conversion_reconstruction_W(Wjade,lfield,alpha_meas,dalpha_meas,msbss_params,scene_params,recon_params,brdf_params);

%% remove clutter from color
%experimental code to adjust color based on the measured clutter's color
if isfield(msbss_params,'color_adjust_fact')

    %construct normed Hjade (each spectral measurement is normed)
    Hjade_normed = Hjade;
    Hpower = [];
    Hsign = [];
    for i=1:size(Hjade,1)
        Hpower(i) = norm(Hjade(i,:));
        Hsign(i) = sign(Hjade(i,3));
        Hjade_normed(i,:) = (abs(Hjade(i,:)))/norm(Hjade(i,:));
    end
    
    %find expected clutter spectrum by assuming the clutter is overwhelming
    %and therefore the main color in the light fields
    clut_spec_est = [];
    for spec=1:length(l_field)
        clut_spec_est(spec) = norm(l_field{spec},'fro');
    end
    clut_spec_est = clut_spec_est/norm(clut_spec_est(:));

    %using the color adjust factor, subtract the color from Hjade_normed
    %and also try to add the total power back in
    mult_fact = msbss_params.color_adjust_fact;
    clut_spec_est = mult_fact*clut_spec_est/norm(clut_spec_est);
    for i=1:size(Hjade,1)
        Hjade_normed(i,:) = Hjade_normed(i,:)-clut_spec_est;
    %     Hjade_normed(i,:) = Hsign(i)*Hjade_normed(i,:)*Hpower(i);
        Hjade_normed(i,:) = Hsign(i)*abs(Hjade_normed(i,:))*Hpower(i);
    end
    Hjade = Hjade_normed;
end

%% remove clutter elements by analyzing residuals (dist)

%for debugging
dist_sum = [];
for i=1:length(dist_jade)
    dist_sum(i) = sum(dist_jade{i});
end

%remove clutter from distance metric and then rerun reconstruction
for k=1:(msbss_params.num_clutter) %repeat for num_clutter
    [Wjade,Hjade] = remove_worseresidual_component(Wjade,Hjade,dist_jade);
    [lfield_jade,~,dalpha_jade,Sinc_vant_jade,Sinc_jade,dist_jade] = perform_conversion_reconstruction_W(Wjade,lfield,alpha_meas,dalpha_meas,msbss_params,scene_params,recon_params,brdf_params);
end

%find best fitting for each ground truth to better compare
[Wjade,Hjade,lfield_jade,Sinc_jade_best,jade_best_i] = find_flip_best_fit(Wjade,Hjade,lfield_jade,Sinc_jade,Sinc_ground,1);

% Hjade = -1*Hjade;

%% create output structure
results_jade = struct();

results_jade.Sinc_jade_best = Sinc_jade_best;
results_jade.jade_best_i = jade_best_i;
results_jade.Wjade = Wjade;
results_jade.Hjade = Hjade;
results_jade.lfield_jade = lfield_jade;
results_jade.Sinc_vant_jade = Sinc_vant_jade;
results_jade.Sinc_jade = Sinc_jade;
results_jade.dalpha_jade = dalpha_jade;
results_jade.alpha_jade = alpha_jade;
results_jade.dist_jade = dist_jade;
