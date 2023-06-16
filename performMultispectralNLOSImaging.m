%Code to perform the MS-BSS and MS-CPA pipelines in the paper "Isolating
%Signals in Passive Non-Line-of-Sight Imaging using Spectral Content." This
%code sets up the entire workflow and keeps track of experimental
%variables. It is meant to be called upon separately after the variables
%have been set.
%
%Make sure proper 
%--------------------------------------------------------------------------
% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - all_params -- struct that holds more structs, each which contain
% experimental and scene parameters
% -- all_params.scene_params -- parameters concerning the scene
% -- all_params.ground_params -- parameters concerning the ground truth
% -- all_params.brdf_params -- parameters concerning the scattering properties
% (BRDF) of the wall
% -- all_params.recon_params -- parameters concerning reconstructing the
% hidden scene using light field reconstruction
% -- all_params.msbss_params -- parameters concerning the MS-BSS algorithm
% -- all_params.mscpa_params -- parameters concerning the MS-CPA algorithm
% -- all_params.lsrecon_params -- parameters concerning an alternative
% reconstruction method (used for the "optimized preconditioning"
% comparison in the paper)
%- lfield_locations -- cell containing strings which are the locations of
%the light fields
%
%--------------------------------------------------------------------------
% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All outputs are dictionaries which contain the results for each method
% - Sinc_dict -- reconstruction results (in radiance vs incident angle)
% - reconmet_dict -- reconstruction metrics results
% - distance_dict -- reconstruction distance/residuals
% - separated_dict -- contains some useful examples to visualize how each
% BSS method is unmixing the light field
% - predictmet_dict -- metrics which try to predict how hard the unmixing is
% - vars_dict -- important variables


function[recmet_dict,dist_dict,predmet_dict,Sinc_dict,separated_dict,vars_dict] = ...
    performMultispectralNLOSImaging(all_params, lfield_locations)

%% add important paths
% addpath(genpath('C:\Users\hashe037\Desktop\DARPA_REVEAL\multispec_work\auxillary_functions'));
% addpath(genpath('C:\Users\hashe037\Desktop\DARPA_REVEAL\multispec_work\perform_multispec_sigsep'));
% addpath(genpath('C:\Users\hashe037\Desktop\DARPA_REVEAL\multispec_work\scripts'))


%% expand all_params
scene_params = all_params.scene_params;
brdf_params = all_params.brdf_params;
recon_params = all_params.recon_params;
msbss_params = all_params.msbss_params;
mscpa_params = all_params.mscpa_params;
lsrecon_params = all_params.lsrecon_params;
ground_params = all_params.ground_params;
meas_params = struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Setup Parameters and Variables %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load in and preprocess light fields
%most light field and similar parameters are in "meas_params" structure
[scene_params,meas_params] = loadinProcessLfields(ground_params,scene_params,lfield_locations);

fprintf('Light fields and parameters loaded in! \n')

%whether to do single non-uniformly-colored ground
%useful for real-life experiments since this allows a single ground source
%but its spectrum can change vantage-by-vantage (non-uniformly colored)
%repeat many of the previous preprocessing
if ground_params.single_ground
    meas_params.lfield_sing_ground = loadinProcessSingleground(ground_params,meas_params,scene_params);
    fprintf('Using single ground parameters \n')
else %dummy variable
    meas_params.lfield_sing_ground = {zeros(5)};
end

%% find BRDF and corresponding SVD components
%useful for reconstruction and preconditioning (in next section)
[scene_params,brdf_params] = defineBrdfVariables(scene_params,brdf_params);

%% perform preconditioning (Section 4.2 of paper)
% Projecting light fields onto the BRDF singular vectors and performing differential
% Also detailed in "Passive non-line-of-sight imaging using plenoptic
% information" (Optics Express 2020)
[meas_params,recon_params.dalpha_occ,recon_params.diff_occ] ...
    = performPreconditioning(scene_params,brdf_params,meas_params,ground_params,recon_params);


%% calculate predictive metrics to estimate the difficulty of NLOS imaging problem
[predmet_dict] = ...
    findPredictiveMetrics(meas_params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Perform Spectral-Agnostic Reconstructions %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Spectral-agnostic light field reconstructions 
[results_agnostic] = ...
    performAgnosticRecons(scene_params,brdf_params,meas_params,ground_params,recon_params);

%% Best single-spectra measurement
%Naive approach, find the spectral reconstruction that has lowest
%reconstruction residual/distance

[Sinc_vant_singfilt,Sinc_singfilt,bestspec] ...
    = findBestSingleSpectra(results_agnostic.Sinc, results_agnostic.Sinc_vant, results_agnostic.dist_agnostic);
results_agnostic.Sinc_vant_singfilt = Sinc_vant_singfilt;
results_agnostic.Sinc_singfilt = Sinc_singfilt;
results_agnostic.bestspec = bestspec;

fprintf('Single spectra with lowest residual: %i \n',bestspec);


%% perform "optimized preconditiong" method
% perform a modification of the "optimized preconditioning" method as
% described by "Fast Computational Periscopy in Challenging Ambient Light 
% Conditions through Optimized Preconditioning" (ICCP 2021).
% This requires prior knowledge or simulation of clutter scattering 
% structures (lfield_clut). In this work, we found the best
% way to implement it is to perform on differential light field
% measurements (diff(lfield,1))

if recon_params.do_optprecon
    [results_precon] = performOptPrecon(meas_params,...
        scene_params,brdf_params,ground_params,lsrecon_params);
    msbss_params.lfield_precon = results_precon.lfield_precon; %for opt precon preconditioning

    fprintf('Done with Optimized Preconditioning \n')
else
    results_precon = struct();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Perform Multispectral Blind Source Separation (MS-BSS) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% run JADE ICA Algorithm
run_jade = 1;

if run_jade
    [results_jade] = performJadeBss(meas_params,scene_params,brdf_params,recon_params,msbss_params,results_agnostic);
else
    results_jade = struct();
end

fprintf('Done with MS-BSS! \n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Perform Multispectral Closest-Points Algorithm (MS-CPA) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% redo certain parameters for MS-CPA
%MS-CPA can handle more SVD vectors and performs better with them. Also
%have different smoothing of alpha parameters
brdf_params.vector_num = mscpa_params.mscpa_vector_num;
brdf_params.vectors_remove = mscpa_params.mscpa_vectors_remove;
[scene_params,brdf_params] = defineBrdfVariables(scene_params,brdf_params);

[meas_params_mscpa,recon_params.dalpha_occ_mscpa,~] ...
    = performPreconditioning(scene_params,brdf_params,meas_params,ground_params,recon_params);

meas_params.alpha_meas_mscpa = meas_params_mscpa.alpha_meas;
meas_params.dalpha_meas_mscpa = meas_params_mscpa.dalpha_meas;
meas_params.nonsmooth_alpha_meas_mscpa = meas_params_mscpa.nonsmooth_alpha_meas;
meas_params.nonsmooth_dalpha_meas_mscpa = meas_params_mscpa.nonsmooth_dalpha_meas;
meas_params.alpha_ground_mscpa = meas_params_mscpa.alpha_ground;
meas_params.dalpha_ground_mscpa = meas_params_mscpa.dalpha_ground;
meas_params.alpha_clutter_mscpa = meas_params_mscpa.alpha_clutter;
meas_params.dalpha_clutter_mscpa = meas_params_mscpa.dalpha_clutter;

clear meas_params_mscpa

%% run MS-CPA pipeline and algorithm
do_mscpa = 1;
if do_mscpa
    [results_mscpa] = performMscpa(ground_params,scene_params,...
        brdf_params,meas_params,recon_params,mscpa_params,results_agnostic);
    fprintf('Done with MS-CPA! \n')
else
    results_mscpa = struct()
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Post-Proccessing and Finding Reconstruction Metrics  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% post-proccessing
%smooth reconstruct Sincs (mostly for looks but also to handle high
%frequency noise)
Sinc_smooth = recon_params.Sinc_smooth_param;%.1; %smoothing factor
if recon_params.smooth_Sincs
    [meas_params,results_agnostic,results_precon,results_jade,results_mscpa] ...
        = SmoothSincs(Sinc_smooth,meas_params,results_agnostic,results_precon,results_jade,results_mscpa);
end

%take absolute value of Sincs (helps account for flipping)
if recon_params.abs_Sincs
    [meas_params,results_agnostic,results_precon,results_jade,results_mscpa] ...
        = absoluteSincs(meas_params,results_agnostic,results_precon,results_jade,results_mscpa);
end


%% Find reconstruction metrics
[recmet_dict]  = findReconstructionMetrics(ground_params,recon_params,msbss_params,...
    results_agnostic,results_precon,results_jade,results_mscpa);

%% Construct other dictionaries in more friendly-final format

[vars_dict,dist_dict,Sinc_dict,separated_dict] = createOutputDictionaries(ground_params,brdf_params,meas_params,msbss_params,...
    results_agnostic,results_precon,results_jade,results_mscpa);





