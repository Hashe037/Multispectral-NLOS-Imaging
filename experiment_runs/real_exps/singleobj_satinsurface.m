%Code for running experiment with single real-life experiment 

%% add paths
addpath('..\..\');
addpath(genpath('..\..\utils'));
addpath(genpath('..\..\important_functions'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File/run information
basefolder = dir('..\..\').folder;
expfolder = strcat(basefolder,"\calibration_and_data\real_singleobj_satinsurface\"); 

obj_color = "blue"; %can be blue, red, or green
ground_params = struct(); %holds information about ground truth

if strcmp(obj_color,"blue")
    expname = "blue_obj_run\";
    ground_params.ground_locations = strcat(expfolder,"blue_obj\light_fields\lfield_3.csv");
elseif strcmp(obj_color,"red")
    expname = "red_obj_run\";
    ground_params.ground_locations = strcat(expfolder,"red_obj\light_fields\lfield_3.csv");
else
    expname = "green_obj_run\";
    ground_params.ground_locations = strcat(expfolder,"green_obj\light_fields\lfield_3.csv");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ground information
load(strcat(expfolder,expname,"obj_spec.mat"))
ground_params.spec_ground = obj_spec;

ground_params.single_ground = 0;
ground_params.single_ground_location = "";

%% clutter info
ground_params.clutter_locations = {};
ground_params.clutter_locations{1} = strcat(expfolder,"hidden_side_clutter_1\light_fields\lfield_2.csv");
ground_params.clutter_locations{2} = strcat(expfolder,"hidden_side_clutter_2\light_fields\lfield_2.csv");
ground_params.clutter_locations{3} = strcat(expfolder,"observer_side_clutter_1\light_fields\lfield_2.csv");
load(strcat(expfolder,expname,"clut_spec.mat"))
ground_params.spec_clut = clut_spec;

%% Experiment Parameters
[scene_params, brdf_params, recon_params, msbss_params, mscpa_params, lsrecon_params] ...
    = get_params_real_singleobj_satin(basefolder);

%change parameters slightly based on which object run
if strcmp(obj_color,"blue")
    mscpa_params.mscpa_vector_num = 3; %number of SVD vectors to use
    mscpa_params.mscpa_vectors_remove = [];
    msbss_params.numcomp = 4; %number of components for BSS
    msbss_params.num_clutter = 3; %how many clutter objects t
elseif strcmp(obj_color,"red")
    brdf_params.vector_num = 5; %number of vectors to use for SVD
    brdf_params.vectors_remove = [4]; %vectors to remove because troublesome
    mscpa_params.mscpa_vector_num = 6; %number of SVD vectors to use
    mscpa_params.mscpa_vectors_remove = [4];
    msbss_params.numcomp = 5; %number of components for BSS
    msbss_params.num_clutter = 4; %how many clutter objects t
    recon_params.dis_method = "residnorm_reconweighted"; %metric to use for DFOV
else
    mscpa_params.nullsize_method = 0; %how to determine the nullsize of the CP algorithm
    mscpa_params.mscpa_vector_num = 6; %6 %number of SVD vectors to use
    mscpa_params.mscpa_vectors_remove = [4];
    msbss_params.numcomp = 3; %number of components for BSS
    msbss_params.num_clutter = 2; %how many clutter objects to remove (loops this many times)
end


%% define the different folder locations for each run
%find filenames
foldername = strcat(expfolder,expname,'light_fields\');
filenames = ["lfield_1.csv","lfield_2.csv","lfield_3.csv","lfield_4.csv","lfield_5.csv"];

%load in each object
for file_ind = 1:length(filenames)
    lfield_locations_list(file_ind) = strcat(foldername,filenames(file_ind));
end

%% run algorithms
all_params = struct();
all_params.scene_params = scene_params;
all_params.brdf_params = brdf_params;
all_params.recon_params = recon_params;
all_params.mscpa_params = mscpa_params;
all_params.msbss_params = msbss_params;
all_params.lsrecon_params = lsrecon_params;
all_params.ground_params = ground_params;

%run
[recmet_results,dist_results,predmet_results, ...
    Sinc_results,separated_results,vars_dict] ...
    = performMultispectralNLOSImaging(all_params, lfield_locations_list);
fprintf('Done \n\n')

%% show results

fprintf('Shape Error: LS: %e, Filt: %e, Precon: %e, JADE: %e , MSCPA: %e \n',...
    recmet_results.total.shape.agnostic,recmet_results.total.shape.singfilt,recmet_results.total.shape.precon, ...
    recmet_results.total.shape.jade,recmet_results.total.shape.mscpa)

fprintf('Shape Log Error: LS: %.2f, Filt: %.2f, Precon: %.2f, JADE: %.2f , MSCPA: %.2f \n',...
    log10(recmet_results.total.shape.agnostic),log10(recmet_results.total.shape.singfilt),log10(recmet_results.total.shape.precon), ...
    log10(recmet_results.total.shape.jade),log10(recmet_results.total.shape.mscpa))

fprintf('Shape Improvement Over LS: LS: %.2f, Filt: %.2f, Precon: %.2f, JADE: %.2f , MSCPA: %.2f \n',...
    1,1/(recmet_results.total.shape.singfilt/recmet_results.total.shape.agnostic),1/(recmet_results.total.shape.precon/recmet_results.total.shape.agnostic), ...
    1/(recmet_results.total.shape.jade/recmet_results.total.shape.agnostic),1/(recmet_results.total.shape.mscpa/recmet_results.total.shape.agnostic))

%% show images
show_coloredimages = 1;

if show_coloredimages
    xvals = 150:490; 
    color1 = spectrumRGB(450);
    color2 = spectrumRGB(500);
    color3 = spectrumRGB(550);
    color4 = spectrumRGB(600);
    color5 = spectrumRGB(650);
    color6 = spectrumRGB(700);
    color_list = [color1;color2;color3;color4;color5;color6];

    % for spec=1:size(spec_all.spec.jade,2)
    run=1;
    obj_multi = [];
    agnostic_multi = [];
    singfilt_multi = [];
    jade_multi = [];
    mscpa_multi = [];
    precon_multi = [];
    for spec=1:5
        obj_multi(:,:,spec) = Sinc_results.each.ground{1}(xvals)'*vars_dict.spec_ground(spec);
        agnostic_multi(:,:,spec) = abs(Sinc_results.each.agnostic{spec}(xvals)');%*spec_all.spec.dfov(run,comp,:);
        singfilt_multi(:,:,spec) = abs(Sinc_results.each.singfilt{spec}(xvals)');%*spec_all.spec.dfov(run,comp,:);
        jade_multi(:,:,spec) = abs(Sinc_results.each.jade{spec}(xvals)');%*spec_all.spec.jade(run,comp,:);
        mscpa_multi(:,:,spec) = abs(Sinc_results.each.mscpa{spec}(xvals)');%*spec_all.spec.mscpa(run,:);
        precon_multi(:,:,spec) = abs(Sinc_results.each.precon{spec}(xvals)');
    end
    mscpa_multi(isnan(mscpa_multi)) = 0;
    
    obj_color = transform_multitocolor(obj_multi,color_list);
    agnostic_color = transform_multitocolor(agnostic_multi,color_list);
    singfilt_color = transform_multitocolor(singfilt_multi,color_list);
    jade_color = transform_multitocolor(jade_multi,color_list);
    mscpa_color = transform_multitocolor(mscpa_multi,color_list);
    precon_color = transform_multitocolor(precon_multi,color_list);
    
    figure,image(rescale(permute(obj_color,[2,1,3]))),title('Ground')
    figure,image(rescale(permute(agnostic_color,[2,1,3]))),title('Agnostic')
    figure,image(rescale(permute(singfilt_color,[2,1,3]))),title('Filt')
    figure,image(rescale(permute(jade_color,[2,1,3]))),title('JADE')
    figure,image(rescale(permute(mscpa_color,[2,1,3]))),title('MS-CPA')
    figure,image(rescale(permute(precon_color,[2,1,3]))),title('Opt Precon')


    do_save = 0;
    if do_save
        imwrite(rescale(obj_color),'obj.png')
        imwrite(rescale(agnostic_color),'ls.png')
        imwrite(rescale(singfilt_color),'filt.png')
        imwrite(rescale(jade_color),'msbss.png')
        imwrite(rescale(mscpa_color),'mscpa.png')
        imwrite(rescale(precon_color),'precon.png')
    end

end






