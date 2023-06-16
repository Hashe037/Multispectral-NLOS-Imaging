%Code for running experiment with single real-life experiment 
%multiple objects

%% add paths
addpath('..\..\');
addpath(genpath('..\..\utils'));
addpath(genpath('..\..\important_functions'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File/run information

%real
basefolder = dir('..\..\').folder;
expfolder = strcat(basefolder,"\calibration_and_data\real_multiobj_satinsurface\"); 
expname = "all_obj_run\";

objname{1} = "blue_obj";
objname{2} = "green_obj";
objname{3} = "red_obj";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ground information and loading in spectral content
ground_params = struct(); %holds information about ground truth
ground_params.ground_locations(1) = strcat(expfolder,objname{1},"\light_fields\lfield_3.csv");
ground_params.ground_locations(2) = strcat(expfolder,objname{2},"\light_fields\lfield_3.csv");
ground_params.ground_locations(3) = strcat(expfolder,objname{3},"\light_fields\lfield_3.csv");

load(strcat(expfolder,objname{1},"\obj_spec.mat"));
obj_spec_temp(1,:) = obj_spec;
load(strcat(expfolder,objname{2},"\obj_spec.mat"));
obj_spec_temp(2,:) = obj_spec;
load(strcat(expfolder,objname{3},"\obj_spec.mat"));
obj_spec_temp(3,:) = obj_spec;
obj_spec = obj_spec_temp;
ground_params.spec_ground = obj_spec;

%combining all ground into one
ground_params.single_ground = 1;
ground_params.single_ground_location = strcat(expfolder,'all_obj_ground');

%% clutter info
ground_params.clutter_locations = {};
ground_params.clutter_locations{1} = strcat(expfolder,"hidden_side_clutter_1\light_fields\lfield_2.csv");
ground_params.clutter_locations{2} = strcat(expfolder,"hidden_side_clutter_2\light_fields\lfield_2.csv");
ground_params.clutter_locations{3} = strcat(expfolder,"observer_side_clutter_1\light_fields\lfield_2.csv");
load(strcat(expfolder,expname,"clut_spec.mat"))
ground_params.spec_clut = clut_spec;

%% Experiment Parameters
[scene_params, brdf_params, recon_params, msbss_params, mscpa_params, lsrecon_params] ...
    = get_params_real_multiobj_satin(basefolder);

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

    xvals = 150:490; %opt5
    color1 = spectrumRGB(450);
    color2 = spectrumRGB(500);
    color3 = spectrumRGB(550);
    color4 = spectrumRGB(600);
    color5 = spectrumRGB(650);
    color6 = spectrumRGB(700);
    color_list = [color1;color2;color3;color4;color5;color6];

    run=1;
    obj_multi = [];
    agnostic_multi = [];
    filt_multi = [];
    jade_multi = [];
    mscpa_multi = [];
    precon_multi = [];
    for spec=1:5
        obj_multi(:,:,spec) = Sinc_results.each.groundone{spec}(xvals)';
        agnostic_multi(:,:,spec) = abs(Sinc_results.each.agnostic{1,spec}(xvals)');%*spec_all.spec.agnostic(run,comp,:);
        filt_multi(:,:,spec) = abs(Sinc_results.each.singfilt{1,spec}(xvals)');%*spec_all.spec.agnostic(run,comp,:);
        jade_multi(:,:,spec) = abs(Sinc_results.each.jade{1,spec}(xvals)');%*spec_all.spec.jade(run,comp,:);
        mscpa_multi(:,:,spec) = abs(Sinc_results.each.mscpa{1,spec}(xvals)');%*spec_all.spec.mscpa(run,:);
        precon_multi(:,:,spec) = abs(Sinc_results.each.precon{1,spec}(xvals)');

        for obji=2:3
            mscpa_multi(:,:,spec) = mscpa_multi(:,:,spec)+abs(Sinc_results.each.mscpa{obji,spec}(xvals)');
            jade_multi(:,:,spec) = jade_multi(:,:,spec)+abs(Sinc_results.each.jade{obji,spec}(xvals)');%*spec_all.spec.jade(run,comp,:);
        end
        
    end
    
    obj_color = transform_multitocolor(obj_multi,color_list);
    agnostic_color = transform_multitocolor(agnostic_multi,color_list);
    filt_color = transform_multitocolor(filt_multi,color_list);
    jade_color = transform_multitocolor(jade_multi,color_list);
    mscpa_color = transform_multitocolor(mscpa_multi,color_list);
    precon_color = transform_multitocolor(precon_multi,color_list);

    figure,image(rescale(permute(obj_color,[2,1,3]))),title('Ground')
    figure,image(rescale(permute(agnostic_color,[2,1,3]))),title('LS')
    figure,image(rescale(permute(filt_color,[2,1,3]))),title('Filt')
    figure,image(rescale(permute(jade_color,[2,1,3]))),title('JADE')
    figure,image(rescale(permute(mscpa_color,[2,1,3]))),title('MS-CPA')
    figure,image(rescale(permute(precon_color,[2,1,3]))),title('Opt Precon')

end

do_save = 0;
if do_save
    imwrite(rescale(obj_color),'obj.png')
    imwrite(rescale(agnostic_color),'ls.png')
    imwrite(rescale(filt_color),'filt.png')
    imwrite(rescale(jade_color),'msbss.png')
    imwrite(rescale(mscpa_color),'mscpa.png')
    imwrite(rescale(precon_color),'precon.png')
end



