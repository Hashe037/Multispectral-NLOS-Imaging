%Code for running experiments where the number of objects is changing but
%not the intensity

%% add paths
addpath('..\..\');
addpath(genpath('..\..\utils'));
addpath(genpath('..\..\important_functions'))

%real measured gamma curve (convert pixel values to normalized radiance)
gamma_curve3 = @(x) -8.4e-15*x.^6 + 5.6e-12*x.^5 - 1.5e-09*x.^4 + 2.1e-7*x.^3+4.7e-7*x.^2+.00013*x-.00093;

%% File/run information
basefolder = dir('..\..\').folder;
expfolder = strcat(basefolder,"\calibration_and_data\singleintensity_multipleobjects_satinsurface\"); 
numObj_list = [1,2,3,4,5];

%% Ground information
ground_params = struct(); %holds information about ground truth
ground_params.clutter_locations = {};
ground_params.clutter_locations{1} = strcat(expfolder,"obj_1\ROI_AR_sp5_spots300.csv");
ground_params.single_ground = 0;
ground_params.single_ground_location = "";

clut_spec_file{1} = csvread(strcat(expfolder,'obj_1\clut1_radiance.csv')); %size numSpecs x numInts 
[clutter_rads,clutter_specs] = findSpecWithPixelFile(clut_spec_file,4);

%% Experiment Parameters
[scene_params, brdf_params, recon_params, msbss_params, mscpa_params, lsrecon_params] ...
    = get_params_multipleobjects_satin(basefolder);

recon_params.do_optprecon = 0;

%% Define the different folder locations for each run
[lfield_locations_list,obj_spectrums] = find_filelocations_multiobj(expfolder,numObj_list,gamma_curve3);

%% debug mode (if stop is put into code)
do_debug=1; 
if do_debug
    for exp_ind = 3
        %some extra info
        numObj = numObj_list(exp_ind);

        ground_params.ground_locations = [""];
        ground_params.spec_ground = [];
        ground_params.spec_clut = [];
        for i=1:numObj
            ground_params.ground_locations(i) = strcat(expfolder,"obj_",num2str(numObj),"\ROI_AR_sp",num2str(i),"_spots300.csv");
        end
        ground_params.spec_ground = obj_spectrums{exp_ind}(:,1:scene_params.num_specs);
        msbss_params.numcomp = numObj+1; %number of components for BSS

        if iscell(clutter_specs)
            for clutter_i = 1:length(clutter_specs)
                ground_params.spec_clut(clutter_i,:) = clutter_specs{clutter_i};
            end
        else %single clut
            ground_params.spec_clut(1,:) = clutter_specs;
        end

        %create color dictionary
        color_dict = {};
        for obj_i = 1:size(ground_params.spec_ground,1)
            color_dict{obj_i} = ground_params.spec_ground(obj_i,:);
        end
        mscpa_params.color_dict = color_dict;

        all_params = struct();
        all_params.scene_params = scene_params;
        all_params.brdf_params = brdf_params;
        all_params.recon_params = recon_params;
        all_params.mscpa_params = mscpa_params;
        all_params.msbss_params = msbss_params;
        all_params.lsrecon_params = lsrecon_params;
        all_params.ground_params = ground_params;

        %run
        [recmet_results{exp_ind},dist_results{exp_ind},predmet_results{exp_ind}, ...
            Sinc_results{exp_ind},separated_results{exp_ind},vars_dict{exp_ind}] ...
            = performMultispectralNLOSImaging(all_params, lfield_locations_list{exp_ind});
        fprintf('Debug done\n')
    end
end

%% run over all files
recmet_results = {}; %reconstruction metrics
dist_results = {}; %distances in dfovs
predmet_results = {}; %preresultsive metrics 
spec_results = {}; %specular metrics/info
ldif_results = {}; %difference in light fields
Sinc_results = {}; %reconstructions
separated_results = {}; %separated light fields and dalphas
vars_dict = {}; %important variables

all_t = 0; %total time
for exp_ind = 1:length(lfield_locations_list)
    %some extra info
    numObj = numObj_list(exp_ind);

    ground_params.ground_locations = [""];
    ground_params.spec_ground = [];
    ground_params.spec_clut = [];
    for i=1:numObj
        ground_params.ground_locations(i) = strcat(expfolder,"obj_",num2str(numObj),"\ROI_AR_sp",num2str(i),"_spots300.csv");
    end
    ground_params.spec_ground = obj_spectrums{exp_ind}(:,1:scene_params.num_specs);
%     msbss_params.numcomp = numObj; %number of components for BSS
    msbss_params.numcomp = numObj+1; %number of components for BSS
%     msbss_params.numcomp = min(numObj+2,scene_params.numSpecs); %number of components for BSS
    if iscell(clutter_specs)
        for clutter_i = 1:length(clutter_specs)
            ground_params.spec_clut(clutter_i,:) = clutter_specs{clutter_i};
        end
    else %single clut
        ground_params.spec_clut(1,:) = clutter_specs;
    end

    %create color dictionary
    color_dict = {};
    for obj_i = 1:size(ground_params.spec_ground,1)
        color_dict{obj_i} = ground_params.spec_ground(obj_i,:);
    end
    mscpa_params.color_dict = color_dict;

    %make all_params structure
    all_params = struct();
    all_params.scene_params = scene_params;
    all_params.brdf_params = brdf_params;
    all_params.recon_params = recon_params;
    all_params.msbss_params = msbss_params;
    all_params.mscpa_params = mscpa_params;
    all_params.lsrecon_params = lsrecon_params;
    all_params.ground_params = ground_params;

    %run
    tic
    [recmet_results{exp_ind},dist_results{exp_ind},predmet_results{exp_ind}, ...
        Sinc_results{exp_ind},separated_results{exp_ind},vars_dict{exp_ind}] ...
        = performMultispectralNLOSImaging(all_params, lfield_locations_list{exp_ind});
    t=toc;
    all_t = all_t+t;
    fprintf('Run %i done with time elapsed %.2f secs \n',exp_ind,t)
end
fprintf('All runs complete with total time taken of %.2f seconds \n',all_t)
beep

%% load data across all the runs
Sinc_results2 = Sinc_results;
recmet_results2 = recmet_results; 
predmet_results2 = predmet_results;
dist_results2 = nan;
separated_results2 = nan;
for i=1:length(Sinc_results2)
    Sinc_results2{i} = rmfield(Sinc_results2{i},'each');
    Sinc_results2{i} = rmfield(Sinc_results2{i},'spec');
    recmet_results2{i} = rmfield(recmet_results2{i},'each');
    predmet_results2{i} =  rmfield(predmet_results2{i},'specsim_mat');
end

[recmet_all,dist_all,predmet_all,Sinc_all,separated_all]...
    = load_across_runs(recmet_results2,dist_results2,predmet_results2,Sinc_results2,separated_results2);

Sinc_all_norm = normalize_Sinc_struct(Sinc_all);

%% show results
%make an entry "nan" if you do not want to show it

recmet_all.total.shape.singfilt = nan;
recmet_all.total.shape.noprecon = nan;
recmet_all.total.shape.precon = nan;
plot_struct_fields(recmet_all.total.shape,numObj_list,1,'Shape Error vs Number of Objects','Number of Objects','Log Shape Error')
grid on
xlim([0,6])
ylim([-5.4,-2])


%Lincs
ints_toshow = [1,2,3];
Sinc_all_norm.total.noprecon = nan;
% L_inc_all_norm.total.precon = nan;
plot_struct_fields_lincs(Sinc_all_norm.total,vars_dict{1}.inc_angles,[-64,-31],numObj_list,ints_toshow, ...
    'Normalized Linc','Incident Angle','Recovered Radiance')

%% show results with color
show_color = 1;
run = 3; %number of objects
xvals = 302:587;
% xvals = 50:336;

color1 = [0 0.4470 0.7410];
color2 = [0.8500 0.3250 0.0980];
color3 = [0.9290 0.6940 0.1250];
color4 = [0.4940 0.1840 0.5560];
color5 = [0.4660 0.6740 0.1880];
color6 = [0.3010 0.7450 0.9330];
color7 = [0.6350 0.0780 0.1840];
color_list = [color1;color2;color3;color4;color5;color6;color7];

if show_color == 1
    nspecs = 6;
    obj_multi = zeros(length(xvals),1,nspecs);
    agnostic_multi = zeros(length(xvals),1,nspecs);
    singfilt_multi = zeros(length(xvals),1,nspecs);
    jade_multi = zeros(length(xvals),1,nspecs);
    precon_multi = zeros(length(xvals),1,nspecs);
    mscpa_multi = zeros(length(xvals),1,nspecs);
    for spec=1:nspecs
        for obj_i = 1:run
            obj_multi(:,:,spec) = obj_multi(:,:,spec)+Sinc_results{run}.each.ground{obj_i}(xvals)'*vars_dict{run}.spec_ground(obj_i,spec);
            singfilt_multi(:,:,spec) = abs(Sinc_results{run}.each.singfilt{obj_i,spec}(xvals)');%*spec_all.spec.dfov(run,comp,:);
            agnostic_multi(:,:,spec) = agnostic_multi(:,:,spec)+(Sinc_results{run}.each.agnostic{obj_i,spec}(xvals)');
            jade_multi(:,:,spec) = jade_multi(:,:,spec)+Sinc_results{run}.each.jade{obj_i,spec}(xvals)';
            mscpa_multi(:,:,spec) = mscpa_multi(:,:,spec)+Sinc_results{run}.each.mscpa{obj_i,spec}(xvals)';
            if recon_params.do_optprecon
                precon_multi(:,:,spec) = abs(Sinc_results{run}.each.precon{obj_i,spec}(xvals)');
            end
        end
    end
    
    mscpa_multi(isnan(mscpa_multi)) = 0;
    
    obj_color = transform_multitocolor(obj_multi,color_list);
    agnostic_color = transform_multitocolor(agnostic_multi,color_list);
    singfilt_color = transform_multitocolor(singfilt_multi,color_list);
    jade_color = transform_multitocolor(jade_multi,color_list);
    mscpa_color = transform_multitocolor(mscpa_multi,color_list);
    if recon_params.do_optprecon
        precon_color = transform_multitocolor(precon_multi,color_list);
    end
    
    figure,image(rescale(permute(obj_color,[2,1,3]))),title('Ground')
    figure,image(rescale(permute(agnostic_color,[2,1,3]))),title('Agnostic')
    figure,image(rescale(permute(singfilt_color,[2,1,3]))),title('Sing Spec')
    if recon_params.do_optprecon
        figure,image(rescale(permute(precon_color,[2,1,3]))),title('Opt Precon')
    end
    figure,image(rescale(permute(jade_color,[2,1,3]))),title('MS-BSS')
    figure,image(rescale(permute(mscpa_color,[2,1,3]))),title('MS-CPA')

    do_save = 0;
    if do_save
        imwrite(rescale(obj_color),'obj.png')
        imwrite(rescale(agnostic_color),'ls.png')
        imwrite(rescale(singfilt_color),'filt.png')
        imwrite(rescale(jade_color),'msbss.png')
        imwrite(rescale(mscpa_color),'mscpa.png')
        if recon_params.do_optprecon
            imwrite(rescale(precon_color),'precon.png')
        end
    end
end




