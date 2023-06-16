%Code for running experiments where the object is changing intensity with
%respect to the clutter which is at same intensity. 
%
%Can run whole code to get outcomes
%% add paths
addpath('..\..\');
addpath(genpath('..\..\utils'));
addpath(genpath('..\..\important_functions'))

%real measured gamma curve (convert pixel values to normalized radiance)
gamma_curve3 = @(x) -8.4e-15*x.^6 + 5.6e-12*x.^5 - 1.5e-09*x.^4 + 2.1e-7*x.^3+4.7e-7*x.^2+.00013*x-.00093;

%% File/run information
basefolder = dir('..\..\').folder;
expfolder = strcat(basefolder,"\calibration_and_data\multipleintensities_oneobject_satinsurface\"); 

%% ground truth parameters and clutter info
ground_params = struct(); %holds information about ground truth (CFOV scene)
ground_params.clutter_locations = {};
ground_params.clutter_locations{1} = strcat(expfolder,"int_10\ROI_AR_sp2_spots300.csv");
ground_params.clutter_locations{2} = strcat(expfolder,"int_10\ROI_AR_sp3_spots300.csv");
ground_params.single_ground = 0; %whether there is just one ground or allow multiple
ground_params.single_ground_location = ""; %location of single ground file

%load in clutter information
clut_spec_file{1} = csvread(strcat(expfolder,'\clut1_radiance.csv')); %size numSpecs x numInts 
clut_spec_file{2} = csvread(strcat(expfolder,'\clut2_radiance.csv')); %size numSpecs x numInts 
[clutter_rads,clutter_specs] = findSpecWithPixelFile(clut_spec_file,4);

%% Experiment Parameters

[scene_params, brdf_params, recon_params, msbss_params, mscpa_params, lsrecon_params] ...
    = get_params_multipleintensities_satin(basefolder);
recon_params.do_optprecon = 0; %whether to perform optimized preconditioning or not

%% define the different folder locations for each run
%load in object radiance values
filenames = ["ROI_AR_sp4_spots300.csv","ROI_AR_sp5_spots300.csv","ROI_AR_sp6_spots300.csv","ROI_AR_sp7_spots300.csv","ROI_AR_sp8_spots300.csv"];
offset = 4;
[lfield_locations_list,obj_spectrums] = find_filelocations(expfolder,filenames,gamma_curve3,offset);

%% debug mode (if stop is put into code)
do_debug=1; 
if do_debug
    for exp_ind = 3    

        %set specific parameters
        ground_params.ground_locations = strcat(expfolder,"int_",num2str(exp_ind),"\ROI_AR_sp1_spots300.csv");
        ground_params.spec_ground = obj_spectrums(exp_ind,:);
        if iscell(clutter_specs)
            for clutter_i = 1:length(clutter_specs)
                ground_params.spec_clut(clutter_i,:) = clutter_specs{clutter_i}(exp_ind,:);
            end
        else %single clut
            ground_params.spec_clut(1,:) = clutter_specs(exp_ind,:);
        end

        %set parameters
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

    %set specific parameters
    ground_params.ground_locations = strcat(expfolder,"int_",num2str(exp_ind),"\ROI_AR_sp1_spots300.csv");
    ground_params.spec_ground = obj_spectrums(exp_ind,:);
    if iscell(clutter_specs)
        for clutter_i = 1:length(clutter_specs)
            ground_params.spec_clut(clutter_i,:) = clutter_specs{clutter_i}(exp_ind,:);
        end
    else %single clut
        ground_params.spec_clut(1,:) = clutter_specs(exp_ind,:);
    end

    %set parameters
    all_params = struct();
    all_params.scene_params = scene_params;
    all_params.brdf_params = brdf_params;
    all_params.recon_params = recon_params;
    all_params.mscpa_params = mscpa_params;
    all_params.msbss_params = msbss_params;
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

%combine all results into single structure across all runs
[recmet_all,dist_all,predmet_all,Sinc_all,separated_all]...
    = load_across_runs(recmet_results,dist_results,predmet_results,Sinc_results,separated_results);

%normalized scattering incident distributions
Sinc_all_norm = normalize_Sinc_struct(Sinc_all);

%% show results
%make an entry "nan" if you do not want to show it

%show normalized MSE results
% recmet_all.total.shape.precon = nan;
recmet_all.total.shape.noprecon = nan;
plot_struct_fields(recmet_all.total.shape,log10(predmet_all.pow_rat),1,'Shape Error vs SCR','Log of SBR','Log Shape Error')
grid on
xlim([1.2,3.5])
ylim([-5.2,-2.5])

%show reconstructing incident scattering distributions Sinc
ints_toshow = [1,5,9];
Sinc_all_norm.total.clutter = nan;
plot_struct_fields_lincs(Sinc_all_norm.total,vars_dict{1}.inc_angles,[-64,-31],predmet_all.pow_rat,ints_toshow, ...
    'Normalized Linc','Incident Angle','Recovered Radiance')

%% show results with color
show_color = 1;
run = 5; %low run means harder
xvals = 252:600;

color1 = [0 0.4470 0.7410];
color2 = [0.8500 0.3250 0.0980];
color3 = [0.9290 0.6940 0.1250];
color4 = [0.4940 0.1840 0.5560];
color5 = [0.4660 0.6740 0.1880];
color6 = [0.3010 0.7450 0.9330];
color7 = [0.6350 0.0780 0.1840];
color_list = [color1;color2;color3;color4;color5;color6;color7];

if show_color == 1
    nspecs = 5;
    obj_multi = zeros(length(xvals),1,nspecs);
    agnostic_multi = zeros(length(xvals),1,nspecs);
    jade_multi = zeros(length(xvals),1,nspecs);
    mscpa_multi = zeros(length(xvals),1,nspecs);
    precon_multi = zeros(length(xvals),1,nspecs);
    for spec=1:5
        obj_multi(:,:,spec) = Sinc_all.each.ground(run,1,xvals)*vars_dict{run}.spec_ground(spec);
        agnostic_multi(:,:,spec) = Sinc_all.each.agnostic(run,spec,xvals);%*spec_all.spec.dfov(run,comp,:);
        jade_multi(:,:,spec) = abs(Sinc_all.each.jade(run,spec,xvals));%*spec_all.spec.jade(run,comp,:);
        mscpa_multi(:,:,spec) = abs(Sinc_all.each.mscpa(run,spec,xvals));%*spec_all.spec.msdfov(run,:);
        if recon_params.do_optprecon
            precon_multi(:,:,spec) = Sinc_all.each.precon(run,spec,xvals);
        end
    end
    %account for MSDFOV only uses 3 spectra
    % msdfov_multi(:,:,1) = msdfov_multi(:,:,2);
    % msdfov_multi(:,:,5) = msdfov_multi(:,:,4);
    
    mscpa_multi(isnan(mscpa_multi)) = 0;
    
    obj_color = transform_multitocolor(obj_multi,color_list);
    agnostic_color = transform_multitocolor(agnostic_multi,color_list);
    jade_color = transform_multitocolor(jade_multi,color_list);
    mscpa_color = transform_multitocolor(mscpa_multi,color_list);
    if recon_params.do_optprecon
        precon_color = transform_multitocolor(precon_multi,color_list);
    end
    
    figure,image(rescale(permute(obj_color,[2,1,3]))),title('Ground')
    figure,image(rescale(permute(agnostic_color,[2,1,3]))),title('Agnostic')
    figure,image(rescale(permute(jade_color,[2,1,3]))),title('JADE')
    figure,image(rescale(permute(mscpa_color,[2,1,3]))),title('MS-CPA')
    if recon_params.do_optprecon
        figure,image(rescale(permute(precon_color,[2,1,3]))),title('Opt Precon')
    end

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






