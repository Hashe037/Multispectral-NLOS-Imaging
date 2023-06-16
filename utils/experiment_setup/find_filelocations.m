

function[lfield_locations_list,obj_spectrums] ...
    = find_filelocations(basename,filenames,gamma_curve,offset)

obj_rad = csvread(strcat(basename,"obj_radiance.csv")); %radiance of objects in size numSpecs x numInts

%define object spectrums
all_obj_specs = [];
for i=1:size(obj_rad,2) %each int
    for j=(offset+1):size(obj_rad,1) %each spectrum value
        all_obj_specs(i,j-offset) = gamma_curve(obj_rad(j,i)); %skip first value
    end
end
max_obj_spec_val = max(all_obj_specs(:)); %maximum intensity value of spec

%find filenames
folder_names = [""];
for i=1:size(obj_rad,2)
    folder_names(end+1) = strcat(basename,"int_",num2str(i));
end
folder_names = folder_names(2:end);


%load in relevant data and spectrums
lfield_locations_list = {};
obj_spectrums = [];
for folder_ind=1:length(folder_names)
    folder = folder_names(folder_ind);

    %normalize to max value (intuitive)
    obj_spectrums(folder_ind,:) = all_obj_specs(folder_ind,:)/max(all_obj_specs(folder_ind,:)); %spectrum
    
    %normalize to last value (what was run in config_2)
    obj_spectrums(folder_ind,:) = all_obj_specs(folder_ind,:)...
        ./repmat(all_obj_specs(folder_ind,end),1,size(all_obj_specs,2)); %spectrum

    lfield_locations_list{folder_ind} = [""];
    for filename_ind =1:length(filenames)
        filename = filenames(filename_ind);
        lfield_locations_list{folder_ind}(filename_ind) = strcat(folder,'\',filename);
    end
end
