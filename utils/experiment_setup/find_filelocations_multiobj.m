
function[lfield_locations_list,obj_spectrums] = find_filelocations_multiobj(basename,numObj_list,gamma_curve)

foldername = {};
obj_rad = {};
for numObj_i = 1:length(numObj_list)
    numObj = numObj_list(numObj_i);
    toskip = 1+1+numObj; %background+clutter+objects
    foldername{numObj_i} = strcat(basename,"obj_",num2str(numObj),"\");

    %load in object radiance values
    for obj=1:numObj
        obj_rad{numObj_i}(obj,:) = csvread(strcat(foldername{numObj_i},"obj",num2str(obj),"_radiance.csv")); %radiance of objects in size numSpecs x numInts
    end

    %define object spectrums
    all_obj_specs{numObj_i} = [];
    for i=1:size(obj_rad{numObj_i},1) %each obj
        for j=(toskip+1):size(obj_rad{numObj_i},2) %each spectrum value
            all_obj_specs{numObj_i}(i,j-toskip) = gamma_curve(obj_rad{numObj_i}(i,j)); %skip first value
        end
    end
    max_obj_spec_val{numObj_i} = max(all_obj_specs{numObj_i}(:)); %maximum intensity value of spec
    
    %find the lfield names
    filenames{numObj_i} = [""];
    maxSpecs = 6;
    for i=1:maxSpecs
        filenames{numObj_i}(i) = strcat("ROI_AR_sp",num2str(i+toskip-1),"_spots300.csv");
    end
end

    
%find lfield locations
lfield_locations_list = {};
obj_spectrums = [];
for numObj_i=1:length(numObj_list)
    folder = foldername{numObj_i};

    %normalize to max spec value
%     obj_spectrums{numObj_i} = all_obj_specs{numObj_i}/max(all_obj_specs{numObj_i}(:)); %spectrum

    %normalize to last value (what experiment ran)
    obj_spectrums{numObj_i} = all_obj_specs{numObj_i} ...
        ./repmat(all_obj_specs{numObj_i}(:,end),1,size(all_obj_specs{numObj_i},2)); %spectrum

    lfield_locations_list{numObj_i} = [""];
    for filename_ind =1:length(filenames{numObj_i})
        filename = filenames{numObj_i}(filename_ind);
        lfield_locations_list{numObj_i}(filename_ind) = strcat(folder,filename);
    end
end