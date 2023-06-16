%find radiance of objects given the pixel intensity file (already loaded in)
%each column is the spectrum for a given run/intensity

function[rad,spec] = findSpecWithPixelFile(spec_file,toskip)

%gamma curve of oled monitor
gamma_curve3 = @(x) -8.4e-15*x.^6 + 5.6e-12*x.^5 - 1.5e-09*x.^4 + 2.1e-7*x.^3+4.7e-7*x.^2+.00013*x-.00093;

%for simulation
% gamma_curve3 =  @(x) (x/255).^(2.3);

if iscell(spec_file) %multiple files
    for file=1:length(spec_file)
        %define clutter spectrums
        rad{file} = [];
        for i=1:size(spec_file{file},2) %each int
            for j=(toskip+1):size(spec_file{file},1) %each spectrum value
                rad{file}(i,j-toskip) = gamma_curve3(spec_file{file}(j,i)); %skip first value
            end
        end
%         max_rad{file} = max(rad{file}(:)); %maximum intensity value of spec
%         spec{file} = rad{file}/max_rad{file};

        last_rad{file} = rad{file}(:,end); %last value (depends on experiment)
        spec{file} = rad{file}./repmat(last_rad{file},1,size(rad{file},2));
    end
else
    %define clutter spectrums
    rad = [];
    for i=1:size(spec_file,2) %each int
        for j=(toskip+1):size(spec_file,1) %each spectrum value
            rad(i,j-toskip) = gamma_curve3(spec_file(j,i)); %skip first value
        end
    end
%     max_rad = max(rad(:)); %maximum intensity value of spec
%     spec = rad/max_rad;

    last_rad = rad(:,end); %last value (depends on experiment)
    spec = rad./repmat(last_rad,1,size(rad,2));
end