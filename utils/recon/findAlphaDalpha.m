%Calculate the alpha (projected SVD vector coefficients) and dalpha (change
%in alpha over vantages) for a given multispectral light field "lfield".
%This is the basis for the hidden scene reconstruction using the light
%field, where the projection is used to limit noise and regularize the
%solution.
%
%--------------------------------------------------------------------------
% Inputs
% lfield_general -- multispectral light field (cell with length num_specs,
% each cell has light field which is size num_vants x num_angs)
% U_vant_t -- truncated left SVD of BRDF scattering angles x vector_num 
% D_tot -- BRDF SVD values
% dalpha_smooth -- smoothing constant of dalpha for the smoothing spline
%
%--------------------------------------------------------------------------
% Outputs
% alpha_meas -- projected SVD vector coefficients onto each spectra of the
% multispectral light field
% dalpha_meas -- change in alpha over vantages
% non_smooth(alpha/dalpha) -- versions that are not smoothed

function[alpha_meas,dalpha_meas,non_smooth_alpha_meas,non_smooth_dalpha_meas] = findAlphaDalpha(lfield_general,U_vant_t,D_tot,alpha_smooth)

num_specs = length(lfield_general); %number of spectral measurements
num_vants = length(U_vant_t); %number of vantages
vector_num = size(U_vant_t{1},2); %number of BRDF components

alpha_meas = {}; %alpha from measurements
dalpha_meas = {}; %change in alpha from measurements
non_smooth_alpha_meas = {}; %non-smoothed version
non_smooth_dalpha_meas = {}; %non-smoothed version

for ii = 1:num_specs
    
    %calculate alpha
    for i=1:num_vants
        alpha_meas{ii}(i,:) = project_lfield_onto_vector(U_vant_t{i},D_tot,lfield_general{ii},i);
    end
    non_smooth_alpha_meas = alpha_meas;
    
    %calculate non-smoothed dalpha
    dalpha_meast = ([alpha_meas{ii}',zeros(vector_num,1)]-[zeros(vector_num,1),alpha_meas{ii}'])';
    non_smooth_dalpha_meas{ii} = dalpha_meast(2:end-1,:);
    
    %smooth alpha
    for i=1:vector_num
%         alpha_meas{ii}(:,i) = smoothdata(alpha_meas{ii}(:,i),'gaussian',40);
%         alpha_meas{ii}(:,i) = smoothdata(alpha_meas{ii}(:,i),'movmean',40);
%         alpha_meas{ii}(:,i) = smoothdata(alpha_meas{ii}(:,i),'lowess',80);
        f = fit((1:size(alpha_meas{ii},1))',alpha_meas{ii}(:,i),'smoothingspline','SmoothingParam',alpha_smooth);
        alpha_meas{ii}(:,i) = f(1:length(alpha_meas{ii}));
%         f = fit((1:size(alpha_meas{ii},1))',alpha_meas{ii}(:,i),'cubicinterp');
%         alpha_meas{ii}(:,i) = f(1:length(alpha_meas{ii}));
    end
    
    %calculate dalpha
    dalpha_meast = ([alpha_meas{ii}',zeros(vector_num,1)]-[zeros(vector_num,1),alpha_meas{ii}'])';
    dalpha_meas{ii} = dalpha_meast(2:end-1,:);
end
