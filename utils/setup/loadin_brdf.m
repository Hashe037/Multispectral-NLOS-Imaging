%loadin the BRDF
function[brdf, brdf_scat_angs, brdf_inc_angs, brdf_scat_range] ...
    = loadin_brdf(brdf_params)

load(brdf_params.brdf_location); %loads in 'brdf' size numScatAngles x numIncAngles
%alongside total_max/total_min (max/min scattering angles) and inc_min/inc_max
%(max/min incident angles)
brdf = brdf/norm(brdf); %normalize BRDF

brdf_scat_angs = linspace(total_min,total_max,size(brdf,1)); %deg, scattering angles of BRDF
brdf_inc_angs = -linspace(inc_min,inc_max,size(brdf,2)); %deg, incident angles of the BRDF
brdf_scat_range = total_max-total_min; %deg, full scattered angles range of BRDF