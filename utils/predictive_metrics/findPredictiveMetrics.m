%Clculate the predictive metrics to estimate how hard the MS-BSS/MS-CPA
%problem will be and how hard the spectrally-agnostic version is.
%
%Note this is very hard to predict, since not only do many of the
%algorithms use different criteria for unmixing (i.e. kurtosis for JADE,
%orthogonality for PCA) but also this is not well-studied for BSS methods.
% 
%--------------------------------------------------------------------------
% Noteable Outputs
% nonnorm_corr_lfield/dalpha -- non-normalized correlation between clutter
% and ground
% cond_num/2/norm -- condition number of matrices, where the "2" signifies
% to the third element rather than to the last. Norm is that each observation is
% normalized first
% pow_rat -- ratio of clutter to object power in multispectral light fields
% spec_sim_mat -- matrix of specular similarities of different elements
%--------------------------------------------------------------------------

function[predmet_dict] = ...
    findPredictiveMetrics(meas_params)

lfield = meas_params.lfield;
lfield_ground = meas_params.lfield_ground;
lfield_clutter = meas_params.lfield_clutter;
dalpha_ground = meas_params.dalpha_ground;
dalpha_clutter = meas_params.dalpha_clutter;
spec_ground = meas_params.spec_ground;
spec_clut = meas_params.spec_clut;

num_ground = length(lfield_ground);
num_clutter = length(lfield_clutter);

%non-normalized correlation of ground to clutter
nonnorm_corr_lfield = find_nonnorm_corr(sum_cell(lfield_clutter),sum_cell(lfield_ground));
nonnorm_corr_dalpha = find_nonnorm_corr(sum_cell(dalpha_clutter),sum_cell(dalpha_ground));

%find condition number of light field across all spectra
[cond_num,cond_num2,cond_num_norm,cond_num2_norm] = find_condnum(lfield);

%find power ratio of object and clutter (USED IN THE PAPER)
obj_power = [];
for ground_i=1:num_ground
    obj_power(ground_i) = sum(norm(lfield_ground{ground_i}(:)).*spec_ground(ground_i,:));
end
clut_power = [];
for clutter_i=1:num_clutter
    clut_power(clutter_i) = sum(norm(lfield_clutter{clutter_i}(:)).*spec_clut(clutter_i,:));
end
pow_rat = sum(clut_power)/sum(obj_power);

%specular similarity between objects and clutters (ALSO USED IN THE PAPER)
specsim = @(x,y) dot(abs(x)/norm(x(:)),abs(y)/norm(y(:)));
specsim_mat = zeros(length(lfield_ground),length(lfield_clutter));
for ground_i=1:num_ground
    for clutter_i=1:num_clutter
        specsim_mat(ground_i,clutter_i) = specsim(spec_ground(ground_i,:),spec_clut(clutter_i,:)); 
    end
end
    
%predict SNR (with the given light field camera)
%shot noise limited system (signal/2 is photon count), 400*3 spatial
%averaging
%surface inhomogenities do affect results but this should be estimate
%assuming Gaussian, averaging should affect by sqrt
snr = sqrt(mean(mean(sum_cell(lfield)))/length(lfield)/2)*sqrt(400*3);
snr_db = 20*log10(snr);

% construct actual dictionary
predmet_dict = struct();
predmet_dict.pow_rat = pow_rat;
predmet_dict.specsim_mat = specsim_mat;
predmet_dict.nonnorm_corr_lfield = nonnorm_corr_lfield;
predmet_dict.nonnorm_corr_dalpha = nonnorm_corr_dalpha;
predmet_dict.cond_num = cond_num;
predmet_dict.cond_num2 = cond_num2;
predmet_dict.cond_num_norm = cond_num_norm;
predmet_dict.cond_num2_norm = cond_num2_norm;
predmet_dict.snr_db = snr_db;

% fprintf('Predictive: Inten = %.1f, Pow_rat = %.1f, Lfield_corr = %.4f, cond_num2 = %.1f \n',1/inten_coeff,pow_rat,nonnorm_corr_lfield,cond_num2)
