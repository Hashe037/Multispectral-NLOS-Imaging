%Smooth resulting Sincs for better comparison and to remove high frequency
%noise. Smoothing is done with a smoothing spline
%
%Don't need to do MS-CPA since it was already smoothed in its code (it
%needs harsher smoothing).

function[meas_params,results_agnostic,results_precon,results_jade,results_mscpa] ...
    = SmoothSincs(smooth_param,meas_params,lsrecon_params,results_agnostic,results_precon,results_jade,results_mscpa)

num_ground = length(results_agnostic.Sinc_ground);
num_specs = length(results_agnostic.Sinc);

%ground truth
Sinc_ground = results_agnostic.Sinc_ground;
for i=1:length(Sinc_ground)
    f = fit((1:length(Sinc_ground{i}))',Sinc_ground{i}','smoothingspline','SmoothingParam',smooth_param);
    Sinc_ground{i} = f(1:length(Sinc_ground{i}))';
end
results_agnostic.Sinc_ground = Sinc_ground;

% Sinc_ground = results_agnostic.Sinc_ground;
% for i=1:length(Sinc_ground)
%     f = fit((1:length(Sinc_ground{i}))',Sinc_ground{i}','smoothingspline','SmoothingParam',smooth_param);
%     Sinc_ground{i} = f(1:length(Sinc_ground{i}))';
% end
% results_agnostic.Sinc_ground = Sinc_ground;

%regular measurements
Sinc = results_agnostic.Sinc;
for i=1:length(Sinc)
    f = fit((1:length(Sinc{i}))',Sinc{i}','smoothingspline','SmoothingParam',smooth_param);
    Sinc{i} = f(1:length(Sinc{i}))';
end
results_agnostic.Sinc = Sinc;

%single spectra best results
Sinc_singfilt = results_agnostic.Sinc_singfilt;
for i=1:length(Sinc_singfilt)
    f = fit((1:length(Sinc_singfilt{i}))',Sinc_singfilt{i}','smoothingspline','SmoothingParam',smooth_param);
    Sinc_singfilt{i} = f(1:length(Sinc_singfilt{i}))';
end
results_agnostic.Sinc_singfilt = Sinc_singfilt;

%MS-BSS with JADE
if isfield(results_jade,'Sinc_jade_best')
    Sinc_jade_best = results_jade.Sinc_jade_best;
    for i=1:length(Sinc_jade_best)
        f = fit((1:length(Sinc_jade_best{i}))',Sinc_jade_best{i}','smoothingspline','SmoothingParam',smooth_param);
        Sinc_jade_best{i} = f(1:length(Sinc_jade_best{i}))';
    end
    results_jade.Sinc_jade_best = Sinc_jade_best;
end

%if single ground
if isfield(results_agnostic,'Sinc_sing_ground')
    Sinc_sing_ground = results_agnostic.Sinc_sing_ground;
    for i=1:length(Sinc_sing_ground)
        f = fit((1:length(Sinc_sing_ground{i}))',Sinc_sing_ground{i}','smoothingspline','SmoothingParam',smooth_param);
        Sinc_sing_ground{i} = f(1:length(Sinc_sing_ground{i}))';
    end
    results_agnostic.Sinc_sing_ground = Sinc_sing_ground;
end

%if doing optimization preconditioning
if isfield(results_precon,'Sinc_precon')
    %have less smoothing for L_inc_precon since it is more sensitive to it
    if isfield(lsrecon_params,'precon_smooth')
        smooth_param_precon = lsrecon_params.precon_smooth;
    else
        smooth_param_precon = smooth_param*10; %*50
    end
    Sinc_precon = results_precon.Sinc_precon;
    Sinc_noprecon = results_precon.Sinc_noprecon;
    Sinc_precon_ground = results_precon.Sinc_precon_ground;
    for i=1:length(Sinc_precon)
        f = fit((1:length(Sinc_precon{i}))',Sinc_precon{i},'smoothingspline','SmoothingParam',smooth_param_precon);
        Sinc_precon{i} = f(1:length(Sinc_precon{i}));
    end
    for i=1:length(Sinc_noprecon)
        f = fit((1:length(Sinc_noprecon{i}))',Sinc_noprecon{i},'smoothingspline','SmoothingParam',smooth_param_precon);
        Sinc_noprecon{i} = f(1:length(Sinc_noprecon{i}));
    end
    for i=1:length(Sinc_precon_ground)
        f = fit((1:length(Sinc_precon_ground{i}))',Sinc_precon_ground{i},'smoothingspline','SmoothingParam',smooth_param_precon);
        Sinc_precon_ground{i} = f(1:length(Sinc_precon_ground{i}));
    end  

    results_precon.Sinc_precon = Sinc_precon;
    results_precon.Sinc_noprecon = Sinc_noprecon;
    results_precon.Sinc_precon_ground = Sinc_precon_ground;

    %upsample Sinc and save the smoothed versions
    Sinc_precon_up = {};
    Sinc_noprecon_up = {};
    num_vants_precon = length(Sinc_precon{1});
    num_vants = length(Sinc_ground{1});
    for spec_c = 1:length(Sinc_precon)
        Sinc_precon_up{spec_c} = interp1(Sinc_precon{spec_c}, linspace(1,num_vants_precon,num_vants));
        Sinc_noprecon_up{spec_c} = interp1(Sinc_noprecon{spec_c}, linspace(1,num_vants_precon,num_vants));
    end
    
    for ground_i = 1:length(Sinc_precon_ground)
        Sinc_precon_ground_up{ground_i} = interp1(Sinc_precon_ground{ground_i}.*sum(meas_params.spec_ground(ground_i,:)), linspace(1,num_vants_precon,num_vants));
    end
    results_precon.Sinc_precon_up = Sinc_precon_up;
    results_precon.Sinc_noprecon_up = Sinc_noprecon_up;
    results_precon.Sinc_precon_ground_up = Sinc_precon_ground_up;



end