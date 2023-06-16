%Take absolute value Lincs to get more consistent results
function[meas_params,results_agnostic,results_precon,results_jade,results_mscpa] ...
    = absoluteSincs(meas_params,results_agnostic,results_precon,results_jade,results_mscpa)


%ground truth
Sinc_ground = results_agnostic.Sinc_ground;
for i=1:length(Sinc_ground)
    Sinc_ground{i} = abs(Sinc_ground{i});
end
results_agnostic.Sinc_ground = Sinc_ground;

%regular measurements
Sinc = results_agnostic.Sinc;
for i=1:length(Sinc)
    Sinc{i} = abs(Sinc{i});
end
results_agnostic.Sinc = Sinc;

%single spectra best results
Sinc_singfilt = results_agnostic.Sinc_singfilt;
Sinc_singfilt = abs(Sinc_singfilt);
results_agnostic.Sinc_singfilt = Sinc_singfilt;

%MS-BSS with JADE
if isfield(results_jade,'Sinc_jade_best')
    Sinc_jade_best = results_jade.Sinc_jade_best;
    for i=1:length(Sinc_jade_best)
        Sinc_jade_best{i} = abs(Sinc_jade_best{i});
    end
    results_jade.Sinc_jade_best = Sinc_jade_best;
end

%if single ground
if isfield(results_agnostic,'Sinc_sing_ground')
    Sinc_sing_ground = results_agnostic.Sinc_sing_ground;
    for i=1:length(Sinc_sing_ground)
        Sinc_sing_ground{i} = abs(Sinc_sing_ground{i});
    end
    results_agnostic.Sinc_sing_ground = Sinc_sing_ground;
end

%if doing optimization preconditioning
if isfield(results_precon,'Sinc_precon')
    %have less smoothing for L_inc_precon since it is more sensitive to it
    Sinc_precon = results_precon.Sinc_precon;
    Sinc_noprecon = results_precon.Sinc_noprecon;
    Sinc_precon_ground = results_precon.Sinc_precon_ground;
    Sinc_precon{i} = abs(Sinc_precon{i});
    Sinc_noprecon{i} = abs(Sinc_noprecon{i});
    Sinc_precon_ground{i} = abs(Sinc_precon_ground{i});
    results_precon.Sinc_precon = Sinc_precon;
    results_precon.Sinc_noprecon = Sinc_noprecon;
    results_precon.Sinc_precon_ground = Sinc_precon_ground;
end
