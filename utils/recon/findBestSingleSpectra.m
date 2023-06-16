%Find single spectral measurement that lends the smallest reconstruction
%residual

function[Sinc_vant_singfilt,Sinc_singfilt,bestspec] ...
    = findBestSingleSpectra(Sinc,Sinc_vant,dist_recon)

dist_recon_list = [];
for i=1:length(dist_recon)
    dist_recon_list(i) = sum(dist_recon{i});
end
[~,bestspec] = min(dist_recon_list); %minimum residual

%make same size as others
Sinc_singfilt = {};
Sinc_vant_singfilt = {};
for i=1:length(Sinc)
    Sinc_singfilt{i} = zeros(size(Sinc{i}));
    Sinc_vant_singfilt{i} = zeros(size(Sinc_vant{i}));
end

Sinc_singfilt{bestspec} = Sinc{bestspec};
Sinc_vant_singfilt{bestspec} = Sinc_vant{bestspec};


