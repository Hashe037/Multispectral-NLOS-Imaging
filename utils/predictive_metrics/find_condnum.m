%find condition number of light field over all spectra


function[cond_num,cond_num2,cond_num_norm,cond_num2_norm] ...
    = find_condnum(lfield_ms)

numSpecs = length(lfield_ms);

total_lfield = zeros(numSpecs,numel(lfield_ms{1})); %hold all light fields
total_lfield_norm = zeros(numSpecs,numel(lfield_ms{1}));
for spec = 1:numSpecs
    total_lfield(spec,:) = lfield_ms{spec}(:);
    total_lfield_norm(spec,:) = total_lfield(spec,:)/norm(total_lfield(spec,:));
end
% total_lfield_norm = (total_lfield-mean(total_lfield(:)))/sqrt(var(total_lfield(:)));
[~,Dmat,~] = svd(total_lfield,'econ');
Dmat = diag(Dmat);
cond_num = Dmat(1)/Dmat(end);
cond_num2 = Dmat(1)/Dmat(3);
[~,Dmat,~] = svd(total_lfield_norm,'econ');
Dmat = diag(Dmat);
cond_num_norm = Dmat(1)/Dmat(end);
cond_num2_norm = Dmat(1)/Dmat(3);