%find the non-normalized correlation between two light fields


function[nonnorm_corr] = find_nonnorm_corr(lfield1,lfield2)

nonnorm_corr = dot(lfield1(:),lfield2(:))*norm(lfield1(:))*norm(lfield2(:));