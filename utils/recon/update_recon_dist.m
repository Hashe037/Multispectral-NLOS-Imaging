%change distance metric for lfield recon (if required)

function[dist] = update_recon_dist(dist,dis_method,Sinc_vant,dalpha,dalpha_occ)

for i=1:length(dist) %each vantage
    if strcmp(dis_method,"residnorm") %norm of residual
        dist(i) = dist(i); %z star p2
    elseif strcmp(dis_method,"residnorm_totweighted") %norm of residual weighted by length of both
        dist(i) = dist(i)/(norm(Sinc_vant(i)*dalpha_occ(i,:))+norm(dalpha(i,:)))*norm(dalpha(i,:))^2;
    elseif strcmp(dis_method,"residnorm_measweighted") %norm of residual weighted just by measured dalpha
        dist(i) = dist(i)/(norm(dalpha(i,:))+.00001);
    elseif strcmp(dis_method,"residnorm_reconweighted") %norm of residual weighted by reconstruction values
        dist(i) = atan(dist(i)/(norm(Sinc_vant(i)*dalpha_occ(i,:))+.1));
    else
        error('Error: dis_method not recognized')
    end
end