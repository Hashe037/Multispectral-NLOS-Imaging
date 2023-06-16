%calculate final metrics with the Sincbest for single ground source

function[shape_mse,mse,pow_dif]  = calculate_metrics_withbest_single(Sinc_best,best_i,H,coeff_multiple,Sinc_g,power_diff)

total_spots = length(Sinc_g);

mse =  nansum((-Sinc_best*sum(H(best_i,:))*coeff_multiple+Sinc_g).^2)/total_spots;
shape_mse = nansum((-Sinc_best*sum(H(best_i,:))/norm(Sinc_best*sum(H(best_i,:)))+Sinc_g/norm(Sinc_g)).^2)/total_spots;
pow_dif = power_diff(Sinc_best*sum(H(best_i,:))*coeff_multiple,Sinc_g);