%calculate final metrics with the Lincbest for single ground source

%L_inc_total--total array weighted according to spectral content as well

%ver==1: go by spectrum
%ver==2: sum over all spectrum

function[shape_mse,mse,pow_dif,L_inc_total]  = calculate_metrics_total(L_inc_best,best_i,spec_c,H,coeff_multiple,L_inc_g,power_diff,ver)

total_spots = length(L_inc_g);

%create total array
if ver==1
    L_inc_total = sum_cell(L_inc_best,H(best_i,spec_c));
else
    L_inc_total = sum_cell(L_inc_best,sum(H(best_i,:),2));
end

mse =  nansum((-L_inc_total*coeff_multiple+L_inc_g).^2)/total_spots;
shape_mse = nansum((-L_inc_total/norm(L_inc_total)+L_inc_g/norm(L_inc_g)).^2)/total_spots;
pow_dif = power_diff(L_inc_total*coeff_multiple,L_inc_g);