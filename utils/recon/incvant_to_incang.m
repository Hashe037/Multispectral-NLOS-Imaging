%convert vantage-based incident radiance to degrees by the portion of the
%occluder that is being revealed by that vantage

function[Sinc] = incvant_to_incang(Sinc_vant,inc_angs,vant_pos,occ_x,occ_y)

Sinc = zeros(1,length(inc_angs)); 

for i=1:size(Sinc_vant,2) %for each vantage except last one
    x_1 = vant_pos(i);
    x_2 = vant_pos(i+1);
    theta_1 = atand((x_1+occ_x)/occ_y);
    theta_2 = atand((x_2+occ_x)/occ_y);
     
    %for inc_angles that are negative (new)
    abo = find(inc_angs<(-1*theta_1)); %find all angles above angles
    belo = find(inc_angs>(-1*theta_2)); %find all angles below second angles
    
    delta_ang = intersect(abo,belo); %find INDICES of angles that correspond
    Sinc(delta_ang) = Sinc(delta_ang)+Sinc_vant(i);
end
