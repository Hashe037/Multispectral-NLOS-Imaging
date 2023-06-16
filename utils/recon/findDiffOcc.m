%find dalpha_occ except not in SVD but in regular angle coordinates
%Useful for using just spatial differential for preconditioner.

function[diff_occ] = findDiffOcc(brdf,vant_pos,inc_angs,vant_angs,cam_angs,brdf_scat_angs,occ_x,occ_y)

occ_x = -occ_x; %make d negative

total_min = min(brdf_scat_angs(:));
total_max = max(brdf_scat_angs(:));
ang_range = total_max-total_min;

beta_ang_cell = {};% = zeros(size(brdf,1),length(vantage_pos)-1); %angular response per vantage
for i=1:(length(vant_pos)-1) %for each vantage except last one
    x_1 = vant_pos(i);
    x_2 = vant_pos(i+1);
    theta_1 = atand((occ_x-x_1)/occ_y); 
    theta_2 = atand((occ_x-x_2)/occ_y); 
    abo = find(inc_angs<theta_1); %find all angles above angles
    belo = find(inc_angs>theta_2); %find all angles below second angles
    delta_ang = intersect(abo,belo); %find INDICES of angles that correspond
    
    
    %truncation values to just the outgoing angles allowed by certain vantage
    %position
    ang_start = round((vant_angs(i,1)-total_min)/ang_range*size(brdf,1))+1; %ratio of vantage min to brdf min
    ang_end = round((vant_angs(i,2)-total_min)/ang_range*size(brdf,1)); %ratio of vantage min to brdf min

    %save beta_ang_cell
    if length(delta_ang)==0 && i ~= 1
        beta_ang{i} = beta_ang{i-1};
    else
        beta_ang{i} = sum(brdf(ang_start:ang_end,delta_ang),2)/(length(delta_ang));% integrate over incident domain %dont divid???
%         beta_ang{i} = sum(brdf(ang_start:ang_end,delta_ang),2)'/(length(delta_ang));% integrate over incident domain %dont divid???
    end
end
    
beta_ang{length(vant_pos)} = beta_ang{length(vant_pos)-1};

%truncate down to camera angles
% [~,Locb] = ismember(cam_angles,angles);
% Locb(Locb==0) = [];
% beta_ang = beta_ang(Locb,:);

%downsample outgoing angles to the camera angles measured in experiment
%(angles of camera movement)
diff_occ = zeros(length(cam_angs),length(vant_pos));
for i=1:(length(vant_pos)) %for each vantage 
    samp_points = linspace(1,length(beta_ang{i}),length(cam_angs));
    diff_occ(:,i) = interp1(1:length(beta_ang{i}),beta_ang{i},samp_points);
%     forward_beta(:,i) = forward_beta(:,i)/norm(forward_beta(:,i))*U_vant_norms{i}(j); %norm should be same as U_vant{i}
end


% [~,Locb] = ismember(cam_angles,scattered_angles);
% Locb(Locb==0) = [];
% beta_ang = beta_ang(:,Locb);
%     