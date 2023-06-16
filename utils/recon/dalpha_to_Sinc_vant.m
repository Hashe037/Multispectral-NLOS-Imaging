%Project each vantage measurement on the occluder constraint (dalpha_occ)
%to get the expected incident radiance as a function of vantage index

function[Sinc_vant,dist] = dalpha_to_Sinc_vant(dalpha,dalpha_occ)


Sinc_vant = zeros(1,size(dalpha,2)); %number of vantage points
dist = zeros(1,size(dalpha,2)-1); %number of vantage points
% dalpha_recon_norm = zeros(1,size(dalpha_meas,2)-1); %for norm
dalpha_recon_norm = zeros(size(dalpha)); %for actual values
for i=1:size(dalpha,1) %each vantage point
    %negative due to direction of gradient in caclulating dalpha
    Sinc_vant(i) = -1*dalpha(i,:)*dalpha_occ(i,:)'/(norm(dalpha_occ(i,:))^2); 
    dist(i) = norm(dalpha(i,:)+Sinc_vant(i)*dalpha_occ(i,:)); %euclidean distance in projection
    % dalpha_recon_norm(:,i) = Sinc_vant(i)*dalpha_occ(:,i); 
end
    