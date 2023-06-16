%return scattered angles for each point on surface
%
% Input:
% 
% xoff -- offset of position to camera
% cam_dist -- distance from sample to lens (cm)
% cam_angs -- camera angles

function[scat_angles] = calc_scat_angs_geometry(xoff,cam_dist,cam_angs)

scat_angles = [];

for ang_ind=1:length(cam_angs)
    ang = cam_angs(ang_ind);
    scat_angles(ang_ind) = atand((cam_dist*sind(ang)-xoff)/(cam_dist*cosd(ang)));
end
