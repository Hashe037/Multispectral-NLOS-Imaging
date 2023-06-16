%return incident angles for each slit given lcd position WRT CENTER OF
%SAMPLE
%
% Input:
% 
% diam -- diameter of LCD screen
% phi -- 
% phi
% lcd_x -- xpos of center of LCD screen wrt wall (negative)
% lcd_y -- ypos of center of LCD wrt sample (negative)
% d
% numSlits -- number of slits

%First slit is going to be the slit most normal to the scattering surface
%while the last is the most glancing

function[inc_angles] = calc_inc_angles_geometry(diam,phi,psi,lcd_x,lcd_y,d,numSlits)

inc_angles = [];

% inc = diam/numSlits; %increment on lcd screen in distance
% x_inc = cosd(phi)*inc; %increment in x position
% y_inc = sind(phi)*inc; %increment in y position
rad = diam/2;

%calculate x and y positions for each slit
lcd_x_pos = linspace(lcd_x+rad*sind(phi)-d,lcd_x-rad*sind(phi)-d,numSlits);
lcd_y_pos = linspace(lcd_y-rad*cosd(phi),lcd_y+rad*cosd(phi),numSlits);
inc_angles = asind(lcd_x_pos./sqrt(lcd_x_pos.^2+lcd_y_pos.^2));
% inc_angles = atand(lcd_x_pos./lcd_y_pos);



