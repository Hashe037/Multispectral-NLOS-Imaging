%lfield, lfield_ground, lfield_clutter, spec_ground, spec_clut, num_ground, num_clutter

function[meas_params] ...
    = loadinLfields(ground_params,lfield_locations,back_sig,num_specs)

%measured experimental light fields
lfield = {};
for ind=1:num_specs
    lfield{ind} = csvread(lfield_locations(ind))-back_sig;
end 

%ground truth
lfield_ground = {};
spec_ground = [];
num_ground = length(ground_params.ground_locations); %number of ground measurements
for ground_i = 1:num_ground
    lfield_ground{ground_i} = csvread(ground_params.ground_locations(ground_i))-back_sig;
    spec_ground(ground_i,:) = ground_params.spec_ground(ground_i,1:num_specs);
end

%clutter 
lfield_clutter = {};
spec_clut = [];
num_clutter = length(ground_params.clutter_locations);
for clutter_i = 1:num_clutter
    lfield_clutter{clutter_i} = csvread(ground_params.clutter_locations{clutter_i})-back_sig;
    spec_clut(clutter_i,:) = ground_params.spec_clut(clutter_i,1:num_specs);
end

meas_params = struct();
meas_params.lfield = lfield;
meas_params.lfield_ground = lfield_ground;
meas_params.lfield_clutter = lfield_clutter;
meas_params.spec_ground = spec_ground;
meas_params.spec_clut = spec_clut;

% 

